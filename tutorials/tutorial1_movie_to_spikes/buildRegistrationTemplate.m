function [template, info] = buildRegistrationTemplate(nwbFile, options)
% BUILDREGISTRATIONTEMPLATE  Rebuild the registration template the pipeline used.
%
% [template, info] = BUILDREGISTRATIONTEMPLATE(nwbFile, options)
%
% The template every frame is registered onto is not in the deposited asset,
% because the source pipeline never saved it. It is, however, deterministically
% reproducible from two things the asset does carry: the raw movie and the
% running speed. This rebuilds it.
%
% The recipe, from large_voltage_FOV_MC_movieSegmentation_batch_v7 lines 26-38
% -----------------------------------------------------------------------------
%   1. Smooth running speed with a 1000-sample moving mean.
%   2. Call a frame "still" when that smoothed speed is below 1e-3.
%   3. Take the FIRST run of 1000 consecutive still frames.
%   4. The template is the pixelwise median of the raw frames over that run.
%
% Steps 3 and 4 are kept. Steps 1 and 2 are NOT the default here; see below.
%
% Why the stillness test was changed
% ----------------------------------
% The published test does not discriminate on this session. Averaging speed
% over 1000 samples flattens even the running bouts below 1e-3, so ALL 1,152,000
% frames qualify as still and "the first run of 1000 consecutive still frames"
% resolves to frames 1-1000 -- what no search at all would return. Those are the
% least photobleached frames in the session and they sit at a chunk start, where
% the session notes record a recording artifact.
%
% The default rule is instead the author's: a frame is still when the encoder
% registered no movement at all,
%
%     abs(loc) < min(abs(loc(loc ~= 0)))
%
% i.e. below the smallest nonzero speed the encoder can report. On this session
% the smallest nonzero |loc| is 6.26e-09, so this selects exactly the frames
% where loc is zero: 95.11% of them. It is written as a magnitude comparison
% rather than `loc == 0` so that it adapts to a real quantisation floor on other
% data, and is robust to floating-point noise around zero.
%
% Note the magnitude. `loc` here runs from -0.0229 to +0.0190 -- 2.92% of samples
% are negative, the treadmill turning backwards -- so `min(loc(loc~=0))` would be
% the most negative value, and comparing against it selects nothing at all. The
% published test has the mirror-image quirk: it compares the SIGNED smoothed
% speed against 1e-3, so every backward sample passes automatically.
%
% Effect on this session: the template moves from frames 1-1000 to frames
% 9678-10677, off the chunk start and onto genuinely motionless frames.
%
% This CHANGES the template relative to the published run, and therefore the
% corrected movie and everything below it. Pass stillRule = 'smoothedThreshold'
% to reproduce the published behaviour exactly.
%
% Why a still period, and why the median
% --------------------------------------
% Registering to a moving frame bakes that frame's own displacement into every
% correction, so the reference has to come from a period when the animal was
% still. The median rather than the mean because 1000 frames at 787 Hz spans 1.3
% seconds, in which these cells spike: the mean would smear action potentials
% into the reference image, the median rejects them as the brief outliers they
% are.
%
% What this is NOT
% ----------------
% `/acquisition/images/segmentation_reference_image` in the asset is a different
% image: the mean of a 10-second excerpt, produced by the Stage 0 segmentation
% script for its regionprops call. Scripts 1 to 5 load it and never use it. Do
% not register to it.
%
% Inputs
% ------
% nwbFile : char
%     The raw-movie excerpt asset.
%
% options : struct, all optional
%     .stillRule         'zeroMotion' (default) or 'smoothedThreshold' (the
%                        published behaviour). See above.
%     .templateFrames    run length required, samples. Default 1000.
%     .speedThreshold    'smoothedThreshold' only: cut on the smoothed signed
%                        speed. Default 1e-3.
%     .smoothingWindow   'smoothedThreshold' only: moving-mean window on speed,
%                        samples. Default 1000.
%     .searchLimit       last sample of the search. Default 288000, the first
%                        acquisition chunk, matching the source, which built the
%                        template from the first acquisition folder only. Pass
%                        Inf to search the whole deposited window.
%     .verbose           default true.
%
% Outputs
% -------
% template : [nRow x nCol] single, in raw camera counts
% info     : struct with .startSample, .stopSample, .nFrames, .meanSpeed,
%            .maxSpeed, .fractionStill, .searchLimit, .options
%
% MATLAB R2019b compatible. Reads through MATLAB's built-in HDF5 support; MatNWB
% is not required.

if nargin < 2, options = struct(); end
options = setDefault(options, 'stillRule',       'zeroMotion');
options = setDefault(options, 'templateFrames',  1000);
options = setDefault(options, 'speedThreshold',  1e-3);
options = setDefault(options, 'smoothingWindow', 1000);
options = setDefault(options, 'searchLimit',     288000);
options = setDefault(options, 'verbose',         true);

assert(isfile(nwbFile), 'Template:MissingFile', 'NWB file not found: %s', nwbFile);

SPEED = '/processing/behavior/behavioral_time_series/running_speed/data';
MOVIE = '/acquisition/raw_voltage_movie/data';

speed = double(h5read(nwbFile, SPEED));
speed = speed(:)';
nSample = numel(speed);

% MATLAB sees the HDF5 dims reversed, so this is [column row time].
movieDims = h5info(nwbFile, MOVIE).Dataspace.Size;
nCol = movieDims(1); nRow = movieDims(2); nFrameTotal = movieDims(3);
assert(nFrameTotal == nSample, 'Template:LengthMismatch', ...
    'Movie has %d frames but running speed has %d samples.', nFrameTotal, nSample);

searchLimit = min(options.searchLimit, nSample);

% ------------------------------------------------- find the first still run
switch lower(options.stillRule)
    case 'zeromotion'
        moving = speed(speed ~= 0);
        assert(~isempty(moving), 'Template:NeverMoves', ...
            'Running speed is zero everywhere; cannot derive an encoder floor.');
        stillThreshold = min(abs(moving));
        isStill  = abs(speed) < stillThreshold;
        % Kept only so the diagnostic plot can show the published statistic
        % alongside the one actually applied.
        smoothed = movmean(speed, options.smoothingWindow);
    case 'smoothedthreshold'
        smoothed = movmean(speed, options.smoothingWindow);
        stillThreshold = options.speedThreshold;
        isStill  = smoothed < stillThreshold;      % signed, as published
    otherwise
        error('Template:BadStillRule', ...
            'stillRule must be ''zeroMotion'' or ''smoothedThreshold''; got ''%s''.', ...
            options.stillRule);
end
startSample = firstRun(isStill(1:searchLimit), options.templateFrames);

if isnan(startSample)
    % Report what the data actually offers rather than just failing, because the
    % usual cause is a threshold that does not suit this session.
    longest = longestRun(isStill(1:searchLimit));
    error('Template:NoStillPeriod', ...
        ['No run of %d consecutive still frames in the first %d samples under ' ...
         'the ''%s'' rule (threshold %g). The longest run is %d and %d of %d ' ...
         'samples qualify. Either raise searchLimit or relax the rule.'], ...
        options.templateFrames, searchLimit, options.stillRule, stillThreshold, ...
        longest, sum(isStill(1:searchLimit)), searchLimit);
end
stopSample = startSample + options.templateFrames - 1;

% The source's search can land inside a plateau of still frames rather than at a
% true run boundary; this implementation cannot, but assert it regardless. A
% template built over frames where the animal was moving would bake that motion
% into the reference, and every later stage would inherit it silently.
assert(all(isStill(startSample:stopSample)), 'Template:RunNotStill', ...
    'Selected template window contains frames above the still threshold.');
assert(stopSample - startSample + 1 == options.templateFrames, ...
    'Template:RunWrongLength', 'Template window is %d frames, expected %d.', ...
    stopSample - startSample + 1, options.templateFrames);

if options.verbose
    fprintf(['template: frames %d-%d (%.1f s in), rule ''%s'' (threshold %.3g), ' ...
             '%d of %d searched frames still (%.1f%%)\n'], ...
        startSample, stopSample, startSample*1.27e-3, options.stillRule, ...
        stillThreshold, sum(isStill(1:searchLimit)), searchLimit, ...
        100*sum(isStill(1:searchLimit))/searchLimit);
end

% ------------------------------------------------------- median of those frames
% Read as [column row time] and permute to [row column time], the frame layout
% the rest of the tutorial uses.
block = h5read(nwbFile, MOVIE, [1 1 startSample], [nCol nRow options.templateFrames]);
block = permute(block, [2 1 3]);
template = single(median(single(block), 3));

info.startSample   = startSample;
info.stopSample    = stopSample;
info.nFrames       = options.templateFrames;
info.meanSpeed     = mean(speed(startSample:stopSample));
info.maxSpeed      = max(speed(startSample:stopSample));
info.fractionStill = sum(isStill(1:searchLimit)) / searchLimit;
info.searchLimit    = searchLimit;
info.smoothedSpeed  = smoothed;
info.stillRule      = options.stillRule;
info.stillThreshold = stillThreshold;
info.isStill        = isStill;
info.fractionStillSearched = sum(isStill(1:searchLimit)) / searchLimit;
info.options        = options;
end

%% ========================================================= local functions

function idx = firstRun(mask, runLength)
% First index at which `mask` is true for runLength consecutive samples, or NaN.
% Done by cumulative sum rather than by a loop so it stays fast on 288,000
% samples: the count of true values in any window is a difference of cumsums,
% and a full window means that difference equals its width.
%
% The source computes this differently, with a while loop over cumsum values
% that would misbehave on plateaus; the intent there is unambiguous and this is
% that intent, implemented directly.
mask = logical(mask(:)');
if numel(mask) < runLength, idx = NaN; return; end
c = [0 cumsum(mask)];
windowCount = c(runLength+1:end) - c(1:end-runLength);
idx = find(windowCount == runLength, 1, 'first');
if isempty(idx), idx = NaN; end
end

function n = longestRun(mask)
% Longest run of true values, for the diagnostic in the error message.
mask = logical(mask(:)');
if ~any(mask), n = 0; return; end
edges = diff([false mask false]);
n = max(find(edges == -1) - find(edges == 1));
end

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = value; end
end
