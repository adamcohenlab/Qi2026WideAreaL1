function result = demixCellNMF(patch, footprints, options)
% DEMIXCELLNMF  Separate one cell's voltage trace from its overlapping neighbours.
%
% result = DEMIXCELLNMF(patch, footprints, options)
%
% Wide-field voltage imaging of a dense population puts many somata within a point
% spread function of each other, so the pixels around any one cell also carry its
% neighbours' signal. Averaging those pixels with a fixed spatial mask gives a trace
% contaminated with the neighbours' spikes. This separates them.
%
% The three ideas
% ---------------
% 1. NEGATE. Voltron is a negative-going indicator: depolarisation makes the cell
%    dimmer. After negating the movie, spikes are positive excursions, and both the
%    spatial footprint and the temporal activity are genuinely non-negative
%    quantities. That is what makes a non-negative factorisation the right model
%    rather than an arbitrary one.
%
% 2. SPLIT INTO BANDS. Spikes live in the fast band, subthreshold voltage in the
%    slow band. They need different treatment, because only one of them carries the
%    landmarks that make cells distinguishable.
%
% 3. LEARN THE FOOTPRINT ON THE FAST BAND, THEN FREEZE IT. A cell's own spikes are
%    what make its footprint identifiable against its neighbours; the slow band has
%    no such landmarks and cannot constrain a spatial filter on its own. So the
%    footprint is estimated from the fast band, then held fixed while the slow band
%    is extracted through it. Using one spatial filter for both bands is also what
%    keeps the subthreshold trace on the same scale as the spiking trace, which is
%    what later normalisation by spike height depends on.
%
% Inputs
% ------
% patch : [nSamples x nPix] single or double
%     The motion-corrected movie restricted to this cell's pixel disk. In the
%     deposited NWB asset this is the mc_pixel_subset columns belonging to the cell.
%     NOTE this is the CORRECTED movie, not the raw one: demixing consumes the
%     output of the motion-correction stage, never the raw frames.
%
% footprints : [nPix x nComponent] double
%     Initial spatial footprints restricted to the same pixels, target cell FIRST,
%     overlapping neighbours after. From c_ftprnt, Gaussian-smoothed.
%
% options : struct
%   Required:
%     .pixelIndices    linear indices of the patch pixels in the full frame, which
%                      ALS_nonneg needs to rebuild images for its spatial smoothing
%     .frameSize       [nRows nColumns] of the full frame
%     .motionFreeFrames  frames free of large motion
%     .badFrames       frames to exclude, typically near chunk boundaries
%   Optional, defaulted to the published values:
%     .baselineWindow  F0 running-median window, samples. Default 24000 (30.5 s).
%     .bandWindow      slow/fast split running-median window. Default 60 (76 ms).
%     .pixelStdLimit   reject pixels whose moving standard deviation exceeds this,
%                      in raw camera counts. Default 5000. NOT PORTABLE: see note.
%     .chunkLength     window for the pixel-rejection moving std. Default 288000.
%     .initSamples     samples used for footprint initialisation. Default 2e5.
%     .maxIterOverlap  ALS iterations when neighbours are present. Default 5.
%     .maxIterAlone    ALS iterations with no neighbours. Default 3.
%     .fitFraction     fraction of usable frames the footprint is fitted on.
%                      Default 1/3.
%
% Note on pixelStdLimit
% ---------------------
% This rejects pixels at the edge of the field of view whose cells translate out of
% frame during large motion, producing enormous apparent variance. It is NOT a
% defective-pixel threshold. The value is in raw camera counts and was hand-picked
% for one recording, so it does not transfer to other data. See
% PIPELINE_BADPIXEL_PROPOSAL.md for the reasoning and a dimensionless alternative.
%
% Outputs
% -------
% result : struct
%     .trace             the demixed trace, fast plus slow band
%     .traceFast, .traceSlow   the two bands separately
%     .footprint         the learned spatial footprint, over patch pixels
%     .background        background spatial components
%     .backgroundTrace   background temporal components
%     .naiveTrace        mask-weighted average of the same pixels, for comparison.
%                        This is what you get WITHOUT demixing, and the point of
%                        the tutorial is that it contains the neighbours' spikes.
%     .rejectedPixels    indices of pixels zeroed by pixelStdLimit
%     .options           the resolved option set
%
% Requires: ALS_nonneg, ALS_nonneg_BackgroundUpdateOnly, cellSpkFtprntInit2, toimg,
% tovec. ALS_nonneg_BackgroundUpdateOnly calls CaImAn's HALS_temporal, so
% CaImAn-MATLAB must be on the path. See PIPELINE_SCRIPTS_AUDIT.md section H.
% MATLAB R2019b compatible.

%% ---------------------------------------------------------------- options
if nargin < 3, options = struct(); end
options = setDefault(options, 'baselineWindow', 8000*3);
options = setDefault(options, 'bandWindow', 60);
options = setDefault(options, 'pixelStdLimit', 5000);
options = setDefault(options, 'chunkLength', 288000);
options = setDefault(options, 'initSamples', 2e5);
options = setDefault(options, 'maxIterOverlap', 5);
options = setDefault(options, 'maxIterAlone', 3);
options = setDefault(options, 'fitFraction', 1/3);
for f = {'pixelIndices','frameSize','motionFreeFrames','badFrames'}
    assert(isfield(options, f{1}), 'Demix:MissingOption', 'options.%s is required.', f{1});
end

patch = double(patch);
[nSample, nPix] = size(patch);
nComponent = size(footprints, 2);
assert(size(footprints, 1) == nPix, 'Demix:FootprintSize', ...
    'footprints has %d rows but the patch has %d pixels.', size(footprints,1), nPix);
frameSize = options.frameSize;
pixelIdx  = options.pixelIndices(:);

isBad = false(1, nSample);
isBad(options.badFrames(options.badFrames >= 1 & options.badFrames <= nSample)) = true;

%% -------------------------------------------------- baseline and bad pixels
% A running median removes slow drift and photobleaching together, leaving the
% fluctuations around the local baseline.
patch = patch - movmedian(patch, options.baselineWindow, 1);

pixelStd = max(movstd(patch(~isBad, :), [options.chunkLength, 1]), [], 1);
rejected = find(pixelStd > options.pixelStdLimit);
patch(:, rejected) = 0;
result.rejectedPixels = rejected;

%% ---------------------------------------------------------- negate and split
% The sign flip. Everything downstream depends on it, including the negative-peak
% threshold calibration in the spike-finding stage.
flipped = -patch';                       % [nPix x nSample]
flipped(~isfinite(flipped)) = 0;

slow = movmedian(flipped, options.bandWindow, 2);
fast = flipped - slow;

%% ------------------------------------------------ the naive comparison
% A fixed mask-weighted average of the same pixels: what you get WITHOUT demixing.
% It is taken after baseline removal and after the negation, so it differs from the
% demixed trace only in how the pixels are combined, not in sign or preprocessing.
% Comparing against the un-negated patch instead would make the two anti-correlated
% and the comparison meaningless.
weights = footprints(:,1) / sum(footprints(:,1));
result.naiveTrace = weights' * flipped;

%% ------------------------------------------------- footprint initialisation
% Each component is initialised by a spike-triggered footprint estimate on the fast
% band, seeded with its own smoothed footprint against the pooled footprints of the
% others as background, then re-run seeded with its own output.
initFrames = intersect(1:min(options.initSamples, nSample), options.motionFreeFrames);
initFrames = initFrames(~isBad(initFrames));
assert(~isempty(initFrames), 'Demix:NoInitFrames', ...
    'No motion-free, non-boundary frames available for initialisation.');

A0 = zeros(nPix, nComponent);
for j = 1:nComponent
    target = footprints(:, j);
    target = target / sum(target);
    others = footprints(:, setdiff(1:nComponent, j));
    if isempty(others)
        others = ones(nPix, 1);
    else
        others = others / sum(others(:));
    end
    A0(:, j) = cellSpkFtprntInit2(fast(:, initFrames), [target others]);
    A0(:, j) = cellSpkFtprntInit2(fast(:, initFrames), [A0(:, j) others]);
end
A0 = max(A0, 0);

% The background basis is the neighbours plus a constant term, which absorbs
% whatever the named components do not explain.
if nComponent > 1
    b0 = [A0(:, 2:end), ones(nPix, 1)];
else
    b0 = ones(nPix, 1);
end

%% ------------------------------- factorise the fast band, then freeze A
% The footprint is fitted on motion-free, non-boundary frames, and only the leading
% fraction of them: enough spikes to pin the spatial filter, without paying for the
% whole recording.
fitFrames = intersect(options.motionFreeFrames, find(~isBad));
fitFrames = fitFrames(1:max(1, floor(numel(fitFrames) * options.fitFraction)));
if nComponent == 1
    maxIter = options.maxIterAlone;
else
    maxIter = options.maxIterOverlap;
end

[A, ~, bFast, ~] = ALS_nonneg(fast(:, fitFrames), single(A0(:,1)), single(b0), ...
    numel(fitFrames), pixelIdx, frameSize, maxIter);

% A is now frozen. Both bands are extracted through it, updating only the
% background, which is what keeps them on a common scale.
[~, traceFast, ~, bgFast] = ALS_nonneg_BackgroundUpdateOnly(fast, ...
    single(A(:,1)), single(bFast), 1, pixelIdx, frameSize);

[~, traceSlow, bSlow, bgSlow] = ALS_nonneg_BackgroundUpdateOnly(slow, ...
    single(A(:,1)), single([bFast(:, 1:end-1), ones(size(bFast,1), 1)]), 3, ...
    pixelIdx, frameSize);

%% ------------------------------------------------------------------ pack
result.traceFast       = traceFast;
result.traceSlow       = traceSlow;
result.trace           = traceFast + traceSlow;
result.footprint       = A(:,1);
result.background      = [bFast, bSlow(:, end)];
result.backgroundTrace = bgFast + bgSlow;
result.initialFootprints = A0;
result.options         = options;
end

%% ========================================================= local functions

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = value; end
end
