function result = detectSpikesMatchedFilter(trace, options)
% DETECTSPIKESMATCHEDFILTER  Optical spike detection with a calibrated false-positive count.
%
% result = DETECTSPIKESMATCHEDFILTER(trace, options)
%
% Detects spikes in one cell's demixed voltage trace using a whitened matched
% filter, and sets the detection threshold from the trace's own negative peaks so
% that the expected number of false positives is a chosen quantity rather than an
% estimated one.
%
% Why negative peaks work
% -----------------------
% The demixing negates the movie, because Voltron is a negative-going indicator.
% After that flip, true spikes are strictly positive excursions, so the *negative*
% peaks of the filtered trace contain no spikes at all: they are a direct, in-band,
% per-cell sample of the noise peak distribution, drawn from the very same trace.
%
% Sorting the negative peak heights and taking the Nth largest therefore gives a
% threshold at which, by symmetry of the noise, exactly N noise peaks are expected
% to exceed it on the positive side. The false-positive count is chosen in advance
% and the threshold adapts to whatever each cell's noise happens to be.
%
% Inputs
% ------
% trace : [1 x nSamples] double
%     One cell's demixed, motion-regressed voltage trace. Detection is invariant
%     to a constant per-cell scale factor, so a trace normalized by spike height
%     yields the same spike times as an unnormalized one.
%
% options : struct
%   Required:
%     .framePeriod        seconds per sample (1.27e-3 for this dataset)
%     .chunkLength        samples per acquisition chunk, for the boundary guards
%                         and the per-chunk noise equalization. Derive it from the
%                         recording rather than hard-coding it.
%     .motionFreeFrames   indices of frames judged free of large motion; the
%                         template and the adaptive thresholds are estimated only
%                         from these. See selectMotionFreeFrames.
%   Optional, defaulted to the published values:
%     .falsePositiveRate  expected false positives per second per cell.
%                         Default 0.136699, which is the published 500 events over
%                         the 3657.6 s ten-chunk session. THIS IS A RATE, NOT A
%                         COUNT: see the note below.
%     .templateHalfWindow template half-width in seconds. Default 10e-3.
%     .clipFirst          adaptive_thresh clip, first pass. Default 10.
%     .clipSecond         adaptive_thresh clip, second pass. Default 100.
%     .pnorm              adaptive_thresh pnorm. Default 0.25.
%     .minSpikes          adaptive_thresh minimum. Default 10.
%     .detectGuard        boundary guard for detection, samples. Default [-100 100].
%     .noiseGuard         boundary guard for the negative-peak sample, samples.
%                         Default [-5 10].
%
% Why falsePositiveRate is a rate
% -------------------------------
% The source pipeline hard-codes nFalsePos = 500, which is a count over a whole
% ~1 hour session, so it is a false-positive *rate* only once divided by duration.
% Running the same literal on a shorter excerpt silently inflates the rate: 500
% over a 24.4 minute excerpt is 2.5x the published rate. Expressing it as a rate
% and deriving the count keeps the science fixed when the window changes.
%
% Outputs
% -------
% result : struct
%     .spikes                 detected spike sample indices
%     .spikeValues            filtered amplitude at each detected spike
%     .threshold              detection threshold, in equalized filter units
%     .negativePeakValues     the expectedFalsePositives largest negative peaks,
%                             an empirical sample of the false-positive
%                             distribution, consumed by fitSpikeHeightMixture
%     .expectedFalsePositives count derived from falsePositiveRate and duration
%     .filtered               matched-filter output after per-chunk noise
%                             equalization; the trace the threshold applies to
%     .matchedFilterOutput    matched-filter output before equalization
%     .template               spike template from the second pass
%     .provisionalSpikes      first-pass spikes used to build the template
%     .options                the fully resolved option set, for the record
%
% Requires: adaptive_thresh, whitened_matched_filter, spike_trace.
% Statistics and Machine Learning Toolbox (prctile via adaptive_thresh).
% MATLAB R2019b compatible.

%% ---------------------------------------------------------------- options
if nargin < 2, options = struct(); end
options = setDefault(options, 'falsePositiveRate', 500 / 3657.6);
options = setDefault(options, 'templateHalfWindow', 10e-3);
options = setDefault(options, 'clipFirst', 10);
options = setDefault(options, 'clipSecond', 100);
options = setDefault(options, 'pnorm', 0.25);
options = setDefault(options, 'minSpikes', 10);
options = setDefault(options, 'detectGuard', [-100 100]);
options = setDefault(options, 'noiseGuard', [-5 10]);
for f = {'framePeriod','chunkLength','motionFreeFrames'}
    assert(isfield(options, f{1}), 'Spikes:MissingOption', ...
        'options.%s is required.', f{1});
end

dt          = options.framePeriod;
chunkLength = options.chunkLength;
trace       = reshape(double(trace), 1, []);
nTotal      = numel(trace);

% Three lengths matter and the source pipeline distinguishes them, so this does
% too. nTotal is the trace as supplied, which in the source is already padded to a
% whole number of chunks. nValid is how much of it is real signal rather than
% padding: the matched filter is built only from that part, because the pad is
% zeros and would corrupt the noise spectrum it whitens against. nPadded is the
% length the per-chunk noise equalization reshapes against.
options = setDefault(options, 'nValidSamples', nTotal);
nValid  = min(options.nValidSamples, nTotal);
nPadded = ceil(nTotal / chunkLength) * chunkLength;

halfWindow = round(options.templateHalfWindow / dt);
result.expectedFalsePositives = max(1, round(options.falsePositiveRate * nPadded * dt));

%% ------------------------------------------------- degenerate trace guard
if std(trace) == 0
    result.spikes = []; result.spikeValues = []; result.threshold = NaN;
    result.negativePeakValues = []; result.filtered = zeros(1, nPadded);
    result.matchedFilterOutput = zeros(1, nValid);
    result.template = zeros(1, 2*halfWindow+1);
    result.provisionalSpikes = []; result.options = options;
    return
end

%% ------------------------------------------------------------ pre-filter
% Chunk boundaries carry recording artifacts. They are interpolated across rather
% than deleted, so the sample grid stays uniform for the filtering that follows.
boundary = boundaryIndices(nPadded, chunkLength, options.detectGuard);
boundary = boundary(boundary <= nTotal);
keep     = setdiff(1:nTotal, boundary);
traceInterp = trace;
traceInterp(boundary) = interp1(keep, trace(keep), boundary, 'pchip', trace(end-5));

% Two successive short running-median high-passes, both over the whole supplied
% trace. The first defines the trace the matched filter runs on; the second is
% used only to find candidate peaks.
traceHP   = traceInterp - movmedian(traceInterp, 20);
traceFind = traceHP - movmedian(traceHP, 15);

%% ------------------------------------------- pass 1: provisional template
% Candidate peaks are taken only from motion-free frames. The published code adds
% the boundary indices back into the kept set here; that is inert, because those
% samples were replaced by interpolation above, and it is preserved rather than
% tidied so the numerics match the published run exactly.
evaluateOn = union(options.motionFreeFrames(:), boundary(:));
provisional = peaksAbove(zeroOutside(traceFind, evaluateOn, nTotal), ...
    options, options.clipFirst);
provisional = trimToWindow(provisional, halfWindow, nTotal);

matched = whitened_matched_filter(traceHP(1:nValid), provisional, -halfWindow:halfWindow);

%% ----------------------------------------------- pass 2: refined template
% The filter output covers only the valid samples, so it is zero-extended back to
% the full length before the frame mask is applied. The source relies on MATLAB
% growing the array implicitly at this point; doing it explicitly is the same
% arithmetic and does not hide where the zeros come from.
matchedFull = [matched, zeros(1, nTotal - numel(matched))];
refined = peaksAbove(zeroOutside(matchedFull, evaluateOn, nTotal), ...
    options, options.clipSecond);
refined = trimToWindow(refined, halfWindow, nTotal);

matched = whitened_matched_filter(traceHP(1:nValid), refined, -halfWindow:halfWindow);
[~, template] = spike_trace(traceHP, refined, halfWindow);

%% --------------------------------------- per-chunk noise equalization
% The indicator bleaches, so the noise amplitude drifts across the session. Each
% chunk is rescaled to the session-wide median absolute deviation of the filter
% output, which makes one threshold valid throughout.
padded = [matched, zeros(1, nPadded - numel(matched))];
byChunk = reshape(padded, chunkLength, []);
byChunk = byChunk ./ median(abs(byChunk), 1) * median(abs(matched));
equalized = reshape(byChunk, 1, []);

%% ------------------------------------- negative-peak threshold calibration
% A tighter guard than detection uses, matching the published code: boundary
% samples must not contaminate the noise sample the threshold is derived from.
noiseBoundary = boundaryIndices(nPadded, chunkLength, options.noiseGuard);
valid = true(1, nPadded);
valid(noiseBoundary) = false;

negativePeaks = findpeaks(-equalized .* valid, 'sortStr', 'descend');
nFalsePos = min(result.expectedFalsePositives, numel(negativePeaks));
assert(nFalsePos >= 1, 'Spikes:NoNegativePeaks', ...
    'The filtered trace has no negative peaks; the threshold cannot be calibrated.');
if nFalsePos < result.expectedFalsePositives
    warning('Spikes:FewNegativePeaks', ...
        ['Only %d negative peaks are available but %d false positives were ' ...
         'requested, so the threshold is the smallest negative peak and the ' ...
         'realised false-positive rate is below the requested one.'], ...
        nFalsePos, result.expectedFalsePositives);
end
threshold = negativePeaks(nFalsePos);

[positivePeaks, positiveLocs] = findpeaks(equalized .* valid);
spikes = positiveLocs(positivePeaks > threshold);

%% ------------------------------------------------------------------ pack
result.spikes              = spikes;
result.spikeValues         = equalized(spikes);
result.threshold           = threshold;
result.negativePeakValues  = negativePeaks(1:nFalsePos);
result.filtered            = equalized;
result.matchedFilterOutput = matched;
result.template            = template;
result.provisionalSpikes   = refined;
result.options             = options;
end

%% ========================================================= local functions

function idx = boundaryIndices(nTotal, chunkLength, guard)
% Samples within [guard(1), guard(2)] of every chunk boundary, including the
% start of the recording, clipped to the valid range.
offsets = guard(1):guard(2);
idx = (0:chunkLength:nTotal)' + offsets;
idx = unique(idx(:));
idx = idx(idx >= 1 & idx <= nTotal)';
end

function out = zeroOutside(signal, keepIdx, nSample)
% Suppress everything outside the frames the threshold may be estimated from, so
% findpeaks cannot pick candidates there.
out = signal;
out(setdiff(1:nSample, keepIdx)) = 0;
end

function locs = peaksAbove(signal, options, clip)
pks = findpeaks(signal);
thresh = adaptive_thresh(pks, clip, options.pnorm, options.minSpikes);
[~, locs] = findpeaks(signal, 'minpeakheight', thresh);
end

function locs = trimToWindow(locs, halfWindow, nSample)
% A spike too close to either end has no full template window around it.
locs = locs(locs - halfWindow > 0 & locs + halfWindow < nSample);
end

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = value; end
end
