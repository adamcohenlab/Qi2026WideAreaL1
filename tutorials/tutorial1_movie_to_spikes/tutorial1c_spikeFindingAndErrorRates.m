function results = tutorial1c_spikeFindingAndErrorRates(nwbFile, cellIndices, doPlot)
% TUTORIAL 1c  Spike finding, and modelling the false positives and false negatives.
%
% results = TUTORIAL1C_SPIKEFINDINGANDERRORRATES(nwbFile, cellIndices, doPlot)
%
% Detects spikes in optical voltage traces and then models the detection error
% rates, using one NWB file and nothing else. No raw movie is needed: this tutorial
% starts from demixed traces, so it is the cheapest entry point into the pipeline.
%
% What this tutorial is about
% ---------------------------
% Optical spike detection has a threshold problem. Set it low and noise is counted
% as spikes; set it high and real spikes are missed. Usually the resulting error
% rates are unknown, and are quoted from simulations or from a handful of
% ground-truth recordings that may not resemble the data at hand.
%
% This pipeline gets both error rates from the data itself, using one observation:
%
%     Because the indicator is negative-going, the movie is negated before
%     demixing. After that flip, real spikes are strictly POSITIVE excursions.
%     The NEGATIVE peaks of the filtered trace therefore contain no spikes at
%     all -- they are a clean, in-band sample of this exact cell's noise.
%
% Two things follow.
%
%   1. Setting the threshold at the Nth largest negative peak means, by symmetry
%      of the noise, that N noise peaks are expected to exceed it on the positive
%      side. The false-positive count is CHOSEN, not estimated, and the threshold
%      adapts to each cell's own noise. (Stage 2 below.)
%
%   2. The retained negative peaks are a measured sample of the false-positive
%      height distribution. Holding that component fixed, the detected heights can
%      be modelled as a mixture of it and an unknown true-spike distribution. That
%      gives a per-spike false-positive probability, and -- as the mass of the
%      fitted true-spike component lying below threshold -- the false-negative
%      rate. (Stage 3 below.)
%
% Inputs
% ------
% nwbFile : char
%     Path to the raw-movie excerpt asset, which carries the demixed traces, the
%     per-cell motion traces and the published spike times in one file.
% cellIndices : numeric, optional
%     One-based cell indices to analyse. Default is the first five cells flagged
%     `tutorial_cell` in the file. Every stage is per-cell, so a handful is enough
%     to follow the method; the full analysis is the same loop over all cells.
% doPlot : logical, optional
%     Draw the diagnostic figures. Default true.
%
% Outputs
% -------
% results : struct array, one element per analysed cell, with fields
%     .cellIndex, .detection, .fit, and .comparison against the published spikes.
%
% Requires: adaptive_thresh, whitened_matched_filter, spike_trace (all local to
% this project), and Statistics and Machine Learning Toolbox.
% MATLAB R2019b compatible.

if nargin < 2, cellIndices = []; end
if nargin < 3 || isempty(doPlot), doPlot = true; end
assert(isfile(nwbFile), 'Tutorial1c:MissingFile', 'NWB file not found: %s', nwbFile);
addpath(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

%% ===================================================================== data
fprintf('Reading %s\n', nwbFile);
session = readTutorialSession(nwbFile);
fprintf('  %d cells, %d samples, %d chunks of %d, %.1f s at %.2f Hz\n', ...
    session.nCells, session.nSamples, session.nChunks, session.chunkLength, ...
    session.nSamples*session.framePeriod, 1/session.framePeriod);

if isempty(cellIndices)
    cellIndices = find(session.isTutorialCell);
    cellIndices = cellIndices(1:min(5, numel(cellIndices)))';
end
fprintf('  analysing cells %s\n', mat2str(cellIndices));

%% ============================================ stage 1: motion-free frames
% Every threshold and template below is estimated only from frames where the
% tissue was reasonably still, because residual motion produces trace excursions
% that a spike detector cannot distinguish from spikes. The multiplier is the
% spike-finding value; see selectMotionFreeFrames for why it differs by stage.
fprintf('\nStage 1: selecting motion-free frames\n');
motion = selectMotionFreeFrames(session.motionX, session.motionY, ...
    struct('multiplier', 2, 'reportMixtureCut', true));
fprintf('  kept %.1f%% of samples (cut = %.4g, published rule: 2 x mean)\n', ...
    100*motion.fractionKept, motion.threshold);
if ~isnan(motion.mixtureThreshold)
    fprintf(['  for comparison only, the mixture-model cut would be %.4g, ' ...
             'keeping %.1f%%\n'], motion.mixtureThreshold, 100*motion.mixtureFraction);
end

%% ================================ stages 2 and 3, per cell
detectOptions = struct( ...
    'framePeriod',      session.framePeriod, ...
    'chunkLength',      session.chunkLength, ...
    'motionFreeFrames', motion.frames, ...
    'falsePositiveRate', 500 / 3657.6);   % the published 500 events per session,
                                          % expressed as a rate so that shortening
                                          % the window does not change the science
fitOptions = struct('framePeriod', session.framePeriod);

results = struct('cellIndex', {}, 'detection', {}, 'fit', {}, 'comparison', {});
for k = 1:numel(cellIndices)
    ci = cellIndices(k);
    fprintf('\n--- cell %d (%d of %d) ---\n', ci, k, numel(cellIndices));
    trace = readTrace(nwbFile, ci, session.nSamples);

    % Stage 2. Whitened matched filter, then the negative-peak threshold.
    detection = detectSpikesMatchedFilter(trace, detectOptions);
    fprintf(['Stage 2: %d spikes; threshold %.3f set at negative peak #%d, ' ...
             'so ~%d false positives are expected over %.0f s (%.3f Hz)\n'], ...
        numel(detection.spikes), detection.threshold, ...
        detection.expectedFalsePositives, detection.expectedFalsePositives, ...
        session.nSamples*session.framePeriod, ...
        detection.expectedFalsePositives/(session.nSamples*session.framePeriod));

    % Stage 3. Mixture of truncated Gaussians, and both error rates.
    fit = fitSpikeHeightMixture(detection, fitOptions);
    fprintf(['Stage 3: mixing weight p = %.4f; true-spike height %.2f +/- %.2f; ' ...
             'mean per-spike P(false positive) = %.3f\n'], ...
        fit.mixingWeight, fit.trueSpikeParams(1), fit.trueSpikeParams(2), ...
        mean(fit.falsePositiveProbability));
    fprintf('         FALSE NEGATIVE RATE = %.1f%% of true spikes fall below threshold\n', ...
        100*fit.falseNegativeRate);
    fprintf('         %d of %d analysed spikes are confident (P(FP) < 0.05)\n', ...
        sum(fit.falsePositiveProbability < 0.05), numel(fit.falsePositiveProbability));

    % Regression check against the spike times stored in the file.
    comparison = compareToPublished(detection.spikes, session.publishedSpikes{ci}, ...
        session.timestamps);
    fprintf(['Check  : published %d spikes, detected %d; %.2f%% matched within ' ...
             '1 sample\n'], comparison.nPublished, comparison.nDetected, ...
        100*comparison.matchedFraction);

    results(k).cellIndex  = ci; %#ok<AGROW>
    results(k).detection  = detection;
    results(k).fit        = fit;
    results(k).comparison = comparison;
end

%% ================================================================== plots
if doPlot
    for k = 1:numel(results)
        plotSpikeFindingDiagnostics(results(k), session);
    end
end

falseNegative = arrayfun(@(r) r.fit.falseNegativeRate, results);
meanPFP       = arrayfun(@(r) mean(r.fit.falsePositiveProbability), results);
matchedPct    = arrayfun(@(r) r.comparison.matchedFraction, results);
fprintf('\nSummary across %d cells\n', numel(results));
fprintf('  false-negative rate  : %.1f%% to %.1f%% (median %.1f%%)\n', ...
    100*min(falseNegative), 100*max(falseNegative), 100*median(falseNegative));
fprintf('  mean P(false positive): %.3f to %.3f\n', min(meanPFP), max(meanPFP));
fprintf('  matched to published  : %.2f%% to %.2f%%\n', ...
    100*min(matchedPct), 100*max(matchedPct));
end

%% ========================================================= local functions

function session = readTutorialSession(nwbFile)
% Reads everything except the traces, which are read one cell at a time so the
% tutorial never holds the full array. Plain HDF5, no MatNWB required.
info = h5info(nwbFile, '/processing/ophys/Fluorescence/roi_response_series/data');
dims = info.Dataspace.Size;              % MATLAB order: [nCells nSamples]
session.nCells   = dims(1);
session.nSamples = dims(2);

t = h5read(nwbFile, '/processing/ophys/Fluorescence/roi_response_series/timestamps');
% The camera period is the modal sample interval. The mean would be wrong: the
% record carries real inter-chunk gaps, which must not be averaged into it.
session.framePeriod = median(diff(double(t(:))));

chunkIds = h5read(nwbFile, '/intervals/recording_chunks/id');
session.nChunks = numel(chunkIds);
assert(mod(session.nSamples, session.nChunks) == 0, 'Tutorial1c:RaggedChunks', ...
    'Sample count %d does not divide into %d chunks.', session.nSamples, session.nChunks);
session.chunkLength = session.nSamples / session.nChunks;

session.motionX = h5read(nwbFile, '/processing/ophys/motion_trace_x/data');
session.motionY = h5read(nwbFile, '/processing/ophys/motion_trace_y/data');

session.isTutorialCell = h5read(nwbFile, '/units/tutorial_cell') == 1;
spikeSeconds = h5read(nwbFile, '/units/spike_times');
spikeIndex   = h5read(nwbFile, '/units/spike_times_index');
session.publishedSpikes = splitRagged(double(spikeSeconds), double(spikeIndex));
session.timestamps = double(t(:))';
end

function trace = readTrace(nwbFile, cellIndex, nSamples)
trace = h5read(nwbFile, '/processing/ophys/Fluorescence/roi_response_series/data', ...
    [cellIndex 1], [1 nSamples]);
trace = double(reshape(trace, 1, []));
end

function cells = splitRagged(values, stopIndex)
% NWB ragged columns store one flat vector plus per-row end offsets.
cells = cell(numel(stopIndex), 1);
start = 1;
for ii = 1:numel(stopIndex)
    cells{ii} = values(start:stopIndex(ii));
    start = stopIndex(ii) + 1;
end
end

function comparison = compareToPublished(detectedSamples, publishedSeconds, timestamps)
% The file stores published spike times in seconds; detection returns sample
% indices. Comparison is at one-sample tolerance.
%
% Seconds are converted back to samples by looking each time up in the actual
% timestamp vector, NOT by dividing by the frame period. The recording clock is not
% uniform: the three chunk boundaries each advance by 20 microseconds rather than a
% full frame period, because the upstream concatenation discarded the physical pause
% between acquisition chunks. Dividing by the median period accumulates 3.75 ms of
% error, about three samples, which silently pushes later spikes outside a
% one-sample tolerance and reports roughly 50% agreement for a detection that is
% actually correct.
publishedSamples = interp1(timestamps, 1:numel(timestamps), ...
    publishedSeconds(:)', 'nearest', 'extrap');
comparison.nDetected  = numel(detectedSamples);
comparison.nPublished = numel(publishedSamples);
if isempty(publishedSamples) || isempty(detectedSamples)
    comparison.matchedFraction = NaN;
    comparison.matched = [];
    return
end
matched = false(1, numel(publishedSamples));
for ii = 1:numel(publishedSamples)
    matched(ii) = any(abs(detectedSamples - publishedSamples(ii)) <= 1);
end
comparison.matched = matched;
comparison.matchedFraction = mean(matched);
end
