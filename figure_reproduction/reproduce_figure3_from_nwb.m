function results = reproduce_figure3_from_nwb(nwbFile27, nwbFile29, outputDirectory, panels)
% REPRODUCE_FIGURE3_FROM_NWB Recreate the data panels of manuscript Figure 3.
%
%   RESULTS = REPRODUCE_FIGURE3_FROM_NWB(NWBFILE27, NWBFILE29, OUTPUTDIRECTORY, PANELS)
%   reads the two documented sparsePulseHad NWB files, writes raw Figure 3
%   panels and a compact reference MAT file to OUTPUTDIRECTORY, and returns
%   the numerical results. It does not edit or execute the legacy scripts.
%
%   Session assignments are preserved from the manuscript workflow:
%       B, D, F      M-YQ0201-27
%       H waveform   M-YQ0201-27; H speed fit pools M-YQ0201-27 and -29
%       E, G          combined M-YQ0201-27 and M-YQ0201-29
%       J, K          M-YQ0201-29
%       M, N          M-YQ0201-27
%   The J/K versus M/N split follows line 1 of fig3_ipsp_amp_timeConst_pre_post.m
%   (2026-09-24): "use M-YQ0201-29 for ipsp amplitude, M-YQ0201-27 for decay
%   time constant".
%   Panels I (hand-drawn cartoons) and L/O (external bootstrap analyses) are
%   intentionally outside this reproduction scope.
%
%   Raw panel arrangement and styling are diagnostic, not a replacement for
%   the externally assembled manuscript figure.
%
%   PANELS is optional. Supply a panel letter, a comma-separated character
%   vector, or a cell array, for example 'F', 'E,G', or {'J','K'}. Omitting
%   PANELS (or supplying 'all') regenerates all implemented panels. Only the
%   NWB sessions required by the requested panel(s) are loaded.

% Placeholders: pass the paths to your local copies of these files instead.
if nargin < 1 || isempty(nwbFile27)
    nwbFile27 = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';
end
if nargin < 2 || isempty(nwbFile29)
    nwbFile29 = 'path/to/M_YQ0201_29_sparsePulseHad.nwb';
end
if nargin < 3 || isempty(outputDirectory)
    outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'figure3_nwb_reproduction');
end
if nargin < 4
    panels = [];
end
selectedPanels = normalizePanelSelection(panels);
if exist(outputDirectory, 'dir') ~= 7
    mkdir(outputDirectory);
end
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

need27 = any(ismember(selectedPanels, {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'M', 'N'}));
need29 = any(ismember(selectedPanels, {'E', 'G', 'H', 'J', 'K'}));
results = struct();
results.nwbFile27 = nwbFile27;
results.nwbFile29 = nwbFile29;
results.outputDirectory = outputDirectory;
results.selectedPanels = selectedPanels;
if need27
    session27 = loadFigure3SessionFromNWB(nwbFile27);
    assert(session27.numberOfCells == 320, 'Figure3NWB:CellCount', ...
        'Figure 3 expects 320 source cells per session.');
    data27 = selectFigureCells(session27);
    results.figureCellId27 = session27.cellId(data27.sourceCellIndex);
end
if need29
    session29 = loadFigure3SessionFromNWB(nwbFile29);
    assert(session29.numberOfCells == 320, 'Figure3NWB:CellCount', ...
        'Figure 3 expects 320 source cells per session.');
    data29 = selectFigureCells(session29);
    results.figureCellId29 = session29.cellId(data29.sourceCellIndex);
end

if any(strcmp(selectedPanels, 'A')), results.exampleA = plotSynapticExample(session27, outputDirectory); end
if any(strcmp(selectedPanels, 'C')), results.exampleC = plotGapJunctionExample(session27, outputDirectory); end
if any(strcmp(selectedPanels, 'B'))
    plotConnectivityMatrix(data27, 'synaptic', outputDirectory, 'figure3B_synaptic_matrix_raw.png');
    results.synapticMatrix27 = data27.metrics.synapticConnection;
end
if any(strcmp(selectedPanels, 'D'))
    plotConnectivityMatrix(data27, 'gap', outputDirectory, 'figure3D_gap_junction_matrix_raw.png');
    results.gapJunctionMatrix27 = data27.metrics.gapJunctionConnection;
end
if any(strcmp(selectedPanels, 'E'))
    results.ipspDistance = distanceSummary(data27, data29, 'ipsp');
    plotDistanceBoxes(results.ipspDistance, 'IPSP amplitude (norm)', outputDirectory, ...
        'figure3E_ipsp_distance_raw.png', []);
end
if any(strcmp(selectedPanels, 'G'))
    results.spikeletDistance = distanceSummary(data27, data29, 'spikelet');
    % fig3_gjConn.m line 464 fixes the y range for the spikelet panel.
    plotDistanceBoxes(results.spikeletDistance, 'Spikelet amplitude (norm)', outputDirectory, ...
        'figure3G_spikelet_distance_raw.png', [-0.05 0.15]);
end
if any(strcmp(selectedPanels, 'F'))
    [results.ipspTimeSeconds, results.ipspDistanceWaveforms] = distanceBinnedIPSP(session27, data27);
    plotDistanceWaveforms(results.ipspTimeSeconds, results.ipspDistanceWaveforms, ...
        'Peri-stimulus time (s)', 'Spike height', outputDirectory, 'figure3F_ipsp_waveforms_raw.png');
end
if any(strcmp(selectedPanels, 'H'))
    [results.spikeletTimeSeconds, results.spikeletDistanceWaveforms, speedFit27] = ...
        distanceBinnedSpikelets(session27, data27);
    [~, results.spikeletDistanceWaveforms29, speedFit29] = ...
        distanceBinnedSpikelets(session29, data29);
    results.speedFit = pooledSpeedFit(speedFit27, speedFit29);
    plotSpikeletDistanceAndSpeed(results.spikeletTimeSeconds, results.spikeletDistanceWaveforms, ...
        results.speedFit, outputDirectory, 'figure3H_spikelet_distance_raw.png');
end
if any(ismember(selectedPanels, {'J', 'K'}))
    results.prePostAmplitude = calculatePrePostComparisons(data29, 'amplitude');
    results.prePostAmplitude.statistics = plotAmplitudeComparisons(results.prePostAmplitude, selectedPanels, outputDirectory);
end
if any(ismember(selectedPanels, {'M', 'N'}))
    results.prePostTau = calculatePrePostComparisons(data27, 'tau');
    results.prePostTau.statistics = plotTauComparisons(results.prePostTau, selectedPanels, outputDirectory);
end

if numel(selectedPanels) == 12
    resultsFile = 'figure3_nwb_reference_results.mat';
else
    resultsFile = sprintf('figure3_nwb_reference_results_%s.mat', strjoin(selectedPanels, '_'));
end
save(fullfile(outputDirectory, resultsFile), 'results', '-v7.3');
fprintf('Figure 3 panel(s) %s written to %s\n', strjoin(selectedPanels, ','), outputDirectory);
end

function selectedPanels = normalizePanelSelection(panels)
implemented = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'M', 'N'};
if isempty(panels)
    selectedPanels = implemented;
    return
end
if ischar(panels)
    requested = regexp(upper(panels), '[A-Z]+', 'match');
elseif isstring(panels)
    requested = cellstr(upper(panels(:)));
elseif iscell(panels) && all(cellfun(@ischar, panels(:)))
    requested = upper(panels(:));
else
    error('Figure3NWB:PanelSelection', ...
        'PANELS must be a character vector, string array, or cell array of panel letters.');
end
requested = cellfun(@strtrim, requested, 'UniformOutput', false);
if any(strcmp(requested, 'ALL'))
    selectedPanels = implemented;
    return
end
if any(~ismember(requested, implemented))
    error('Figure3NWB:PanelSelection', ...
        'Supported panels are: %s.', strjoin(implemented, ', '));
end
selectedPanels = implemented(ismember(implemented, requested));
end

function data = selectFigureCells(session)
idx = session.figureCellIndices;
data = struct();
data.sourceCellIndex = idx(:);
data.coordinates = session.cellCoordinatePixels(idx, :);
data.pixelSizeMicrometers = session.pixelSizeMicrometers;
metricNames = fieldnames(session.metrics);
data.metrics = struct();
for iMetric = 1:numel(metricNames)
    name = metricNames{iMetric};
    data.metrics.(name) = session.metrics.(name)(idx, idx);
end
end

function example = plotSynapticExample(session, outputDirectory)
sorted = session.sortByX(1:session.numberOfCells);
pre = sorted(32);
post = sorted(24);
blue = session.readCellTimeSeries('blue', pre, []);
% fig3_synapConn.m lines 42-46: tBlueOn = find(diff(blue_all(ii,:) == blueVs(end)) > 0),
% i.e. the last frame BEFORE the pulse reaches its top level, with no +1.
% blueVs(end) is cell 1's top level; every cell shares it, so the cell's own
% maximum is used here to avoid reading a second full-length trace.
allOnsets = find(diff(blue == max(blue)) > 0);
numberOfTrials = numel(allOnsets);   % source line 193 prints this; caption n = 283
fprintf('  panel A: %d blue-pulse trials for the presynaptic cell\n', numberOfTrials);
nBack = 140; nFront = 140;
onsets = allOnsets(allOnsets > nBack & allOnsets + nFront <= numel(blue));
assert(numel(onsets) >= 7, 'Figure3NWB:ExampleAStimuli', ...
    'The panel-A candidate cell has fewer than seven blue-pulse onsets.');
center = onsets(7);
window = center - nBack:center + nFront;
% Source line 204 plots the trial-averaged blue waveform (get_sta_mat_2).
blueAverage = mean(blue(onsets(:) + (-nBack:nFront)), 1);
% Source line 218 plots trial 7 of the presynaptic normalized voltage.
preVoltage = session.readCellTimeSeries('voltage', pre, [window(1) window(end)]);
postIPSP = session.readPairwiseSeries('bluePulseIPSP', pre, post, []);
timeSeconds = session.pairwiseTimeSeconds;

h = figure('Color', 'w', 'Position', [100 100 680 650]);
ax(1) = subplot(3, 1, 1); plot(timeSeconds, blueAverage, 'b', 'LineWidth', 1.2); ylabel('Blue (a.u.)');
ax(2) = subplot(3, 1, 2); plot(timeSeconds, preVoltage, 'b'); ylabel('Cell i voltage (norm)');
ax(3) = subplot(3, 1, 3); plot(timeSeconds, squeeze(postIPSP), 'k', 'LineWidth', 1.2); ylabel('Cell j voltage (norm)'); xlabel('Peri-stimulus time (s)');
title(ax(3), sprintf('Trial average, n = %d trials', numberOfTrials));
% Source line 224.
linkaxes(ax, 'x'); xlim([-0.05 0.14]);
saveFigure(h, outputDirectory, 'figure3A_synaptic_example_raw.png');

example = struct('sourcePreCellIndex', pre, 'sourcePostCellIndex', post, ...
    'preCellId', session.cellId(pre), 'postCellId', session.cellId(post), ...
    'numberOfTrials', numberOfTrials, 'selectedPulseIndex', 7, 'pulseFrame', center);
end

function example = plotGapJunctionExample(session, outputDirectory)
sorted = session.sortByX(1:session.numberOfCells);
pre = sorted(3);
post = sorted(9);
selfSTA = session.readPairwiseSeries('hadamardCrossSpikeSTA', pre, pre, []);
postSTA = session.readPairwiseSeries('hadamardCrossSpikeSTA', pre, post, []);
timeSeconds = session.pairwiseTimeSeconds;
% Source line 189 prints length(spk_t_had_ortho{ci,cj}), the Hadamard-orthogonal
% spikes of cell i that enter this cross-STA (caption n = 7215). The ragged
% column follows the metric serialization: row (post-1)*nCells + pre.
spikeIndex = double(h5read(session.nwbFile, ...
    '/analysis/pairwise_connectivity_metrics/hadamard_orthogonal_spike_times_index'));
spikeCounts = diff([0; spikeIndex(:)]);
numberOfSpikes = spikeCounts((post - 1) * session.numberOfCells + pre);
fprintf('  panel C: %d Hadamard-orthogonal spikes of cell i\n', numberOfSpikes);

h = figure('Color', 'w', 'Position', [100 100 680 500]);
ax(1) = subplot(2, 1, 1); plot(timeSeconds, squeeze(selfSTA), 'b', 'LineWidth', 1.2); ylabel('Cell i STA (norm)');
ax(2) = subplot(2, 1, 2); plot(timeSeconds, squeeze(postSTA), 'k', 'LineWidth', 1.2); ylabel('Cell j STA (norm)'); xlabel('Time from spike (s)');
title(ax(2), sprintf('Spike triggered average, n = %d spikes', numberOfSpikes));
linkaxes(ax, 'x'); xlim([-0.04 0.04]);
saveFigure(h, outputDirectory, 'figure3C_gap_junction_example_raw.png');

example = struct('sourcePreCellIndex', pre, 'sourcePostCellIndex', post, ...
    'preCellId', session.cellId(pre), 'postCellId', session.cellId(post), ...
    'numberOfSpikes', numberOfSpikes);
end

function plotConnectivityMatrix(data, kind, outputDirectory, fileName)
[~, order] = sort(data.coordinates(:, 1));
switch kind
    case 'synaptic'
        amplitude = -data.metrics.ipspAmplitude(order, order);
        connection = data.metrics.synapticConnection(order, order);
        label = 'IPSP amplitude (norm)';
        cmap = winter(256);
    case 'gap'
        amplitude = data.metrics.spikeletAmplitude(order, order);
        connection = data.metrics.gapJunctionConnection(order, order);
        label = 'Spikelet amplitude (norm)';
        cmap = autumn(256);
    otherwise
        error('Figure3NWB:MatrixKind', 'Unknown connectivity-matrix kind.');
end

amplitude(1:size(amplitude, 1)+1:end) = nan;
h = figure('Color', 'w');
imageHandle = imagesc(amplitude);
axis image; xlabel('Presynaptic cell'); ylabel('Postsynaptic cell');
title(label);
colormap(cmap); colorbar;
if strcmp(kind, 'gap')
    % fig3_gjConn.m lines 257-259 map spikelet amplitude onto a fixed 0-0.05
    % colour range ('nearest','extrap'); the synaptic matrix uses its full range.
    caxis([0 0.05]);
end
% The NWB carries the binary calls (pMat/pMat_spklet), but not the original
% per-pair q/p values used for graded manuscript contrast. Binary alpha retains
% the scientific connectivity call without inventing unavailable p-values.
set(imageHandle, 'AlphaData', 0.12 + 0.88 * double(connection));
set(gca, 'Color', 'w', 'TickDir', 'in');
saveFigure(h, outputDirectory, fileName);
end

function summary = distanceSummary(data27, data29, kind)
% Since the [distfix 2026-09-24] edit both panels use the same six 56.7-um
% bins over 60-400 um:
%   E  fig3_synapConn.m line 624: distEdge = linspace(60,400,7)/dx (was 0-400)
%   G  fig3_gjConn.m    line 424: distEdge = linspace(60,400,7)/dx
% Do not replace this with 60:60:400: the colon operator stops at 360 because
% the next step would overshoot 400, silently dropping every 360-400 um pair.
switch kind
    case {'ipsp', 'spikelet'}
        edgesUm = linspace(60, 400, 7);
    otherwise
        error('Figure3NWB:DistanceKind', 'Unknown distance-summary kind.');
end
summary = struct('edgesUm', edgesUm, 'centersUm', mean([edgesUm(1:end-1); edgesUm(2:end)], 1), ...
    'mouse1', [], 'mouse2', [], 'scatter', [], 'kind', kind);
[summary.mouse1, summary.scatter27] = valuesByDistance(data27, kind, edgesUm);
[summary.mouse2, summary.scatter29] = valuesByDistance(data29, kind, edgesUm);
summary.scatter = struct('distanceUm', [summary.scatter27.distanceUm; summary.scatter29.distanceUm], ...
    'value', [summary.scatter27.value; summary.scatter29.value]);
end

function [byBin, raw] = valuesByDistance(data, kind, edgesUm)
distanceUm = data.metrics.distancePixels * data.pixelSizeMicrometers;
switch kind
    case 'ipsp'
        values = -data.metrics.ipspAmplitude;
        valid = data.metrics.ipspAmplitude > 1e-3;
    case 'spikelet'
        values = data.metrics.spikeletAmplitude;
        valid = true(size(values));
    otherwise
        error('Figure3NWB:DistanceKind', 'Unknown distance-summary kind.');
end
% fig3_synapConn.m line 595 and fig3_gjConn.m line 388 ([distfix 2026-09-24]):
% the scatter, like the boxes, now excludes pairs closer than 60 um (was > 0).
valid = valid & distanceUm > 60 & distanceUm < 400;
raw = struct('distanceUm', distanceUm(valid), 'value', values(valid));
byBin = cell(numel(edgesUm)-1, 1);
for iBin = 1:numel(byBin)
    inBin = valid & distanceUm >= edgesUm(iBin) & distanceUm < edgesUm(iBin+1);
    byBin{iBin} = values(inBin);
end
end

function plotDistanceBoxes(summary, yLabelText, outputDirectory, fileName, yLimits)
% Both panels keep the combined raw density scatter visible. The source's
% showScatter toggle (fig3_synapConn.m line 608, fig3_gjConn.m line 408) exists
% to export a box-only version for cosmetic figure assembly, not to change the
% analysis, so the reproduction always draws the underlying data.
h = figure('Color', 'w'); hold on;
centers = summary.centersUm;
boxWidth = 10;
% Both source scripts build the combined raw density scatter first, then place
% the two mouse box groups at distCtr*dx and distCtr*dx + boxW*2, i.e. the
% second group sits to the RIGHT of the bin centre rather than straddling it.
scatter_kde(summary.scatter.distanceUm, summary.scatter.value, 'filled', 'MarkerSize', 8);
% Source colours the density on a log scale with parula and adds a colorbar.
set(gca, 'ColorScale', 'log'); colormap(parula); colorbar;
plotBoxGroups(summary.mouse1, centers, boxWidth, [0.00 0.60 0.45]);
plotBoxGroups(summary.mouse2, centers + boxWidth * 2, boxWidth, [0.95 0.40 0.05]);
xlim([0 420]); xlabel('Inter-cell distance (um)'); ylabel(yLabelText);
if ~isempty(yLimits)
    ylim(yLimits);
end
set(gca, 'XScale', 'linear', 'YScale', 'linear', ...
    'XTick', 0:100:400, 'XTickLabel', {'0', '100', '200', '300', '400'});
% Source restyles the boxplot lines: uniform width, and dashed whisker and
% quartile lines redrawn solid.
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1);
set(findobj(gca, 'LineStyle', '--'), 'LineStyle', '-');
legend27 = plot(nan, nan, '-', 'Color', [0.00 0.60 0.45], 'LineWidth', 1.2);
legend29 = plot(nan, nan, '-', 'Color', [0.95 0.40 0.05], 'LineWidth', 1.2);
legend([legend27 legend29], {'M-YQ0201-27', 'M-YQ0201-29'}, 'Location', 'best');
saveFigure(h, outputDirectory, fileName);
end

function plotBoxGroups(valuesByBin, positions, boxWidth, color)
% The source pads each distance bin into one column of a NaN matrix and makes a
% single boxplot call (fig3_synapConn.m lines 628-656). Calling boxplot once per
% bin instead lets each call reset the axes it shares with the others.
nBin = numel(valuesByBin);
counts = cellfun(@numel, valuesByBin(:)');
if all(counts == 0)
    return
end
values = nan(max(counts), nBin);
for iBin = 1:nBin
    values(1:counts(iBin), iBin) = valuesByBin{iBin};
end
boxplot(values, 'positions', positions, 'widths', boxWidth, 'colors', color, ...
    'symbol', '', 'boxstyle', 'outline', 'plotstyle', 'traditional', ...
    'medianstyle', 'line');
end

function [timeSeconds, traces] = distanceBinnedIPSP(session, data)
n = session.numberOfCells;
wave = session.readPairwiseSeries('bluePulseIPSP', 1:n, 1:n, []);
crosstalk = session.readPairwiseSeries('blueLightCrosstalk', 1:n, 1:n, []);
idx = data.sourceCellIndex;
wave = wave(idx, idx, :);
crosstalk = crosstalk(idx, idx, :);
distanceUm = data.metrics.distancePixels * data.pixelSizeMicrometers;
edgesUm = linspace(60, 400, 11);
traces = nan(numel(session.pairwiseTimeSeconds), numel(edgesUm)-1);
for iBin = 1:numel(edgesUm)-1
    distanceMask = distanceUm >= edgesUm(iBin) & distanceUm < edgesUm(iBin+1);
    % This is the active F calculation in fig3_synapConn.m lines 779-795:
    % nanmean((crossBluePulseStaFNCorFull-fittedCrosstalkN).*distMask.*pMat).
    % Multiplication by pMat sets unconnected in-bin pairs to zero, which
    % remain in nanmean's denominator; do not average only connected pairs.
    numeratorMask = double(distanceMask) .* double(data.metrics.synapticConnection);
    traces(:, iBin) = meanOverPairMask(wave - crosstalk, numeratorMask, double(distanceMask));
end
timeSeconds = session.pairwiseTimeSeconds;
end

function [timeSeconds, traces, speedFit] = distanceBinnedSpikelets(session, data)
n = session.numberOfCells;
wave = session.readPairwiseSeries('hadamardCrossSpikeSTA', 1:n, 1:n, []);
idx = data.sourceCellIndex;
wave = wave(idx, idx, :);
distanceUm = data.metrics.distancePixels * data.pixelSizeMicrometers;
edgesUm = 60:30:420;
centersUm = mean([edgesUm(1:end-1); edgesUm(2:end)], 1);
traces = nan(numel(session.pairwiseTimeSeconds), numel(centersUm));
for iBin = 1:numel(centersUm)
    pairMask = double(distanceUm >= edgesUm(iBin) & distanceUm < edgesUm(iBin+1));
    traces(:, iBin) = meanOverPairMask(wave, pairMask);
end
traces = traces - mean(traces, 1);
% fig3_gjConn.m line 709 hard-codes dt = 1.27e-3 s; the NWB frame period must agree.
framePeriodSeconds = 1.27e-3;
assert(abs(session.pairwiseFramePeriodSeconds - framePeriodSeconds) < 1e-6, ...
    'Figure3NWB:FramePeriod', 'Pairwise frame period differs from the source 1.27 ms.');
speedFit = legacySpeedFit(traces, centersUm, framePeriodSeconds);
timeSeconds = session.pairwiseTimeSeconds;
end

function trace = meanOverPairMask(values, numeratorMask, denominatorMask)
if nargin < 3
    denominatorMask = numeratorMask;
end
nTime = size(values, 3);
trace = nan(nTime, 1);
for iTime = 1:nTime
    frame = values(:, :, iTime);
    % Match nanmean: nonfinite source values are omitted from both numerator
    % and denominator, while finite unconnected F pairs contribute zeros.
    valid = logical(denominatorMask) & isfinite(frame);
    denominator = sum(valid(:));
    if denominator == 0
        continue
    end
    weighted = frame .* numeratorMask;
    weighted(~valid) = 0;
    trace(iTime) = sum(weighted(:)) / denominator;
end
end

function speedFit = legacySpeedFit(traces, centersUm, framePeriodSeconds)
% Mirrors fig3_gjConn.m lines 735-780.  In particular, the cubic is fitted
% on tau(tPk + fitWin), not on fitWin itself. Since 2026-09-24 the source
% converts the peak offset to seconds (line 752, tPkAll = xTest(idxMax)*dt),
% fits distance in um against seconds and reports fitCoeffs(1)/1e6 m/s. The
% earlier frames-based /1000 overstated the speed by 1/1.27 (0.35 vs 0.28 m/s).
nTime = size(traces, 1);
centerFrame = ceil(nTime / 2);
coarseWindow = max(1, centerFrame-16):min(nTime, centerFrame+16);
localOffsets = -2:4;
tau = ((1:nTime) - centerFrame)';
peakOffsetFrames = nan(1, size(traces, 2));
for iBin = 1:size(traces, 2)
    local = zeros(nTime, 1);
    local(coarseWindow) = traces(coarseWindow, iBin);
    [~, peakFrame] = max(local);
    fitFrames = peakFrame + localOffsets;
    if any(fitFrames < 1 | fitFrames > nTime)
        continue
    end
    % Equivalent to: fit(tau(tPk + fitWin)', tr_distBin(tPk + fitWin,ii), 'poly3')
    coefficients = polyfit(tau(fitFrames)', traces(fitFrames, iBin)', 3);
    testOffsets = linspace(localOffsets(1), localOffsets(end), 10000);
    [~, maxIndex] = max(polyval(coefficients, testOffsets));
    peakOffsetFrames(iBin) = testOffsets(maxIndex);
end
peakDelaySeconds = peakOffsetFrames * framePeriodSeconds;
valid = isfinite(peakDelaySeconds);
linearCoefficients = polyfit(peakDelaySeconds(valid), centersUm(valid), 1);
speedFit = struct('distanceBinCentersUm', centersUm, ...
    'legacyPeakOffsetFrames', peakOffsetFrames, ...
    'peakDelaySeconds', peakDelaySeconds, ...
    'linearCoefficientsUmPerSecond', linearCoefficients, ...
    'speedMetersPerSecond', linearCoefficients(1) / 1e6, ...
    'framePeriodSeconds', framePeriodSeconds);
end

function speedFit = pooledSpeedFit(speedFit27, speedFit29)
% Mirrors fig3_gjConn.m lines 796-840: concatenate tPkAll/rBinsCtr from both
% sessions, then fit one line to the pooled points.
delay27 = speedFit27.peakDelaySeconds(:);
delay29 = speedFit29.peakDelaySeconds(:);
distance27 = speedFit27.distanceBinCentersUm(:);
distance29 = speedFit29.distanceBinCentersUm(:);
allDelays = [delay27; delay29];
allDistances = [distance27; distance29];
valid = isfinite(allDelays) & isfinite(allDistances);
linearCoefficients = polyfit(allDelays(valid), allDistances(valid), 1);
speedFit = struct('distanceBinCentersUm', allDistances, ...
    'peakDelaySeconds', allDelays, ...
    'linearCoefficientsUmPerSecond', linearCoefficients, ...
    'speedMetersPerSecond', linearCoefficients(1) / 1e6, ...
    'framePeriodSeconds', speedFit27.framePeriodSeconds, ...
    'mouse27PeakDelaysSeconds', delay27, 'mouse29PeakDelaysSeconds', delay29, ...
    'mouse27DistancesUm', distance27, 'mouse29DistancesUm', distance29, ...
    'mouse27SpeedMetersPerSecond', speedFit27.speedMetersPerSecond, ...
    'mouse29SpeedMetersPerSecond', speedFit29.speedMetersPerSecond);
end

function plotDistanceWaveforms(timeSeconds, traces, xLabelText, yLabelText, outputDirectory, fileName)
h = figure('Color', 'w');
plot(timeSeconds, traces, 'LineWidth', 1.1);
xlabel(xLabelText); ylabel(yLabelText); xlim([-0.05 0.2]);
saveFigure(h, outputDirectory, fileName);
end

function plotSpikeletDistanceAndSpeed(timeSeconds, traces, speedFit, outputDirectory, fileName)
h = figure('Color', 'w', 'Position', [100 100 900 370]);
subplot(1, 2, 1);
plot(timeSeconds, traces, 'LineWidth', 1.1); xlim([-0.02 0.02]);
xlabel('Peri-spike time (s)'); ylabel('Voltage (norm)');
subplot(1, 2, 2);
% Source lines 821-838: delays in ms, fit line drawn from 0 to the largest
% pooled delay, fixed axis limits.
scatter(speedFit.mouse27PeakDelaysSeconds * 1e3, speedFit.mouse27DistancesUm / 1000, 35, ...
    [0.00 0.60 0.45], 'filled'); hold on;
scatter(speedFit.mouse29PeakDelaysSeconds * 1e3, speedFit.mouse29DistancesUm / 1000, 35, ...
    [0.95 0.40 0.05], 'filled');
x = linspace(0, max(speedFit.peakDelaySeconds), 100);
plot(x * 1e3, polyval(speedFit.linearCoefficientsUmPerSecond, x) / 1000, 'r', 'LineWidth', 1.1);
xlim([0 1.5]); ylim([-0.1 0.5]);
xlabel('Peak delay (ms)'); ylabel('Distance (mm)');
title(sprintf('Speed = %.2f m/s', speedFit.speedMetersPerSecond));
legend({'M-YQ0201-27', 'M-YQ0201-29', 'Pooled fit'}, 'Location', 'best');
saveFigure(h, outputDirectory, fileName);
end

function prePost = calculatePrePostComparisons(data, kind)
% Mirrors fig3_ipsp_amp_timeConst_pre_post.m (2026-09-24). KIND is
% 'amplitude' (panels J/K, lines 135-153) or 'tau' (panels M/N, lines
% 260-284). Both blocks now select pairs with the SAME mask,
%   pMat & 60 < d < 400 um & tauSNR < 90th percentile & tau inside tauLim,
% and use it, NaN-coded as maskNan, as the leave-one-out weight, so a pair
% that fails any criterion is dropped from both the scatter and every
% average it would otherwise enter. The earlier source gated J/K without the
% tau-range terms and weighted the averages by pMat (J/K) or by the tau range
% alone (M/N); the n and R values changed accordingly.
pMat = data.metrics.synapticConnection;
tau = data.metrics.ipspDecay;
tauSE = data.metrics.ipspDecaySE;
distanceUm = data.metrics.distancePixels * data.pixelSizeMicrometers;
n = size(pMat, 1);
pMatNan = double(pMat); pMatNan(~pMat) = nan;
% seTau is NaN outside the significant pairs, so tauSNR is too, exactly as
% the source's nan-initialised tauSNR (lines 47 and 57).
tauSNR = tauSE ./ tau;
tauLimits = [min(tau(tau > 0)) max(tau(:))] + [1e-3 -1e-3];
mask = pMat & distanceUm < 400 & distanceUm > 60 & ...
    tauSNR < prctile(tauSNR(:), 90) & tau > tauLimits(1) & tau < tauLimits(2);
maskNan = nan(n); maskNan(mask) = 1;

switch kind
    case 'amplitude'
        values = data.metrics.ipspAmplitude;
        weight = maskNan;
    case 'tau'
        values = tau;
        weight = maskNan .* pMatNan;
    otherwise
        error('Figure3NWB:PrePostKind', 'Unknown pre/post comparison kind.');
end

% <x_ik>_{k~=j}: fixed presynaptic cell i, averaged over its other targets.
% <x_kj>_{k~=i}: fixed postsynaptic cell j, averaged over its other inputs.
presynapticAverage = nan(n); postsynapticAverage = nan(n);
for iCell = 1:n
    other = setdiff(1:n, iCell);
    presynapticAverage(:, iCell) = mean(weight(:, other) .* values(:, other), 2, 'omitnan');
    postsynapticAverage(iCell, :) = mean(weight(other, :) .* values(other, :), 1, 'omitnan');
end
fprintf('  panels %s: n = %d pairs\n', ternary(strcmp(kind, 'amplitude'), 'J/K', 'M/N'), sum(mask(:)));

prePost = struct();
prePost.kind = kind;
prePost.values = values;
prePost.presynapticAverage = presynapticAverage;
prePost.postsynapticAverage = postsynapticAverage;
prePost.mask = mask;
prePost.numberOfPairs = sum(mask(:));
prePost.tauSNR = tauSNR;
prePost.tauLimits = tauLimits;
end

function out = ternary(condition, a, b)
if condition
    out = a;
else
    out = b;
end
end

function statistics = plotAmplitudeComparisons(prePost, selectedPanels, outputDirectory)
showJ = any(strcmp(selectedPanels, 'J'));
showK = any(strcmp(selectedPanels, 'K'));
if showJ && showK
    h = figure('Color', 'w', 'Position', [100 100 900 400]);
    subplot(1, 2, 1); statistics.J = plotAmplitudePanel(prePost, true);
    subplot(1, 2, 2); statistics.K = plotAmplitudePanel(prePost, false);
    saveFigure(h, outputDirectory, 'figure3J_K_ipsp_amplitude_raw.png');
elseif showJ
    h = figure('Color', 'w'); statistics.J = plotAmplitudePanel(prePost, true);
    saveFigure(h, outputDirectory, 'figure3J_ipsp_amplitude_raw.png');
else
    h = figure('Color', 'w'); statistics.K = plotAmplitudePanel(prePost, false);
    saveFigure(h, outputDirectory, 'figure3K_ipsp_amplitude_raw.png');
end
end

function statistics = plotAmplitudePanel(prePost, isJ)
% Source lines 174-198 (J) and 206-224 (K): leave-one-out average on x, the
% pair's own amplitude on y.
if isJ
    average = prePost.presynapticAverage;
    xLabelText = '<w_{ik}>_{k\neq j}';
else
    average = prePost.postsynapticAverage;
    xLabelText = '<w_{kj}>_{k\neq i}';
end
statistics = plotLogComparison(average(prePost.mask), prePost.values(prePost.mask), ...
    prePost.numberOfPairs, xLabelText, 'w_{ij}', [7e-3 2e-1]);
end

function statistics = plotTauComparisons(prePost, selectedPanels, outputDirectory)
showM = any(strcmp(selectedPanels, 'M'));
showN = any(strcmp(selectedPanels, 'N'));
if showM && showN
    h = figure('Color', 'w', 'Position', [100 100 900 400]);
    subplot(1, 2, 1); statistics.M = plotTauPanel(prePost, true);
    subplot(1, 2, 2); statistics.N = plotTauPanel(prePost, false);
    saveFigure(h, outputDirectory, 'figure3M_N_ipsp_decay_raw.png');
elseif showM
    h = figure('Color', 'w'); statistics.M = plotTauPanel(prePost, true);
    saveFigure(h, outputDirectory, 'figure3M_ipsp_decay_raw.png');
else
    h = figure('Color', 'w'); statistics.N = plotTauPanel(prePost, false);
    saveFigure(h, outputDirectory, 'figure3N_ipsp_decay_raw.png');
end
end

function statistics = plotTauPanel(prePost, isM)
% Source lines 307-328 (M) and 338-358 (N): leave-one-out average on x.
if isM
    average = prePost.presynapticAverage;
    xLabelText = '<\tau_{ik}>_{k\neq j}';
else
    average = prePost.postsynapticAverage;
    xLabelText = '<\tau_{hj}>_{h\neq i}';
end
statistics = plotLogComparison(average(prePost.mask), prePost.values(prePost.mask), ...
    prePost.numberOfPairs, xLabelText, '\tau_{ij}', [5e-3 5e-1]);
end

function statistics = plotLogComparison(x, y, sourceCount, xLabelText, yLabelText, axisLimits)
% Mirrors the J/K/M/N source plotting blocks: density scatter, log-space
% fitlm of y on x, and the fitted line plus its default 95% confidence bounds.
% The source drops NaN pairs (a cell with no other qualifying pair has no
% leave-one-out average); requiring positive values as well only guards the
% logarithm and removes no pair in the current data.
valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;
x = x(valid); y = y(valid);
scatter_kde(x, y, 'filled', 'MarkerSize', 5); hold on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xticks([1e-2 1e-1]); yticks([1e-2 1e-1]);
[correlation, pValues] = corrcoef(log(x), log(y));
mdl = fitlm(log(x), log(y));
xFit = linspace(axisLimits(1), axisLimits(2), 200)';
[yFit, yCI] = predict(mdl, log(xFit));
plot(xFit, exp([yFit yCI]), 'k', 'LineWidth', 1.1);
xlabel(xLabelText); ylabel(yLabelText);
title(sprintf('R = %.2f, p = %.3e', correlation(1, 2), ...
    pValues(1, 2)));
fprintf('    %s vs %s: R = %.4f, p = %.4g, mask n = %d, plotted n = %d\n', ...
    yLabelText, xLabelText, correlation(1, 2), pValues(1, 2), sourceCount, numel(x));
xlim(axisLimits); ylim(axisLimits); daspect([1 1 1]);
colorbar;
% Keep source-mask count available in the axis metadata without replacing the
% legacy title, which reports R and p rather than n. The caption's n is the
% mask count (source line 159 / 296 prints sum(mask(:))).
statistics = struct('sourceMaskCount', sourceCount, ...
    'plottedFinitePositiveCount', numel(x), 'R', correlation(1, 2), ...
    'p', pValues(1, 2), 'slope', mdl.Coefficients.Estimate(2));
set(gca, 'UserData', statistics);
end

function saveFigure(figureHandle, outputDirectory, fileName)
set(figureHandle, 'PaperPositionMode', 'auto');
print(figureHandle, fullfile(outputDirectory, fileName), '-dpng', '-r300');
set(figureHandle, 'Visible', 'on');
drawnow;
end
