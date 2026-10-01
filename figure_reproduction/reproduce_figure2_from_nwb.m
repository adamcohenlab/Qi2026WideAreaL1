function results = reproduce_figure2_from_nwb(nwbFile, outputDirectory)
% REPRODUCE_FIGURE2_FROM_NWB Recreate raw Figure 2 panels from an NWB session.
%
%   RESULTS = REPRODUCE_FIGURE2_FROM_NWB(NWBFILE, OUTPUTDIRECTORY) loads the
%   M-YQ0201-29 sparseOpto NWB file, recreates the analytical content of Figure
%   2, saves raw PNG panels and a MAT reference-results file in OUTPUTDIRECTORY,
%   and returns the numerical results. The original scripts in figure scripts/v25
%   are not read from or modified.
%
%   The default NWBFILE is the documented M-YQ0201-29 asset. OUTPUTDIRECTORY
%   defaults to figure2_nwb_reproduction in this directory.
%
%   The manuscript was assembled outside MATLAB. These raw outputs reproduce the
%   analysis and data relationships, not its pixel-level panel arrangement.

if nargin < 1 || isempty(nwbFile)
    % Placeholder: pass the path to your local copy of this file instead.
    nwbFile = 'path/to/M_YQ0201_29_sparseOpto.nwb';
end
if nargin < 2 || isempty(outputDirectory)
    outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'figure2_nwb_reproduction');
end
if exist(outputDirectory, 'dir') ~= 7
    mkdir(outputDirectory);
end

session = loadFigure2SessionFromNWB(nwbFile);
assert(session.numberOfCells == 320, 'Figure2NWB:UnexpectedCellCount', ...
    'Figure 2 expects 320 cells; this NWB file contains %d.', session.numberOfCells);

% This is the exact five-property order used in fig2_ephys.m.
properties = session.propertyMatrix;
propertyLabels = session.propertyLabels;

%% Figure 2C: Ward clustering and optional legacy UMAP
linkageTree = linkage(zscore(properties, [], 1), 'ward');
clusterId = cluster(linkageTree, 'maxclust', 3);
clusterCounts = arrayfun(@(iCluster) sum(clusterId == iCluster), 1:3);

umapCoordinates = [];
if exist('UMAP', 'class') || exist('UMAP', 'file')
    umap = UMAP('n_neighbors', 20, 'n_components', 2);
    umapCoordinates = umap.fit_transform(zscore(properties, [], 1));
else
    warning('Figure2NWB:MissingUMAP', [ ...
        'UMAP is not on the MATLAB path. The dendrogram and all non-UMAP panels ', ...
        'will be written, but Figure 2C will not include the UMAP scatter.']);
end
plotClustering(linkageTree, clusterId, umapCoordinates, outputDirectory);

%% Figure 2A: the legacy panel source starts at fig2_ephys.m line 157.
exampleCellIndex = 176;
% Read blue and voltage in complete protocol trials for all cells. This follows
% the NWB chunking (time x all cells) and avoids the prohibitively slow pattern
% of requesting one cell across the full recording repeatedly.
[stimulatedTrialMask, offStepSTA] = protocolStimulusSummary(session, 140, 390);
stimulatedTrials = find(stimulatedTrialMask(exampleCellIndex, :));
assert(~isempty(stimulatedTrials), 'Figure2NWB:MissingExampleStimulus', ...
    'The Figure 2A example cell has no blue-stimulated trials.');
exampleTrialIndex = stimulatedTrials(1);
[exampleVoltage, ~, exampleBlue, ~] = ...
    session.readTrial(exampleCellIndex, exampleTrialIndex);
[exampleRate, exampleStimulatedTrials] = trialAverageFiringRateFromMask( ...
    session, exampleCellIndex, stimulatedTrialMask(exampleCellIndex, :));
plotExamplePanel(session, exampleBlue, exampleVoltage, ...
    exampleRate, exampleCellIndex, exampleStimulatedTrials, outputDirectory);

%% Figure 2D: stored blue-evoked spike-triggered averages
sta = session.blueEvokedSpikeTriggeredAverage;
assert(size(sta, 2) == 281, 'Figure2NWB:UnexpectedSTAWindow', ...
    'Expected a 281-sample blue-evoked STA; found %d samples.', size(sta, 2));
staTimeSeconds = ((1:size(sta, 2)) - 141) * session.framePeriodSeconds;
clusterSpikeSTA = clusterMean(sta, clusterId);
plotClusterSpikeSTA(staTimeSeconds, clusterSpikeSTA, outputDirectory);

%% Figure 2E: blue-off-triggered averages from stored voltage and blue traces
offStepTimeSeconds = (-140:390) * session.framePeriodSeconds;
clusterOffStepSTA = clusterMean(offStepSTA, clusterId);
plotClusterOffStepSTA(offStepTimeSeconds, clusterOffStepSTA, outputDirectory);

%% Figure 2F: trial-average firing rates during blue-stimulated trials
clusterFiringRate = nan(3, session.framesPerTrial);
for iCluster = 1:3
    cellRates = nan(sum(clusterId == iCluster), session.framesPerTrial);
    memberIndices = find(clusterId == iCluster);
    for iMember = 1:numel(memberIndices)
        cellRates(iMember, :) = trialAverageFiringRateFromMask( ...
            session, memberIndices(iMember), stimulatedTrialMask(memberIndices(iMember), :))';
    end
    clusterFiringRate(iCluster, :) = mean(cellRates, 1);
end
trialTimeSeconds = (1:session.framesPerTrial) * session.framePeriodSeconds;
plotClusterFiringRate(trialTimeSeconds, clusterFiringRate, ...
    session.framePeriodSeconds, outputDirectory);

%% Figure 2G-K: property comparisons in manuscript display order
figureOrder = [4 2 3 5 1];
manuscriptLabels = {'ADPT', 'V_R_h (norm)', 'F_m_a_x (Hz)', 'tau_m (ms)', 'V_a_d_p (norm)'};
statistics = propertyStatistics(properties, clusterId);
plotPropertyComparisons(properties, clusterId, figureOrder, manuscriptLabels, statistics, outputDirectory);

results = struct();
results.nwbFile = nwbFile;
results.outputDirectory = outputDirectory;
results.cellId = session.cellId;
results.propertyMatrix = properties;
results.propertyLabels = propertyLabels;
results.clusterId = clusterId;
results.clusterCounts = clusterCounts;
results.umapCoordinates = umapCoordinates;
results.exampleCellIndex = exampleCellIndex;
results.exampleTrialIndex = exampleTrialIndex;
results.exampleStimulatedTrials = exampleStimulatedTrials;
results.clusterSpikeSTA = clusterSpikeSTA;
results.clusterOffStepSTA = clusterOffStepSTA;
results.clusterFiringRate = clusterFiringRate;
results.propertyStatistics = statistics;
save(fullfile(outputDirectory, 'figure2_nwb_reference_results.mat'), 'results', '-v7.3');
fprintf('Figure 2 raw reproduction written to %s\n', outputDirectory);
fprintf('Cluster counts: [%d %d %d]\n', clusterCounts);
end

function means = clusterMean(values, clusterId)
means = nan(3, size(values, 2));
for iCluster = 1:3
    means(iCluster, :) = mean(values(clusterId == iCluster, :), 1, 'omitnan');
end
end

function [stimulatedTrialMask, offStepSTA] = protocolStimulusSummary(session, numberBefore, numberAfter)
windowLength = numberBefore + numberAfter + 1;
stimulatedTrialMask = false(session.numberOfCells, session.numberOfTrials);
staSum = zeros(session.numberOfCells, windowLength);
staCount = zeros(session.numberOfCells, 1);
allCells = 1:session.numberOfCells;

for trialIndex = 1:session.numberOfTrials
    firstFrame = (trialIndex - 1) * session.framesPerTrial + 1;
    frameRange = [firstFrame, firstFrame + session.framesPerTrial - 1];
    blue = session.readCellTimeSeries('blue', allCells, frameRange);
    voltage = session.readCellTimeSeries('voltage', allCells, frameRange);
    stimulatedTrialMask(:, trialIndex) = any(blue > 0, 2);

    for cellIndex = 1:session.numberOfCells
        offFrames = find(diff(blue(cellIndex, :)) < 0) + 1;
        offFrames = offFrames(offFrames > numberBefore & ...
            offFrames <= session.framesPerTrial - numberAfter);
        for iEvent = 1:numel(offFrames)
            sampleIndices = offFrames(iEvent) + (-numberBefore:numberAfter);
            staSum(cellIndex, :) = staSum(cellIndex, :) + voltage(cellIndex, sampleIndices);
            staCount(cellIndex) = staCount(cellIndex) + 1;
        end
    end
    if mod(trialIndex, 25) == 0 || trialIndex == session.numberOfTrials
        fprintf('Processed blue protocol trial %d of %d\n', trialIndex, session.numberOfTrials);
    end
end

offStepSTA = staSum ./ staCount;
offStepSTA(staCount == 0, :) = nan;
end

function [firingRatePerFrame, trialIndices] = trialAverageFiringRateFromMask(session, cellIndex, trialMask)
trialIndices = find(trialMask);
assert(~isempty(trialIndices), 'Figure2NWB:NoStimulatedTrials', ...
    'Cell %d has no blue-stimulated trials.', cellIndex);
globalSpikeFrames = spikeTimesToFrames(session.spikeTimes{cellIndex}, session.voltageTimestamps);
spikeTrial = floor((globalSpikeFrames - 1) / session.framesPerTrial) + 1;
spikeFrameWithinTrial = mod(globalSpikeFrames - 1, session.framesPerTrial) + 1;
keep = ismember(spikeTrial, trialIndices);
uniqueTrialFrames = unique([spikeTrial(keep), spikeFrameWithinTrial(keep)], 'rows');
firingRatePerFrame = accumarray(uniqueTrialFrames(:, 2), 1, ...
    [session.framesPerTrial 1], @sum, 0) ./ numel(trialIndices);
end

function frameIndex = spikeTimesToFrames(spikeTimes, voltageTimestamps)
if isempty(spikeTimes)
    frameIndex = zeros(0, 1);
    return
end
frameIndex = interp1(voltageTimestamps, (1:numel(voltageTimestamps))', ...
    spikeTimes(:), 'nearest', 'extrap');
frameIndex = round(frameIndex);
frameIndex = frameIndex(frameIndex >= 1 & frameIndex <= numel(voltageTimestamps));
end

function statistics = propertyStatistics(properties, clusterId)
statistics = struct();
statistics.pairLabels = {'1 vs 2', '1 vs 3', '2 vs 3'};
statistics.pValues = nan(5, 3);
statistics.outlierMask = false(size(properties));
for iProperty = 1:5
    [~, outlierMask] = rmoutliers(properties(:, iProperty), 'mean', 'ThresholdFactor', 6);
    statistics.outlierMask(:, iProperty) = outlierMask;
    values = properties(~outlierMask, iProperty);
    groups = clusterId(~outlierMask);
    statistics.pValues(iProperty, :) = [ ...
        ranksum(values(groups == 1), values(groups == 2)), ...
        ranksum(values(groups == 1), values(groups == 3)), ...
        ranksum(values(groups == 2), values(groups == 3))];
end
end

function plotExamplePanel(session, blue, voltage, firingRate, cellIndex, stimulatedTrials, outputDirectory)
timeSeconds = (1:session.framesPerTrial) * session.framePeriodSeconds;
figureHandle = figure('Color', 'w', 'Position', [100 100 800 700]);
axesList(1) = subplot(4, 1, 1);
plot(timeSeconds, blue, 'b', 'LineWidth', 1.2);
ylabel('Blue (a.u.)');
axesList(2) = subplot(4, 1, 2);
plot(timeSeconds, voltage, 'k');
ylabel('Voltage (norm)');
axesList(3) = subplot(4, 1, 3);
numberRasterRows = max(10, numel(stimulatedTrials));
for iTrial = 1:numberRasterRows
    if iTrial <= numel(stimulatedTrials)
        [~, ~, ~, trialSpikes] = session.readTrial(cellIndex, stimulatedTrials(iTrial));
        plot(trialSpikes * session.framePeriodSeconds, iTrial + zeros(size(trialSpikes)), 'k.', 'MarkerSize', 6);
        hold on
    end
end
ylim([0 numberRasterRows + 1]);
ylabel('Trial');
axesList(4) = subplot(4, 1, 4);
plot(timeSeconds, movmean(firingRate / session.framePeriodSeconds, 20), 'k', 'LineWidth', 1.2);
ylabel('Rate (Hz)');
xlabel('Time (s)');
linkaxes(axesList, 'x');
xlim([0 timeSeconds(end)]);
savePanel(figureHandle, outputDirectory, 'figure2A_raw.png');
end

function plotClustering(linkageTree, clusterId, umapCoordinates, outputDirectory)
figureHandle = figure('Color', 'w', 'Position', [100 100 850 400]);
subplot(1, 2, 1);
dendrogram(linkageTree, 0, 'Orientation', 'left');
xlabel('Cluster distance');
set(gca, 'YTick', [], 'YTickLabel', []);
title('Ward hierarchical clustering');
subplot(1, 2, 2);
if isempty(umapCoordinates)
    axis off
    text(0.5, 0.5, 'UMAP implementation not on MATLAB path', ...
        'HorizontalAlignment', 'center');
else
    scatter(umapCoordinates(:, 1), umapCoordinates(:, 2), 24, clusterId, 'filled');
    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('Five intrinsic properties');
    colormap(lines(3));
end
savePanel(figureHandle, outputDirectory, 'figure2C_raw.png');
end

function plotClusterSpikeSTA(timeSeconds, clusterSTA, outputDirectory)
figureHandle = figure('Color', 'w');
baselineFrames = 126:138;
plot(timeSeconds, clusterSTA - mean(clusterSTA(:, baselineFrames), 2), 'LineWidth', 1.3);
xlabel('Time (s)'); ylabel('Voltage (norm)');
xlim([-0.02 0.02]);
legend({'Cluster 1', 'Cluster 2', 'Cluster 3'}, 'Location', 'northeast');
savePanel(figureHandle, outputDirectory, 'figure2D_raw.png');
end

function plotClusterOffStepSTA(timeSeconds, clusterSTA, outputDirectory)
figureHandle = figure('Color', 'w');
normalizer = max(clusterSTA, [], 2);
plot(timeSeconds, clusterSTA ./ normalizer, 'LineWidth', 1.3);
xlabel('Time from blue off step (s)'); ylabel('Voltage (norm.)');
xlim([-0.01 0.05]);
legend({'Cluster 1', 'Cluster 2', 'Cluster 3'}, 'Location', 'northeast');
savePanel(figureHandle, outputDirectory, 'figure2E_raw.png');
end

function plotClusterFiringRate(timeSeconds, clusterRate, framePeriodSeconds, outputDirectory)
figureHandle = figure('Color', 'w');
plot(timeSeconds, movmean(clusterRate / framePeriodSeconds, 25, 2)', 'LineWidth', 1.3);
xlabel('Time (s)'); ylabel('Firing rate (Hz)');
xlim([0 timeSeconds(end)]);
legend({'Cluster 1', 'Cluster 2', 'Cluster 3'}, 'Location', 'northwest');
savePanel(figureHandle, outputDirectory, 'figure2F_raw.png');
end

function plotPropertyComparisons(properties, clusterId, figureOrder, labels, statistics, outputDirectory)
figureHandle = figure('Color', 'w', 'Position', [100 100 1300 280]);
for iPanel = 1:numel(figureOrder)
    propertyIndex = figureOrder(iPanel);
    outlierMask = statistics.outlierMask(:, propertyIndex);
    values = properties(~outlierMask, propertyIndex);
    groups = clusterId(~outlierMask);
    subplot(1, 5, iPanel);
    boxplot(values, groups, 'Symbol', '');
    hold on
    for iCluster = 1:3
        groupValues = values(groups == iCluster);
        jitter = linspace(-0.15, 0.15, numel(groupValues));
        plot(iCluster + jitter, groupValues, '.', 'Color', [0.25 0.25 0.25], 'MarkerSize', 7);
    end
    ylabel(labels{iPanel});
    xlabel('Cluster');
    addSignificanceBars(statistics.pValues(propertyIndex, :), values);
end
savePanel(figureHandle, outputDirectory, 'figure2G_to_K_raw.png');
end

function addSignificanceBars(pValues, values)
yRange = range(values);
if yRange == 0
    yRange = 1;
end
yBase = max(values) + 0.08 * yRange;
pairs = [1 2; 1 3; 2 3];
for iPair = 1:3
    y = yBase + (iPair - 1) * 0.12 * yRange;
    plot(pairs(iPair, :), [y y], 'k-', 'LineWidth', 0.8);
    if pValues(iPair) <= 0.001
        label = '***';
    elseif pValues(iPair) <= 0.01
        label = '**';
    elseif pValues(iPair) <= 0.05
        label = '*';
    else
        label = 'n.s.';
    end
    text(mean(pairs(iPair, :)), y, label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end
ylim([min(values) - 0.05 * yRange, yBase + 0.45 * yRange]);
end

function savePanel(figureHandle, outputDirectory, fileName)
set(figureHandle, 'PaperPositionMode', 'auto');
print(figureHandle, fullfile(outputDirectory, fileName), '-dpng', '-r300');
% Keep the raw panel available for inspection after the reproduction finishes.
set(figureHandle, 'Visible', 'on');
drawnow;
end
