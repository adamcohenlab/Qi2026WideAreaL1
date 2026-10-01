function results = reproduce_figure4_from_nwb(nwbFile27, nwbFile29, outputDirectory, panels)
% REPRODUCE_FIGURE4_FROM_NWB Recreate the data panels of manuscript Figure 4.
%
%   RESULTS = REPRODUCE_FIGURE4_FROM_NWB(NWBFILE27, NWBFILE29, OUTPUTDIRECTORY, PANELS)
%   reads the two documented behavior (spontaneous-activity) NWB files, writes
%   raw Figure 4 panels and a compact reference MAT file to OUTPUTDIRECTORY,
%   and returns the numerical results. It does not edit or execute the legacy
%   scripts.
%
%   Session assignments, confirmed by the author on 2026-09-18:
%       A, B, C, D, E, F, H, I    M-YQ0201-27
%       G                          combined M-YQ0201-27 and M-YQ0201-29
%       STATS                      combined, as the source statistics block does
%   Panel J is a hand-drawn schematic and is intentionally outside this
%   reproduction scope.
%
%   Raw panel arrangement and styling are diagnostic, not a replacement for
%   the externally assembled manuscript figure.
%
%   PANELS is optional. Supply a panel letter, a comma-separated character
%   vector, or a cell array, for example 'H', 'B,C,I', or {'E','F'}. The extra
%   selector 'STATS' reproduces the reported length constants, short-range
%   correlations, correlation widths, and cross-STA depolarization size from
%   both sessions. Omitting PANELS (or supplying 'all')
%   regenerates everything implemented. Only the NWB sessions required by the
%   requested panel(s) are loaded.
%
%   Where the analysis starts
%   -------------------------
%   fig4_xcorr.m has two stages. Stage 1 (its lines 144-248) computes the
%   pairwise voltage cross-correlations and spike-triggered averages and saves
%   the normalized arrays; that stage costs hours and its outputs are exactly
%   what the NWB files store. This script therefore starts at stage 2, the
%   distance binning, for panels B, C, D, H, and I. Panels A, E, F, and G are
%   recomputed from the stored subthreshold voltage, spike times, and
%   behavioral traces.
%
%   Two quiet-state definitions, now unified
%   ----------------------------------------
%   Historically fig4_xcorr.m used a two-component Gaussian-mixture boundary on
%   log(whisking motion) while fig4_Pspk_V_def_VShuffle_GLMShuffle.m used a
%   median threshold, so panels B/C/D/H/I and panels E/F/G rested on different
%   quiet masks. The author unified both scripts on the Gaussian-mixture
%   boundary on 2026-09-18 and reported no noticeable change in the figures.
%   This script uses the single stored quiet_mask throughout, which is that
%   definition. Note that the stored mask comes from an independent unseeded
%   fitgmdist draw in dataPrepsForMAT.m, so it is the same definition as, but
%   not bit-identical to, the mask behind the stored correlation arrays.
%
%   Preserved source behavior
%   -------------------------
%   The following are copied deliberately and must not be "cleaned up":
%     * 31 distance columns from rEdge = [0 1 linspace(60,1800,30)]/dx, with
%       column 2 (6.5-60 um) blanked as the caption's omitted 0-60 um bin and
%       column 1 (self pairs) blanked in the zero-lag summary vectors only.
%     * Panel B and I plot 7 distance traces (including column 1); panel C
%       plots 6.
%     * The chunk-boundary exclusion guard is +/-1000 frames for panels E and F
%       but +/-5000 frames for panel G. The widths differ in the source and are
%       kept per panel.
%     * Panel G excludes pairs closer than 60 um, which the caption does not
%       state.
%     * Panel F seeds rng(0,'twister') before randperm, so its seven
%       highlighted cells are reproducible.
%     * The STATS length constants fit log values with polyfit over the
%       non-blank bins. Only the synchrony fit is restricted to r < 800 um;
%       the 2.5 ms-lag spike fit uses every bin, as in the 2026-09-22 source.
%     * Panels C and D show the spike-rate cross-correlogram in Hz (divided by
%       dt), as in the 2026-09-22 source. The stored values stay per frame.
%
%   Line references to the legacy scripts are to their 2026-09-22 versions.
%
%   See also LOADFIGURE4SESSIONFROMNWB, REPRODUCE_FIGURE3_FROM_NWB.

% Placeholders: pass the paths to your local copies of these files instead.
if nargin < 1 || isempty(nwbFile27)
    nwbFile27 = 'path/to/M_YQ0201_27_behavior.nwb';
end
if nargin < 2 || isempty(nwbFile29)
    nwbFile29 = 'path/to/M_YQ0201_29_behavior.nwb';
end
if nargin < 3 || isempty(outputDirectory)
    outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'figure4_nwb_reproduction');
end
if nargin < 4
    panels = [];
end
selectedPanels = normalizePanelSelection(panels);
if exist(outputDirectory, 'dir') ~= 7
    mkdir(outputDirectory);
end
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'));

need27 = any(ismember(selectedPanels, {'A','B','C','D','E','F','G','H','I','STATS'}));
need29 = any(ismember(selectedPanels, {'G','STATS'}));
needStage2Session27 = any(ismember(selectedPanels, {'B','C','D','I','STATS'}));
needStage2Session29 = any(strcmp(selectedPanels, 'STATS'));
needVoltage27 = any(ismember(selectedPanels, {'E','F','G'}));
needVoltage29 = any(strcmp(selectedPanels, 'G'));

results = struct();
results.nwbFile27 = nwbFile27;
results.nwbFile29 = nwbFile29;
results.outputDirectory = outputDirectory;
results.selectedPanels = selectedPanels;
results.generatedOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');

if need27
    session27 = loadFigure4SessionFromNWB(nwbFile27);
    results.cellCount27 = session27.numberOfCells;
    results.correlationCellCount27 = session27.correlationCellCount;
    results.framePeriodSeconds = session27.framePeriodSeconds;
    results.pixelSizeMicrometers = session27.pixelSizeMicrometers;
end
if need29
    session29 = loadFigure4SessionFromNWB(nwbFile29);
    results.cellCount29 = session29.numberOfCells;
    results.correlationCellCount29 = session29.correlationCellCount;
end

% ---------------------------------------------------- stage 2 distance binning
if needStage2Session27
    stage2_27 = distanceBinnedCorrelations(session27);
    results.distanceBins = stage2_27.bins;
    results.distanceAveraged27 = stage2_27.averaged;
    results.zeroLagSummary27 = stage2_27.summary;
end
if needStage2Session29
    stage2_29 = distanceBinnedCorrelations(session29);
    results.distanceAveraged29 = stage2_29.averaged;
    results.zeroLagSummary29 = stage2_29.summary;
end

% ------------------------------------------------------------------- panel A
if any(strcmp(selectedPanels, 'A'))
    results.panelA = plotVoltageAndSpikes(session27, outputDirectory);
end

% ------------------------------------------------------------ panels B, C, I
if any(strcmp(selectedPanels, 'B'))
    plotDistanceStateMaps(stage2_27, 'xcorrV', [1 4:5:30], 1, ...
        'Lag (s)', 'Correlation', outputDirectory, ...
        'figure4B_voltage_xcorr_distance_raw.png');
end
if any(strcmp(selectedPanels, 'C'))
    % Per-frame spike probability shown in Hz (fig4_xcorr.m lines 398-415).
    plotDistanceStateMaps(stage2_27, 'staRate', 4:5:30, ...
        1 / stage2_27.framePeriodSeconds, 'Peri-spike time (s)', ...
        'P(post spk | pre spk) (Hz)', outputDirectory, ...
        'figure4C_spike_xcorr_distance_raw.png');
end
if any(strcmp(selectedPanels, 'I'))
    plotDistanceStateMaps(stage2_27, 'staV', [1 4:5:30], 1, ...
        'Peri-spike time (s)', 'Voltage (norm)', outputDirectory, ...
        'figure4I_sta_voltage_distance_raw.png');
end

% ------------------------------------------------------------------- panel D
if any(strcmp(selectedPanels, 'D'))
    results.panelD = plotDistanceDecayByLag(stage2_27, outputDirectory);
end

% --------------------------------------------------------- panels E, F and G
% One streaming pass over the subthreshold voltage serves all three panels.
if needVoltage27
    voltage27 = accumulateVoltageStatistics(session27, any(strcmp(selectedPanels, 'G')));
end
if any(ismember(selectedPanels, {'E','F'}))
    histograms27 = voltageHistograms(voltage27);
    results.panelE = histograms27;
end
if any(strcmp(selectedPanels, 'E'))
    plotVoltageHistograms(histograms27, outputDirectory);
end
if any(strcmp(selectedPanels, 'F'))
    results.panelF = plotExcitabilityCurves(histograms27, session27, outputDirectory);
end

if any(strcmp(selectedPanels, 'G'))
    joint27 = neighborVoltageJoint(voltage27);
    joint29 = [];
    if needVoltage29
        voltage29 = accumulateVoltageStatistics(session29, true);
        joint29 = neighborVoltageJoint(voltage29);
        clear voltage29
    end
    results.panelG = combineNeighborJoint(joint27, joint29);
    plotNeighborVoltageJoint(results.panelG, outputDirectory);
end
if needVoltage27
    clear voltage27
end

% ------------------------------------------------------------------- panel H
if any(strcmp(selectedPanels, 'H'))
    results.panelH = plotStaExamples(session27, outputDirectory);
end

% -------------------------------------------------------- reported statistics
if any(strcmp(selectedPanels, 'STATS'))
    results.statistics = reportedStatistics(stage2_27, stage2_29);
end

if numel(selectedPanels) == 10
    resultsFile = 'figure4_nwb_reference_results.mat';
else
    resultsFile = sprintf('figure4_nwb_reference_results_%s.mat', ...
        strjoin(selectedPanels, '_'));
end
resultsPath = fullfile(outputDirectory, resultsFile);
writeWithRetry(@() saveResults(resultsPath, results), resultsPath);
fprintf('Figure 4 panel(s) %s written to %s\n', strjoin(selectedPanels, ','), outputDirectory);
end

% =========================================================================
% Panel selection
% =========================================================================

function selectedPanels = normalizePanelSelection(panels)
implemented = {'A','B','C','D','E','F','G','H','I','STATS'};
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
    error('Figure4NWB:PanelSelection', ...
        'PANELS must be a character vector, string array, or cell array of panel letters.');
end
requested = cellfun(@strtrim, requested, 'UniformOutput', false);
if any(strcmp(requested, 'ALL'))
    selectedPanels = implemented;
    return
end
if any(strcmp(requested, 'J'))
    error('Figure4NWB:PanelSelection', ...
        'Panel J is a hand-drawn schematic and is outside this reproduction scope.');
end
if any(~ismember(requested, implemented))
    error('Figure4NWB:PanelSelection', ...
        'Supported selectors are: %s.', strjoin(implemented, ', '));
end
selectedPanels = implemented(ismember(implemented, requested));
end

% =========================================================================
% Time masks
% =========================================================================

function masks = buildTimeMasks(session)
% Reproduces the masking used by both Figure 4 source scripts.
%
% base      : the analyzed window, frames 5000 to round(2000/dt)
% tRun      : locomotion guard, |running speed| > 1e-3 dilated by +/-1000 frames
% guardN    : frames within +/-N of an acquisition-chunk boundary, to be excluded
% analysisEF: base, both guards at +/-1000, no running, quiet     (panels E, F)
% analysisG : the same with the +/-5000 chunk guard               (panel G)
%
% analysisG is a subset of analysisEF because the wider chunk guard removes
% strictly more frames, which lets panel G reuse the panel E/F voltage matrix.
dt = session.framePeriodSeconds;
nFramesTotal = session.numberOfVoltageFrames;

nWin = 1000;
kernel = zeros(1, nWin * 2 + 1);
kernel(nWin + 1 + (-1000:1000)) = 1;
tRun = conv(abs(session.runningSpeed) > 1e-3, kernel, 'same') > 0;

base = false(1, nFramesTotal);
base(5e3:round(2000 / dt)) = true;

masks = struct();
masks.base = base;
masks.tRun = tRun;
masks.quiet = session.quietMask(:)';
masks.whisking = session.whiskingMask(:)';
masks.guard1000 = chunkBoundaryGuard(session, 1000);
masks.guard5000 = chunkBoundaryGuard(session, 5000);
masks.analysisEF = base & ~masks.guard1000 & ~tRun & masks.quiet;
masks.analysisG = base & ~masks.guard5000 & ~tRun & masks.quiet;
assert(~any(masks.analysisG & ~masks.analysisEF), 'Figure4NWB:MaskNesting', ...
    'The panel G mask must be a subset of the panel E/F mask.');
end

function guard = chunkBoundaryGuard(session, halfWidth)
% Source form: t_noise = (1:nFrame:nFramesTotal+1) + (-halfWidth:halfWidth)',
% then clamped into range. Every acquisition-chunk boundary is excluded.
nFramesTotal = session.numberOfVoltageFrames;
boundaries = 1:session.framesPerChunk:(nFramesTotal + 1);
offsets = (-halfWidth:halfWidth)';
indices = boundaries + offsets;
indices = max(min(indices(:), nFramesTotal), 1);
guard = false(1, nFramesTotal);
guard(indices) = true;
end

% =========================================================================
% Stage 2: distance-binned population averages
% =========================================================================

function stage2 = distanceBinnedCorrelations(session)
% Bins every stored pairwise array by pair distance, following
% fig4_xcorr.m lines 250-287, then forms the zero-lag summary vectors used by
% panels D and by the reported statistics. The cross-STA depolarization size
% (lines 578-582) is also taken here, because it needs the full pairwise quiet
% STA voltage array and that array is already in memory.
dx = session.pixelSizeMicrometers;
distanceMatrix = session.correlationDistancePixels;
nLag = session.numberOfLagFrames;
nLagSamples = session.numberOfLags;

rEdge = [0 1 linspace(60, 1800, 30)] / dx;
rCtr = mean(rEdge([1:end-1; 2:end]), 1);
nDistBin = numel(rEdge) - 1;

seriesFields = {'xcorrV', 'staV', 'staRate'};
seriesNames = {'voltage_cross_correlation', 'spike_triggered_voltage', ...
    'spike_triggered_spike_rate'};

averaged = struct();
for iSeries = 1:numel(seriesFields)
    binned = zeros(nLagSamples, nDistBin, 2);
    for state = 1:2
        % One [nCells x nCells x nLagSamples] array at a time keeps the
        % resident footprint near 140 MB instead of loading all six.
        values = session.readCorrelationSeries(seriesNames{iSeries}, state);
        for iBin = 1:nDistBin
            distanceMask = distanceMatrix >= rEdge(iBin) & distanceMatrix < rEdge(iBin + 1);
            maskNan = nan(size(distanceMask));
            maskNan(distanceMask) = 1;
            binned(:, iBin, state) = squeeze(nanmean(values .* maskNan, [1 2]));
        end
        if strcmp(seriesFields{iSeries}, 'staV') && state == 1
            % Range of each pair's quiet STA over lags -100 to -2 frames, that
            % is the depolarization before the trigger spike, averaged over
            % pairs 60-400 um apart. The STA is already in units of the target
            % cell's spike height.
            pairMaskNan = nan(size(distanceMatrix));
            pairMaskNan(distanceMatrix > 60 / dx & distanceMatrix < 400 / dx) = 1;
            crossStaDepolarization = nanmean(range(values(:, :, nLag + 1 + (-100:-2)), 3) ...
                .* pairMaskNan, [1 2]);
        end
        clear values
    end
    averaged.(seriesFields{iSeries}) = binned;
end

% Source order matters: the zero-lag self-pair rate is zeroed first, then the
% 6.5-60 um column is blanked as the caption's omitted 0-60 um bin.
averaged.staRate(nLag + 1, 1, :) = 0;
averaged.xcorrV(:, 2, :) = nan;
averaged.staV(:, 2, :) = nan;
averaged.staRate(:, 2, :) = nan;

zeroLag = session.zeroLagIndex;
summary = struct();
summary.xcorrVT0Val = blankSelfColumn(squeeze(averaged.xcorrV(zeroLag, :, :))');
summary.staVT0Val = blankSelfColumn(squeeze(averaged.staV(zeroLag, :, :))');
summary.staFRT0Val = blankSelfColumn(squeeze(averaged.staRate(zeroLag, :, :))');
summary.staFRT2Val = blankSelfColumn(squeeze(averaged.staRate(zeroLag + 2, :, :))');

% Synchrony: the zero-lag peak over +/-1 frames minus the larger of the two
% flanking baselines at 3-5 frames on either side (fig4_xcorr.m lines 425-427).
peak = max(averaged.staRate(zeroLag + (-1:1), :, :), [], 1);
flankBefore = mean(averaged.staRate(zeroLag + (-5:-3), :, :), 1);
flankAfter = mean(averaged.staRate(zeroLag + (3:5), :, :), 1);
baseline = max(cat(1, flankBefore, flankAfter), [], 1);
summary.staSync = blankSelfColumn(squeeze(peak - baseline)');

stage2 = struct();
stage2.bins = struct('rEdgePixels', rEdge, 'rCenterPixels', rCtr, ...
    'rCenterMicrometers', rCtr * dx, 'numberOfBins', nDistBin, ...
    'pixelSizeMicrometers', dx);
stage2.averaged = averaged;
stage2.summary = summary;
stage2.crossStaDepolarizationQuiet = crossStaDepolarization;
stage2.lagSeconds = session.lagSeconds;
stage2.framePeriodSeconds = session.framePeriodSeconds;
stage2.zeroLagIndex = zeroLag;
stage2.stateLabel = {'Quiet', 'Whisk'};
end

function values = blankSelfColumn(values)
% The zero-lag summary vectors additionally drop the near-zero-distance
% self-pair column, as the source does before plotting or fitting.
values(:, 1) = nan;
end

% =========================================================================
% Panel A: subthreshold voltage and spiking
% =========================================================================

function panelA = plotVoltageAndSpikes(session, outputDirectory)
% fig4_Pspk_V_def_VShuffle_GLMShuffle.m lines 175-215: the whisking trace with
% quiet and whisking state bars (row 1), the population subthreshold voltage
% (rows 2-3), and the spike raster (rows 4-5). The state bars use the stored
% quiet_mask and whisking_mask, the source's tQuiet and tWhisk.
dt = session.framePeriodSeconds;
frameRange = [round(160 / dt) round(170 / dt)];
windowFrames = frameRange(1):frameRange(2);
cellIndex = session.correlationUnitIndex;
nCells = numel(cellIndex);

% Read all units over the short window, then keep the analysis cells. The
% analysis cells are not contiguous in the units table, and the window is only
% a few thousand frames, so one contiguous read is cheaper than 247 reads.
block = session.readVoltageRows('subthreshold', 1:session.numberOfCells, frameRange);
subthreshold = double(block(cellIndex, :));
timeSeconds = windowFrames * dt;

spikeTimes = cell(nCells, 1);
for iCell = 1:nCells
    frames = session.spikeFrameIndex{cellIndex(iCell)};
    spikeTimes{iCell} = frames(frames >= frameRange(1) & frames <= frameRange(2));
end

whiskInterp = whiskingOnVoltageClock(session);
% Frames outside each state are NaN, so each bar is drawn only over its state.
quietSeconds = timeSeconds;
quietSeconds(~session.quietMask(windowFrames)) = nan;
whiskSeconds = timeSeconds;
whiskSeconds(~session.whiskingMask(windowFrames)) = nan;
ybar = 5;

h = figure('Name', 'Figure 4A', 'Color', 'w', 'Position', [100 100 900 850]);
ax = gobjects(3, 1);

ax(1) = subplot(5, 1, 1);
plot(timeSeconds, whiskInterp(windowFrames));
hold on
plot(quietSeconds, zeros(1, numel(windowFrames)) + ybar, 'color', ones(1, 3) * .5, ...
    'linewidth', 4);
plot(whiskSeconds, zeros(1, numel(windowFrames)) + ybar, 'k', 'linewidth', 4);
ylabel('Whisking');
title('Whisking motion; bars: quiet (gray), whisking (black)');

ax(2) = subplot(5, 1, 2:3);
imagesc(timeSeconds, 1:nCells, subthreshold);
caxis([-1.25 1.25]);
colormap(gca, colorcet('d1'));
cb = colorbar('location', 'east');
cb.Label.String = 'Voltage (norm)';
ylim([0 nCells + 1]);
ylabel('Cell');
title('Subthreshold voltage');

ax(3) = subplot(5, 1, 4:5);
raster_plot(spikeTimes, dt, 'k', timeSeconds([1 end]));
ylim([0 nCells + 1]);
ylabel('Cell');
title('Spikes');

linkaxes(ax, 'x');
xlim(ax(1), timeSeconds([1 end]));
set(ax, 'ydir', 'normal');
% The source creates the whisking axis last, so its xlabel(ax(end), ...) lands
% on the top row. The label is put on the bottom row here instead.
xlabel(ax(end), 'Time (s)');
saveFigure(h, outputDirectory, 'figure4A_voltage_and_spikes_raw.png');

panelA = struct();
panelA.frameRange = frameRange;
panelA.timeRangeSeconds = timeSeconds([1 end]);
panelA.numberOfCells = nCells;
panelA.spikeCountInWindow = cellfun(@numel, spikeTimes);
panelA.whiskingMotionInWindow = whiskInterp(windowFrames);
panelA.quietFrameFraction = mean(session.quietMask(windowFrames));
panelA.whiskingFrameFraction = mean(session.whiskingMask(windowFrames));
panelA.framesPerChunk = session.framesPerChunk;
panelA.stateBarLevel = ybar;
end

function whiskInterp = whiskingOnVoltageClock(session)
% Source lines 66-68. The face-camera whisking trace, without the 8-sample
% moving median used for the state masks, is interpolated onto the voltage
% timebase chunk by chunk with bilinear imresize.
nChunks = session.numberOfChunks;
faceFramesPerChunk = session.numberOfFaceFrames / nChunks;
assert(faceFramesPerChunk == round(faceFramesPerChunk), ...
    'Figure4NWB:FaceChunkStructure', ...
    'The face-video frame count is not a whole multiple of the chunk count.');
% The source hard-codes nFrame = 288000; the loader derives the chunk length
% from recording_chunks. Both behavior sessions give 288000.
if session.framesPerChunk ~= 288000
    warning('Figure4NWB:ChunkFrames', ...
        ['recording_chunks gives %d voltage frames per chunk, but the source ', ...
         'hard-codes nFrame = 288000.'], session.framesPerChunk);
end

byChunk = reshape(session.whiskingMotion, faceFramesPerChunk, nChunks);
interpolated = imresize(byChunk, [session.framesPerChunk nChunks], 'bilinear');
whiskInterp = reshape(interpolated, 1, []);
assert(numel(whiskInterp) == session.numberOfVoltageFrames, ...
    'Figure4NWB:WhiskInterpolation', ...
    'Interpolated whisking trace does not match the voltage frame count.');
end

% =========================================================================
% Panels B, C, I: distance-binned state maps
% =========================================================================

function plotDistanceStateMaps(stage2, seriesField, traceColumns, valueScale, ...
    xLabelText, yLabelText, outputDirectory, fileName)
% fig4_xcorr.m lines 288-320 (B), 338-370 (I), and 390-419 (C). VALUESCALE
% multiplies the plotted traces, image, and colour ceiling: 1 for B and I, 1/dt
% for C, which the source now shows in Hz. The colour floor is taken from the
% unscaled array, as in the source; for C it is 0 either way, because the
% zero-lag self-pair rate was zeroed.
averaged = stage2.averaged.(seriesField);
lagSeconds = stage2.lagSeconds;
centersUm = stage2.bins.rCenterMicrometers;
nDistBin = stage2.bins.numberOfBins;
zeroLag = stage2.zeroLagIndex;

colorMap = max(jet(nDistBin) - .2, 0);
traceColors = colorMap(traceColumns, :);
colorFloor = min(averaged(:));

h = figure('Name', ['Figure 4 ' seriesField], 'Color', 'w', ...
    'Position', [100 100 700 900]);
for state = 1:2
    colorCeiling = averaged(zeroLag, 3, state) * valueScale;

    subplot(4, 1, (state - 1) * 2 + 1);
    plot(lagSeconds, squeeze(averaged(:, traceColumns, state)) * valueScale);
    colororder(gca, traceColors);
    xlabel(xLabelText);
    ylabel(yLabelText);
    xlim(lagSeconds([1 end]));
    title(sprintf('%s, %d distance traces', stage2.stateLabel{state}, ...
        numel(traceColumns)));

    subplot(4, 1, (state - 1) * 2 + 2);
    imagesc(lagSeconds, centersUm, squeeze(averaged(:, :, state))' * valueScale);
    cb = colorbar('location', 'east');
    cb.Label.String = yLabelText;
    if strcmp(seriesField, 'xcorrV')
        caxis([0 colorCeiling]);
    else
        caxis([colorFloor colorCeiling]);
    end
    xlabel(xLabelText);
    ylabel('Distance (um)');
    xlim(lagSeconds([1 end]));
end
saveFigure(h, outputDirectory, fileName);
end

% =========================================================================
% Panel D: distance decay at three lags
% =========================================================================

function panelD = plotDistanceDecayByLag(stage2, outputDirectory)
% fig4_xcorr.m lines 443-474. The right axis shows the spike-rate
% cross-correlogram in Hz, that is divided by dt.
lagOffsets = [0 2 10];
zeroLag = stage2.zeroLagIndex;
centersUm = stage2.bins.rCenterMicrometers;
dt = stage2.framePeriodSeconds;

lineColors = colorcet('bmy');
lineColors = lineColors(round(linspace(1, 100, numel(lagOffsets))), :);

h = figure('Name', 'Figure 4D', 'Color', 'w', 'Position', [100 100 700 500]);
voltageByLag = zeros(numel(lagOffsets), stage2.bins.numberOfBins);
spikeByLag = zeros(numel(lagOffsets), stage2.bins.numberOfBins);
for iLag = 1:numel(lagOffsets)
    voltageByLag(iLag, :) = squeeze(stage2.averaged.xcorrV(zeroLag + lagOffsets(iLag), :, 1));
    spikeByLag(iLag, :) = squeeze(stage2.averaged.staRate(zeroLag + lagOffsets(iLag), :, 1));

    yyaxis left
    plot(centersUm, voltageByLag(iLag, :), 'color', lineColors(iLag, :), 'linestyle', '-');
    hold on

    yyaxis right
    plot(centersUm, spikeByLag(iLag, :) / dt, 'color', lineColors(iLag, :), 'linestyle', '-');
    hold on
end
xlabel('Distance (um)');
yyaxis left
ylabel('V Correlation');
ylim([0.2 0.7]);
yyaxis right
ylabel('P(post spk | pre spk) (Hz)');
ylim([1 5]);
legend(arrayfun(@(x) sprintf('lag = %.2f ms', x), lagOffsets * dt * 1e3, ...
    'uniformoutput', false), 'location', 'northeast');
saveFigure(h, outputDirectory, 'figure4D_distance_decay_by_lag_raw.png');

panelD = struct();
panelD.lagFrames = lagOffsets;
panelD.lagMilliseconds = lagOffsets * dt * 1e3;
panelD.voltageCorrelationByLag = voltageByLag;
panelD.spikeCorrelationByLag = spikeByLag;
panelD.spikeCorrelationByLagHz = spikeByLag / dt;
panelD.distanceMicrometers = centersUm;
end

% =========================================================================
% Panels E, F, G: subthreshold voltage statistics
% =========================================================================

function voltage = accumulateVoltageStatistics(session, needJoint)
% Streams the subthreshold voltage once and accumulates every quantity panels
% E, F, and G need. All three are histogram accumulations over time, so summing
% per-block contributions is exact.
%
% The read must be frame blocks across all cells: the HDF5 chunks span all
% cells, so a per-cell loop would re-read the whole dataset once per cell. See
% the performance note in loadFigure4SessionFromNWB.
masks = buildTimeMasks(session);
cellIndex = session.correlationUnitIndex;
nCells = numel(cellIndex);
nUnits = session.numberOfCells;
nFramesTotal = session.numberOfVoltageFrames;

vEdge = linspace(-2, 2, 101);
nV = numel(vEdge) - 1;

vHist = zeros(nCells, nV);
vSpkHist = zeros(nCells, nV);
N = zeros(nV, nV);

% Restrict each cell's spike frames to the analyzed frames once, up front. The
% source writes intersect(spk_t{ii}, idxMask); source spike frames are unique
% and sorted, so filtering is equivalent.
spikeFrames = cell(nCells, 1);
for iCell = 1:nCells
    frames = session.spikeFrameIndex{cellIndex(iCell)};
    spikeFrames{iCell} = frames(masks.analysisEF(frames));
end

% Panel G neighbor sets: 60-200 um. The 60 um floor is in the source but not in
% the caption.
rLimUm = 200;
distanceUm = session.correlationDistancePixels * session.pixelSizeMicrometers;
neighborOf = cell(nCells, 1);
neighborCount = zeros(nCells, 1);
for iCell = 1:nCells
    neighborOf{iCell} = find(distanceUm(iCell, :) > 60 & distanceUm(iCell, :) < rLimUm);
    neighborCount(iCell) = numel(neighborOf{iCell});
end

blockFrames = max(session.voltageChunkFrames, 1) * 20;
blockStarts = 1:blockFrames:nFramesTotal;
for iBlock = 1:numel(blockStarts)
    firstFrame = blockStarts(iBlock);
    lastFrame = min(firstFrame + blockFrames - 1, nFramesTotal);
    keepEF = masks.analysisEF(firstFrame:lastFrame);
    keepG = masks.analysisG(firstFrame:lastFrame);
    if ~any(keepEF)
        continue
    end

    block = session.readVoltageRows('subthreshold', 1:nUnits, [firstFrame lastFrame]);
    block = double(block(cellIndex, :));

    % ---- panels E and F
    for iCell = 1:nCells
        vHist(iCell, :) = vHist(iCell, :) + ...
            histcounts(block(iCell, keepEF), vEdge);

        frames = spikeFrames{iCell};
        frames = frames(frames >= firstFrame & frames <= lastFrame);
        if ~isempty(frames)
            vSpkHist(iCell, :) = vSpkHist(iCell, :) + ...
                histcounts(block(iCell, frames - firstFrame + 1), vEdge);
        end
    end

    % ---- panel G
    if needJoint && any(keepG)
        gBlock = block(:, keepG);
        ownBin = discretize(gBlock, vEdge);
        for iCell = 1:nCells
            if neighborCount(iCell) == 0
                continue
            end
            own = ownBin(iCell, :);
            valid = ~isnan(own);
            if ~any(valid)
                continue
            end
            neighborBin = ownBin(neighborOf{iCell}, valid);
            ownValid = repmat(own(valid), neighborCount(iCell), 1);
            pairOk = ~isnan(neighborBin);
            if ~any(pairOk(:))
                continue
            end
            % Equivalent to the source's per-own-bin histcounts of the
            % neighbors' voltage, accumulated as a 2-D histogram. The (:)
            % calls matter: with a single neighbor these are row vectors, and
            % accumarray needs one column per subscript dimension.
            neighborSubscript = neighborBin(pairOk);
            ownSubscript = ownValid(pairOk);
            N = N + accumarray([neighborSubscript(:), ownSubscript(:)], 1, [nV nV]);
        end
    end
end

voltage = struct();
voltage.cellIndex = cellIndex;
voltage.numberOfCells = nCells;
voltage.masks = masks;
voltage.voltageEdges = vEdge;
voltage.voltageCenters = mean(vEdge([1:end-1; 2:end]), 1);
voltage.allTimeCounts = vHist;
voltage.spikeTimeCounts = vSpkHist;
voltage.jointCounts = N;
voltage.neighborCountPerCell = neighborCount;
voltage.neighborRangeMicrometers = [60 rLimUm];
voltage.frameCountEF = sum(masks.analysisEF);
voltage.frameCountG = sum(masks.analysisG);
voltage.framePeriodSeconds = session.framePeriodSeconds;
end

function histograms = voltageHistograms(voltage)
% fig4_Pspk_V_def_VShuffle_GLMShuffle.m lines 419-465: the distribution of
% subthreshold voltage at all analyzed times, and at spike times only.
% V(t_s) is the single subthreshold sample at each spike frame, as in the
% source's intersect(spk_t{ii}, idxMask). The caption describes it as the
% "10 ms interval prior to spikes"; the author confirmed on 2026-09-26 that
% the wording stands and the computation is unchanged.
histograms = struct();
histograms.voltageEdges = voltage.voltageEdges;
histograms.voltageCenters = voltage.voltageCenters;
histograms.allTimeCounts = voltage.allTimeCounts;
histograms.spikeTimeCounts = voltage.spikeTimeCounts;
histograms.numberOfCells = voltage.numberOfCells;
histograms.analyzedFrameCount = voltage.frameCountEF;
histograms.framePeriodSeconds = voltage.framePeriodSeconds;
end

function plotVoltageHistograms(histograms, outputDirectory)
h = figure('Name', 'Figure 4E', 'Color', 'w', 'Position', [100 100 600 700]);
ax = gobjects(2, 1);

ax(1) = subplot(2, 1, 1);
plot(histograms.voltageCenters, histograms.allTimeCounts');
title('V');
ylabel('Counts');

ax(2) = subplot(2, 1, 2);
plot(histograms.voltageCenters, histograms.spikeTimeCounts');
title('V(t_s)');
xlabel('V_{sub} (SH)');
ylabel('Counts');

set(ax, 'fontsize', 12);
saveFigure(h, outputDirectory, 'figure4E_voltage_histograms_raw.png');
end

function panelF = plotExcitabilityCurves(histograms, session, outputDirectory)
% The excitability curve is the ratio of the two panel E distributions,
% converted to Hz by dividing by the frame period.
dt = histograms.framePeriodSeconds;
nCells = histograms.numberOfCells;
rate = histograms.spikeTimeCounts ./ histograms.allTimeCounts / dt;

nHighlight = 7;
rng(0, 'twister');
highlightIndex = randperm(nCells, nHighlight);

h = figure('Name', 'Figure 4F', 'Color', 'w', 'Position', [100 100 650 550]);
semilogy(histograms.voltageCenters, rate', 'color', ones(1, 3) * .5);
hold on
semilogy(histograms.voltageCenters, rate(highlightIndex, :)', 'linewidth', 2);
xlabel('V_{sub} (SH)');
ylabel('Hz');
yticks(10 .^ (-2:2:2));
xlim([-.5 2]);
ylim(10 .^ [-2 2.5]);
set(gca, 'fontsize', 12);
saveFigure(h, outputDirectory, 'figure4F_excitability_curves_raw.png');

panelF = struct();
panelF.spikeRateGivenVoltageHz = rate;
panelF.voltageCenters = histograms.voltageCenters;
panelF.highlightCorrelationIndex = highlightIndex;
panelF.highlightUnitId = session.cellId(session.correlationUnitIndex(highlightIndex));
panelF.randomSeed = 'rng(0,''twister'') then randperm(nCells,7)';
end

function joint = neighborVoltageJoint(voltage)
% fig4_Pspk_V_def_VShuffle_GLMShuffle.m lines 842-876. For every analysis cell,
% histogram its neighbors' subthreshold voltage conditioned on the cell's own
% voltage bin. The accumulation happens in accumulateVoltageStatistics so that
% the recording is streamed only once.
joint = struct();
joint.voltageEdges = voltage.voltageEdges;
joint.voltageCenters = voltage.voltageCenters;
joint.counts = voltage.jointCounts;
joint.numberOfCells = voltage.numberOfCells;
joint.neighborCountPerCell = voltage.neighborCountPerCell;
joint.neighborRangeMicrometers = voltage.neighborRangeMicrometers;
joint.analyzedFrameCount = voltage.frameCountG;
end

function panelG = combineNeighborJoint(joint27, joint29)
% The published panel sums the count matrix over both mice before normalizing.
panelG = struct();
panelG.voltageCenters = joint27.voltageCenters;
panelG.voltageEdges = joint27.voltageEdges;
panelG.neighborRangeMicrometers = joint27.neighborRangeMicrometers;
panelG.counts27 = joint27.counts;
panelG.numberOfCells27 = joint27.numberOfCells;
if nargin < 2 || isempty(joint29)
    panelG.counts = joint27.counts;
    panelG.numberOfCells = joint27.numberOfCells;
    panelG.numberOfMice = 1;
else
    assert(isequal(joint27.voltageEdges, joint29.voltageEdges), ...
        'Figure4NWB:JointEdges', 'The two sessions used different voltage bins.');
    panelG.counts29 = joint29.counts;
    panelG.numberOfCells29 = joint29.numberOfCells;
    panelG.counts = joint27.counts + joint29.counts;
    panelG.numberOfCells = joint27.numberOfCells + joint29.numberOfCells;
    panelG.numberOfMice = 2;
end

% N is symmetric by construction, because every neighbor pair is accumulated in
% both orderings. Normalizing each row and transposing therefore makes each
% displayed column sum to 1, so the x axis is the conditioning variable.
panelG.conditionalProbability = (panelG.counts ./ sum(panelG.counts, 2))';
panelG.countMatrixSymmetryError = max(abs(panelG.counts - panelG.counts'), [], 'all') ...
    / max(panelG.counts(:));
end

function plotNeighborVoltageJoint(panelG, outputDirectory)
h = figure('Name', 'Figure 4G', 'Color', 'w', 'Position', [100 100 620 560]);
imagesc(panelG.voltageCenters, panelG.voltageCenters, panelG.conditionalProbability);
hold on
plot([-2 2], [-2 2], 'r', 'linewidth', 1);

% Whiten the lowest probabilities so the occupied region stands out, as the
% source does with a graded ramp over the first 50 colormap entries.
cMap = parula(2 ^ 12);
nGrad = 50;
cMap(1:nGrad, :) = interp1([1 nGrad], [ones(1, 3); cMap(nGrad, :)], ...
    linspace(1, nGrad, nGrad));
colormap(cMap);
cb = colorbar;
cb.Label.String = 'Probability';
daspect([1 1 1]);
xlabel('V_i (SH)');
ylabel('V_{neighbor} (SH)');
title(sprintf('n = %d cells, N = %d mice', panelG.numberOfCells, panelG.numberOfMice));
set(gca, 'ydir', 'normal', 'fontsize', 12);
saveFigure(h, outputDirectory, 'figure4G_neighbor_voltage_conditional_raw.png');
end

% =========================================================================
% Panel H: example spike-triggered averages
% =========================================================================

function panelH = plotStaExamples(session, outputDirectory)
% fig4_xcorr.m lines 543-577. The two triggering cells are indices into the
% correlation cell axis, not unit IDs.
triggerIndex = [209 126];
nNeighbor = 5;
rLimUm = 200;
yOffsetStep = [.2 .2];
publishedDistancesUm = {[88 115 130 183 198], [134 136 136 177 177]};

dx = session.pixelSizeMicrometers;
distanceUm = session.correlationDistancePixels * dx;
lagSeconds = session.lagSeconds;
assert(all(triggerIndex <= session.correlationCellCount), 'Figure4NWB:TriggerIndex', ...
    'The example trigger indices exceed the %d correlation cells in this session.', ...
    session.correlationCellCount);

h = figure('Name', 'Figure 4H', 'Color', 'w', 'Position', [100 100 900 500]);
ax = gobjects(numel(triggerIndex), 1);
panelH = struct();
panelH.triggerCorrelationIndex = triggerIndex;
panelH.triggerUnitId = session.cellId(session.correlationUnitIndex(triggerIndex));
panelH.publishedDistancesUm = publishedDistancesUm;
panelH.neighborDistancesUm = cell(numel(triggerIndex), 1);
panelH.neighborUnitId = cell(numel(triggerIndex), 1);
panelH.traces = cell(numel(triggerIndex), 1);

for iExample = 1:numel(triggerIndex)
    trigger = triggerIndex(iExample);
    candidates = find(distanceUm(trigger, :) < rLimUm & distanceUm(trigger, :) > 60);
    [~, order] = sort(distanceUm(trigger, candidates));
    candidates = candidates(order);
    assert(~isempty(candidates), 'Figure4NWB:NoNeighbors', ...
        'Correlation cell %d has no neighbor between 60 and %d um.', trigger, rLimUm);
    selected = candidates(round(linspace(1, numel(candidates), nNeighbor)));

    staRows = session.readCorrelationTrigger('spike_triggered_voltage', 'quiet', trigger);
    traces = staRows([trigger selected], :);
    traces = traces - mean(traces, 2);
    yOffsets = (0:size(traces, 1) - 1)' * yOffsetStep(iExample);

    ax(iExample) = subplot(1, numel(triggerIndex), iExample);
    plot(lagSeconds, (traces - yOffsets)');
    colororder(gca, winter(nNeighbor + 1));
    text(zeros(nNeighbor, 1) + .1, -yOffsets(2:end), ...
        arrayfun(@(x) sprintf('%.0f um', x), distanceUm(trigger, selected), ...
        'uniformoutput', false));
    xlabel('Peri-spike time (s)');
    ylabel('Voltage (norm)');
    title(sprintf('Trigger correlation index %d', trigger));

    panelH.neighborDistancesUm{iExample} = distanceUm(trigger, selected);
    panelH.neighborUnitId{iExample} = ...
        session.cellId(session.correlationUnitIndex(selected));
    panelH.traces{iExample} = traces;
end
linkaxes(ax, 'y');
saveFigure(h, outputDirectory, 'figure4H_sta_examples_raw.png');

% The published panel prints these distances as its trace labels, so they are a
% direct check that the example selection still resolves to the same cells.
panelH.matchesPublishedLabels = true;
for iExample = 1:numel(triggerIndex)
    reproduced = round(panelH.neighborDistancesUm{iExample});
    expected = publishedDistancesUm{iExample};
    if ~isequal(reproduced, expected)
        panelH.matchesPublishedLabels = false;
        warning('Figure4NWB:PanelHLabels', ...
            ['Panel H example %d gave distances %s um but the published panel ', ...
             'is labelled %s um. The example selection may no longer resolve to ', ...
             'the same cells.'], iExample, mat2str(reproduced), mat2str(expected));
    end
end
if panelH.matchesPublishedLabels
    fprintf('Panel H: both examples reproduce the published distance labels.\n');
end
end

% =========================================================================
% Reported statistics
% =========================================================================

function statistics = reportedStatistics(stage2_27, stage2_29)
% fig4_xcorr.m lines 656-726 (cross-mouse short-range correlation and
% exponential decay length constants), plus the per-session widths at lines
% 479-530 and the cross-STA depolarization size at lines 578-582. The source
% runs the per-session blocks in one session folder at a time, so those are
% reported per mouse.
centersUm = stage2_27.bins.rCenterMicrometers;
centersPixels = stage2_27.bins.rCenterPixels;
dx = stage2_27.bins.pixelSizeMicrometers;
stage2All = {stage2_27, stage2_29};
summaries = {stage2_27.summary, stage2_29.summary};
sessionNames = {'M-YQ0201-27', 'M-YQ0201-29'};
nSessions = numel(summaries);

statistics = struct();
statistics.distanceMicrometers = centersUm;
statistics.sessionNames = sessionNames;

% Short-range voltage correlation. Bin 3 spans 60-120 um; the source labels it
% with the midpoints of the neighboring bin centers, which is where the
% manuscript's "60 um < r < 120 um" comes from.
shortRange = zeros(2, nSessions);
for iSession = 1:nSessions
    shortRange(:, iSession) = summaries{iSession}.xcorrVT0Val(:, 3);
end
statistics.shortRangeBin = 3;
statistics.shortRangeWindowUm = [mean(centersUm([2 3])) mean(centersUm([3 4]))];
statistics.shortRangeVoltageCorrelationPerSession = shortRange;
statistics.shortRangeVoltageCorrelationMean = mean(shortRange, 2);

% Length constants. The voltage and 2.5 ms-lag spike fits use every non-blank
% bin; only the synchrony fit is restricted to r < 800 um (source lines
% 668-672, 689-697, and 715-719). The 2.5 ms-lag fit takes its bin mask from
% staSync rather than staFRT2Val, as the source does; both have the same blank
% columns.
fits = {
    'voltageCrossCorrelation', 'xcorrVT0Val', 'xcorrVT0Val', inf
    'spikeCorrelationLag2',    'staFRT2Val',  'staSync',     inf
    'spikeSynchrony',          'staSync',     'staSync',     800 / dx
    };
for iFit = 1:size(fits, 1)
    field = fits{iFit, 1};
    source = fits{iFit, 2};
    maskSource = fits{iFit, 3};
    centerLimit = fits{iFit, 4};
    perSession = zeros(nSessions, 2);
    for iSession = 1:nSessions
        values = summaries{iSession}.(source);
        maskValues = summaries{iSession}.(maskSource);
        for state = 1:2
            fitMask = ~isnan(maskValues(state, :)) & centersPixels < centerLimit;
            selected = values(state, fitMask);
            assert(all(selected > 0), 'Figure4NWB:LogFit', ...
                ['%s contains a non-positive value for %s in state %d, so the ', ...
                 'source log-space fit cannot be reproduced.'], source, ...
                sessionNames{iSession}, state);
            coefficients = polyfit(centersUm(fitMask), log(selected), 1);
            perSession(iSession, state) = 1 / coefficients(1);
        end
    end
    statistics.(field) = struct( ...
        'lengthConstantPerSessionUm', perSession, ...
        'lengthConstantMeanUm', mean(perSession, 1), ...
        'centerLimitUm', centerLimit * dx, ...
        'stateOrder', {{'quiet', 'whisk'}});
end

% Correlation widths and depolarization size, per session.
voltageFwhmMs = zeros(nSessions, numel(centersUm));
spikeHalfWidthMs = zeros(nSessions, 1);
synchronyFwhmMs = zeros(nSessions, 2);
synchronyHalfWidthFrames = zeros(nSessions, 2, 2);
depolarization = zeros(nSessions, 1);
for iSession = 1:nSessions
    widths = correlationWidths(stage2All{iSession});
    voltageFwhmMs(iSession, :) = widths.voltageFwhmSeconds * 1e3;
    spikeHalfWidthMs(iSession) = widths.spikeHalfWidthSeconds * 1e3;
    synchronyFwhmMs(iSession, :) = widths.synchronyFwhmSeconds' * 1e3;
    synchronyHalfWidthFrames(iSession, :, :) = widths.synchronyHalfWidthFrames;
    depolarization(iSession) = stage2All{iSession}.crossStaDepolarizationQuiet;
end
statistics.voltageCrossCorrelationFwhm = struct( ...
    'perDistanceBinMs', voltageFwhmMs, ...
    'medianPerSessionMs', median(voltageFwhmMs, 2), ...
    'medianMeanMs', mean(median(voltageFwhmMs, 2)), ...
    'state', 'quiet', ...
    'sourceLines', 'fig4_xcorr.m 479-492');
statistics.spikeRateHalfWidth = struct( ...
    'perSessionMs', spikeHalfWidthMs, ...
    'meanMs', mean(spikeHalfWidthMs), ...
    'distanceBin', 3, ...
    'distanceRangeUm', stage2_27.bins.rEdgePixels([3 4]) * dx, ...
    'state', 'whisk', ...
    'sourceLines', 'fig4_xcorr.m 494-504');
statistics.synchronyFwhm = struct( ...
    'perSessionMs', synchronyFwhmMs, ...
    'meanMs', mean(synchronyFwhmMs, 1), ...
    'halfWidthFramesBeforeAfter', synchronyHalfWidthFrames, ...
    'distanceCentersUm', [60 400], ...
    'stateOrder', {{'quiet', 'whisk'}}, ...
    'sourceLines', 'fig4_xcorr.m 514-530');
statistics.crossStaDepolarization = struct( ...
    'perSession', depolarization, ...
    'mean', mean(depolarization), ...
    'units', 'fraction of the target cell spike height', ...
    'pairDistanceUm', [60 400], ...
    'lagFrames', [-100 -2], ...
    'state', 'quiet', ...
    'sourceLines', 'fig4_xcorr.m 578-582');

fprintf('\nFigure 4 reported statistics\n');
fprintf('  short-range V correlation (%.1f-%.1f um), mean of %d mice: %.3f quiet, %.3f whisk\n', ...
    statistics.shortRangeWindowUm, nSessions, statistics.shortRangeVoltageCorrelationMean);
fprintf('  V xcorr length constant, mean: %.1f um quiet, %.1f um whisk\n', ...
    statistics.voltageCrossCorrelation.lengthConstantMeanUm);
fprintf('  spike xcorr at 2.5 ms lag, mean: %.1f um quiet, %.1f um whisk\n', ...
    statistics.spikeCorrelationLag2.lengthConstantMeanUm);
for iSession = 1:nSessions
    fprintf('    %s alone: %.1f um quiet, %.1f um whisk\n', sessionNames{iSession}, ...
        statistics.spikeCorrelationLag2.lengthConstantPerSessionUm(iSession, :));
end
fprintf('  spike synchrony length constant, mean: %.1f um quiet, %.1f um whisk\n', ...
    statistics.spikeSynchrony.lengthConstantMeanUm);
for iSession = 1:nSessions
    fprintf(['  %s: quiet V xcorr FWHM %.1f ms (median over bins); ', ...
        'whisk spike xcorr half-width %.2f ms (60-120 um); ', ...
        'synchrony FWHM %.2f ms quiet, %.2f ms whisk; ', ...
        'cross-STA depolarization %.3f spike heights\n'], sessionNames{iSession}, ...
        statistics.voltageCrossCorrelationFwhm.medianPerSessionMs(iSession), ...
        spikeHalfWidthMs(iSession), synchronyFwhmMs(iSession, :), depolarization(iSession));
end
end

function widths = correlationWidths(stage2)
% Each block is copied from its own fig4_xcorr.m section.
averaged = stage2.averaged;
dt = stage2.framePeriodSeconds;
nLag = stage2.zeroLagIndex - 1;
centersUm = stage2.bins.rCenterMicrometers;
widths = struct();

% Voltage cross-correlation FWHM, lines 479-492: quiet state, every distance
% column. The half level lies between the value at the edge of the lag window
% (-nLag) and the zero-lag peak, not between zero and the peak. Column 2 is all
% NaN, so min() returns index 1 there and that column contributes the full
% 2*nLag*dt window to the median. The source keeps that column, and so does
% this.
tr = squeeze(averaged.xcorrV(1:nLag + 1, :, 1));
hm = tr(1, :) + (tr(nLag + 1, :) - tr(1, :)) / 2;
[~, halfIndex] = min(abs(tr - hm), [], 1);
widths.voltageFwhmSeconds = (nLag + 1 - halfIndex) * dt * 2;

% Spike-train cross-correlogram half-width, lines 494-504: 60-120 um bin,
% whisking state, the same edge-relative half level. Not doubled, so this is a
% half-width, not a full width.
tr = squeeze(averaged.staRate(1:nLag + 1, 3, 2));
hm = tr(1, :) + (tr(nLag + 1, :) - tr(1, :)) / 2;
[~, halfIndex] = min(abs(tr - hm));
widths.spikeHalfWidthSeconds = median((nLag + 1 - halfIndex) * dt);

% Synchrony FWHM, lines 514-530: the baseline-subtracted spike-rate
% cross-correlogram averaged over bins centered in 60-400 um, with the half
% level at half the zero-lag value, interpolated separately on each side.
idxUse = centersUm > 60 & centersUm <= 400;
tr = squeeze(mean(averaged.staRate(:, idxUse, :) - ...
    max(cat(1, mean(averaged.staRate(nLag + 1 + (-5:-3), idxUse, :), 1), ...
    mean(averaged.staRate(nLag + 1 + (3:5), idxUse, :), 1)), [], 1), 2));
halfWidthFrames = zeros(2, 2);
for state = 1:2
    halfWidthFrames(state, 1) = nLag + 1 - ...
        interp1(tr(1:nLag + 1, state), (1:nLag + 1)', tr(nLag + 1, state) / 2);
    halfWidthFrames(state, 2) = ...
        interp1(tr(nLag + 1:end, state), (1:nLag + 1)', tr(nLag + 1, state) / 2) - 1;
end
widths.synchronyHalfWidthFrames = halfWidthFrames;
widths.synchronyFwhmSeconds = sum(halfWidthFrames, 2) * dt;
end

% =========================================================================

function saveFigure(figureHandle, outputDirectory, fileName)
set(figureHandle, 'PaperPositionMode', 'auto');
target = fullfile(outputDirectory, fileName);
writeWithRetry(@() print(figureHandle, target, '-dpng', '-r300'), target);
set(figureHandle, 'Visible', 'on');
drawnow;
end

function saveResults(resultsPath, results)
% save() resolves variable names in the calling workspace, so the results struct
% has to be a named local here rather than captured by an anonymous function.
save(resultsPath, 'results', '-v7.3');
end

function writeWithRetry(writeFcn, target)
% The default output directory lives inside a Dropbox folder. Dropbox locks a
% file for a moment while it syncs, and it dehydrates older files into
% online-only placeholders that must be rehydrated before they can be
% overwritten. Either can make an otherwise valid write fail with
% "PNG library failed: Could not open file" or a similar open error, and the
% failure is transient. Retry briefly, then explain the likely cause instead of
% surfacing the bare library error.
maxAttempts = 5;
pauseSeconds = 2;
for attempt = 1:maxAttempts
    try
        writeFcn();
        if attempt > 1
            fprintf('Wrote %s on attempt %d.\n', target, attempt);
        end
        return
    catch ME
        if attempt == maxAttempts
            error('Figure4NWB:WriteFailed', ...
                ['Could not write this file after %d attempts over %d s:\n' ...
                 '  %s\n' ...
                 'Last error: %s\n\n' ...
                 'This path is inside a Dropbox folder. Dropbox locks files ' ...
                 'while syncing and stores older ones as online-only ' ...
                 'placeholders. Either pause Dropbox syncing and re-run, or ' ...
                 'send the output somewhere local, for example:\n' ...
                 '  reproduce_figure4_from_nwb([], [], ''D:\\scratch\\fig4'')'], ...
                maxAttempts, maxAttempts * pauseSeconds, target, ME.message);
        end
        pause(pauseSeconds);
    end
end
end
