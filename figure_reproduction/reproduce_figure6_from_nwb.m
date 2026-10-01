function results = reproduce_figure6_from_nwb(behaviorNwbFiles, optoNwbFile, ...
    connNwbFile, outputDirectory, panels)
% REPRODUCE_FIGURE6_FROM_NWB Recreate the data panels of manuscript Figure 6.
%
%   RESULTS = REPRODUCE_FIGURE6_FROM_NWB(BEHAVIORNWBFILES, OPTONWBFILE, CONNNWBFILE, OUTPUTDIRECTORY, PANELS)
%   reads the four NWB files Figure 6 needs, writes raw panels and a reference
%   MAT file to OUTPUTDIRECTORY, and returns the numerical results. It does not
%   edit or execute the legacy scripts.
%
%   Figure 6 spans two mice and three protocols. Every assignment below was
%   confirmed numerically during the Phase 1 audit; see FIGURE6_PHASE1_HANDOFF.md.
%
%       B, C              M_YQ0201_29_sparseOpto.nwb       intrinsic properties
%       D                 M_YQ0201_27_sparsePulseHad.nwb   connectivity
%       E, F, H, I, J, K  M_YQ0201_29_behavior.nwb         spontaneous activity
%       G                 both behavior files              pooled, 2 mice
%
%   Panel A is the spinning-disk confocal image. It is not MATLAB-reproducible
%   from NWB (its source TIFFs are not deposited) and the author ruled that it
%   need not be reproduced, so it is not a valid selector.
%
%   All five arguments are optional and follow the same convention as
%   reproduce_figure2/3/4/5_from_nwb.m: an omitted or empty argument falls back
%   to its default, so REPRODUCE_FIGURE6_FROM_NWB() regenerates the whole figure
%   from the documented session paths.
%
%     behaviorNwbFiles  cell array; the first entry is the session panels E, F,
%                       H, I, J and K display, and the whole list is what panel
%                       G pools. Default {M-YQ0201-29; M-YQ0201-27}. A single
%                       character vector is accepted and treated as a one-file
%                       list, which makes panel G single-mouse.
%     optoNwbFile       default M_YQ0201_29_sparseOpto.nwb
%     connNwbFile       default M_YQ0201_27_sparsePulseHad.nwb
%     outputDirectory   default figure6_nwb_reproduction/ beside this file
%     panels            default all of B through K
%
%   PANELS accepts a panel letter, a comma-separated character vector, or a cell
%   array, for example 'D', 'I,J', or {'E','F','G'}. Only the NWB sessions the
%   request needs are opened, so 'D' never touches a behavior file and 'B,C'
%   never touches the connectivity file.
%
%   Each panel is drawn in its own named, visible figure window, printed to PNG
%   at 300 dpi, and left open for inspection, as in the other reproduction
%   scripts.
%
%   Examples
%   --------
%     reproduce_figure6_from_nwb();                          % everything
%     reproduce_figure6_from_nwb([], [], [], [], 'D');       % one cheap panel
%     reproduce_figure6_from_nwb([], [], [], [], 'I,J,K');   % behavior panels
%     reproduce_figure6_from_nwb([], [], [], 'C:\tmp\fig6');
%
%   Source provenance
%   -----------------
%   B   fig6_NPYpm_optoAvgTraces.m           lines 227-287  (figure 12)
%   C   fig6_NPYpm_ephys.m                   lines 139-159  (figure 20)
%   D   fig6_NPYpm_conn.m                    lines 160-271  (figure 1)
%   E   fig6_NPYpm_avgDynamics_raster.m      lines 120-230  (figures 1 and 3)
%   F   fig6_NPY_umap_confusionmat.m         lines  96-116  (figure 2)
%   G   fig6_NPY_umap_confusionmat.m         lines 121-161  (figure 99)
%   H   fig6_NPYpm_avgDynamics_raster.m      lines 231-284  (figure 2)
%   I   fig6_NPY_osc_whisk_TrigAvg_ranksumStats.m lines 192-227 (figure 98)
%   J   fig6_NPY_osc_whisk_TrigAvg_ranksumStats.m lines 343-430 (figure 99)
%   K   fig6_NPYpm_pSpkV_whisk.m             lines 354-426  (figure 11)
%   Full panel-level evidence is in FIGURE_CODE_MAP.md.
%
%   Preserved source behavior
%   -------------------------
%   These are copied deliberately and must not be "cleaned up":
%     * Two different quiet-state rules coexist. Panels I, J and K use the
%       Gaussian-mixture masks stored in the NWB. Panel H's state bar uses the
%       source's median-threshold rule, recomputed here from whisking_motion.
%     * The simultaneous-spike rejection is > 15 here, not the > 20 of
%       fig5_ephys_vs_spiking.m, and it tests sum(spk_t_mat) while assigning
%       into spk_t_mat_c.
%     * Panel B builds its voltage from the **unsmoothed** subthreshold trace,
%       linear interpolation across a [-3:5] frame window around every spike,
%       not the sgolay-filtered subthreshold_voltage the other panels use.
%     * Panel D's blue-stimulation correction is a three-stage procedure: a
%       pw_doubleExp template fitted to 60-120 um unconnected pairs, regressed
%       out of every 60-400 um unconnected pair over the first 18 post-stimulus
%       frames only, keeping the slope coefficient and discarding the intercept.
%     * Panel D averages with a plain nanmean over all pairs in the 60-400 um
%       band, connected and unconnected alike, with no spike-count weighting.
%     * Panels E, F and G bin firing rate at 500 ms over the first 2000 s, with
%       a chunk-boundary guard that is +/-100 frames in the raster script and
%       +/-5000 frames in the confusion-matrix script. Both are computed and
%       compared; see results.panelF.guardComparison.
%     * Panel E orders its raster by a Fisher linear discriminant on the 2-D
%       UMAP coordinates, not by the 1-component UMAP the source computes first
%       and then overwrites.
%     * Panel I's rate STA uses the raw spike matrix while panel J's uses the
%       > 15-rejected one. The author confirmed this is immaterial for
%       event-triggered averages; both are reproduced as written.
%     * Panel J masks the voltage average at chunk boundaries but not the rate
%       average.
%     * Panel K's "whisking" state is the 0-635 ms window after a selected
%       whisking onset, not the whole whisking epoch, and its curves pool counts
%       across cells before dividing.
%
%   Documented deviations from the legacy scripts
%   ---------------------------------------------
%   1. Motion-corrected spike times. The source scripts call
%      spike_times_motion_correction_aggressive, which needs mcTrace_all. That
%      variable is not in the NWB and the author waived the step on 2026-09-18,
%      so the stored raw spike times are used. Note that the Figure 6 behavior
%      scripts already have this call commented out.
%   2. Panel B's cell selection. The source adds a 2-sigma spike-height term
%      that /units/passed_qc does not carry. The author waived it on 2026-09-19
%      and it was measured to exclude no cell in M-YQ0201-29 sparseOpto. This
%      script re-runs that check where it can and records the result in
%      results.panelB.selectionCheck.
%   3. UMAP and k-means are unseeded in the source and no rng call exists
%      anywhere in the Figure 6 scripts. Panels E, F and G are therefore not
%      bit-reproducible. The k-means labels are aligned to the NPY subtype
%      averages after clustering so that panel G's columns cannot swap between
%      runs; see alignGroupsToCellTypes. That fixes the orientation of the
%      comparison only, never its strength.
%   4. Third-party colour maps (colorcet) are replaced with built-in colours.
%      Cosmetic only.
%
%   Outputs
%   -------
%   One PNG per panel plus figure6_nwb_reference_results.mat in
%   OUTPUTDIRECTORY, and a RESULTS struct with one field per panel, a
%   sampleSizes field carrying the caption counts, and a deviations field.
%
%   See also LOADFIGURE6SESSIONFROMNWB, REPRODUCE_FIGURE5_FROM_NWB.

% Placeholders: pass the paths to your local copies of these files instead.
% The first behavior file is the session panels E, F, H-K display.
if nargin < 1 || isempty(behaviorNwbFiles)
    behaviorNwbFiles = { ...
        'path/to/M_YQ0201_29_behavior.nwb'; ...
        'path/to/M_YQ0201_27_behavior.nwb'};
end
if ischar(behaviorNwbFiles)
    behaviorNwbFiles = {behaviorNwbFiles};
end
if nargin < 2 || isempty(optoNwbFile)
    optoNwbFile = 'path/to/M_YQ0201_29_sparseOpto.nwb';
end
if nargin < 3 || isempty(connNwbFile)
    connNwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';
end
if nargin < 4 || isempty(outputDirectory)
    outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
        'figure6_nwb_reproduction');
end
if nargin < 5
    panels = [];
end
selected = normalizePanelSelection(panels);
if exist(outputDirectory, 'dir') ~= 7
    mkdir(outputDirectory);
end

needOpto = any(ismember({'B', 'C'}, selected));
needConn = ismember('D', selected);
needDisplayBehavior = any(ismember({'E', 'F', 'H', 'I', 'J', 'K'}, selected));
needAllBehavior = ismember('G', selected);

if needOpto
    assert(exist(optoNwbFile, 'file') == 2, 'Figure6:MissingOptoFile', ...
        'Panels B and C need the sparseOpto NWB file, not found: %s', optoNwbFile);
end
if needConn
    assert(exist(connNwbFile, 'file') == 2, 'Figure6:MissingConnFile', ...
        'Panel D needs the sparsePulseHad NWB file, not found: %s', connNwbFile);
end
if needDisplayBehavior || needAllBehavior
    assert(~isempty(behaviorNwbFiles), 'Figure6:MissingBehaviorFile', ...
        'Panels E to K need at least one behavior NWB file.');
    nBehaviorNeeded = 1;
    if needAllBehavior
        nBehaviorNeeded = numel(behaviorNwbFiles);
    end
    for iFile = 1:nBehaviorNeeded
        assert(exist(behaviorNwbFiles{iFile}, 'file') == 2, ...
            'Figure6:MissingBehaviorFile', ...
            'Behavior NWB file not found: %s', behaviorNwbFiles{iFile});
    end
end

results = struct();
results.generatedOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');
results.panelsRequested = selected;
results.deviations = documentedDeviations();
results.sampleSizes = struct();
results.skipped = {};

% ---------------------------------------------------------------- panels B, C
if needOpto
    fprintf('loading sparseOpto session ...\n');
    opto = loadFigure6SessionFromNWB(char(optoNwbFile));
    assert(strcmp(opto.sessionType, 'sparseOpto'), 'Figure6:WrongOptoFile', ...
        '%s is not a sparseOpto file.', optoNwbFile);
    results.optoSession = sessionSummary(opto);

    if ismember('B', selected)
        fprintf('panel B: streaming sparseOpto trial averages ...\n');
        results.panelB = plotSubtypeStimulusResponses(opto, outputDirectory);
        results.sampleSizes.panelB_cells = results.panelB.cellCount;
    end
    if ismember('C', selected)
        fprintf('panel C ...\n');
        results.panelC = plotIntrinsicPropertyScatter(opto, outputDirectory);
        results.sampleSizes.panelC_cells = results.panelC.cellCount;
    end
end

% ------------------------------------------------------------------- panel D
if needConn
    fprintf('loading sparsePulseHad session ...\n');
    conn = loadFigure6SessionFromNWB(char(connNwbFile));
    assert(strcmp(conn.sessionType, 'sparsePulseHad'), 'Figure6:WrongConnFile', ...
        '%s is not a sparsePulseHad file.', connNwbFile);
    results.connSession = sessionSummary(conn);

    fprintf('panel D: cell-type average IPSP waveform ...\n');
    results.panelD = plotSubtypeIpspWaveforms(conn, outputDirectory);
    results.sampleSizes.panelD_pairs = results.panelD.pairCount;
end

% ------------------------------------------------- panels E, F, H, I, J, K, G
behaviorGrouping = struct('raster', [], 'confusion', []);
if needDisplayBehavior || needAllBehavior
    fprintf('loading behavior session %s ...\n', behaviorNwbFiles{1});
    beh = loadFigure6SessionFromNWB(char(behaviorNwbFiles{1}));
    assert(strcmp(beh.sessionType, 'behavior'), 'Figure6:WrongBehaviorFile', ...
        '%s is not a behavior file.', behaviorNwbFiles{1});
    results.behaviorSession = sessionSummary(beh);
end

if needDisplayBehavior
    needStream = any(ismember({'H', 'I', 'J', 'K'}, selected));
    if needStream
        fprintf('streaming behavior voltage ...\n');
        stream = accumulateBehaviorStatistics(beh);
        results.behaviorStream = streamSummary(stream);
    end

    if ismember('E', selected)
        fprintf('panel E: UMAP-sorted raster (chunk guard +/-100) ...\n');
        behaviorGrouping.raster = umapGrouping(beh, 100);
        if isempty(behaviorGrouping.raster)
            results.skipped{end+1} = 'E (no UMAP implementation on the path)';
        else
            results.panelE = plotUmapSortedDynamics(beh, behaviorGrouping.raster, ...
                outputDirectory);
        end
    end
    if ismember('F', selected)
        fprintf('panel F: UMAP embedding (chunk guard +/-5000) ...\n');
        behaviorGrouping.confusion = umapGrouping(beh, 5000);
        if isempty(behaviorGrouping.confusion)
            results.skipped{end+1} = 'F (no UMAP implementation on the path)';
        else
            results.panelF = plotUmapEmbedding(behaviorGrouping.confusion, ...
                behaviorGrouping.raster, outputDirectory);
        end
    end
    if ismember('H', selected)
        fprintf('panel H: subtype average dynamics ...\n');
        results.panelH = plotSubtypeSpontaneousDynamics(beh, stream, outputDirectory);
        results.sampleSizes.panelH_cells = results.panelH.cellCount;
    end
    if ismember('I', selected)
        fprintf('panel I: oscillation-triggered average ...\n');
        results.panelI = plotOscillationTriggeredAverage(beh, stream, outputDirectory);
        results.sampleSizes.panelI_events = results.panelI.eventCount;
    end
    if ismember('J', selected)
        fprintf('panel J: whisking-triggered average ...\n');
        results.panelJ = plotWhiskingTriggeredAverage(beh, stream, outputDirectory);
        results.sampleSizes.panelJ_events = results.panelJ.eventCount;
    end
    if ismember('K', selected)
        fprintf('panel K: state-dependent excitability ...\n');
        results.panelK = plotStateDependentExcitability(beh, stream, outputDirectory);
        results.sampleSizes.panelK_cells = results.panelK.cellCount;
    end
end

if needAllBehavior
    fprintf('panel G: pooled NPY vs UMAP-group confusion matrix ...\n');
    results.panelG = plotConfusionMatrix(behaviorNwbFiles, beh, ...
        behaviorGrouping.confusion, outputDirectory);
    if isempty(results.panelG)
        results = rmfield(results, 'panelG');
        results.skipped{end+1} = 'G (no UMAP implementation on the path)';
    else
        results.sampleSizes.panelG_cells = results.panelG.cellCount;
    end
end

resultsPath = fullfile(outputDirectory, 'figure6_nwb_reference_results.mat');
saveResults(resultsPath, results);
fprintf('Figure 6 panel(s) %s written to %s\n', strjoin(selected, ','), outputDirectory);
end

% =========================================================================
function selectedPanels = normalizePanelSelection(panels)
valid = {'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'};
if isempty(panels) || (ischar(panels) && strcmpi(strtrim(panels), 'all'))
    selectedPanels = valid;
    return
end
if ischar(panels)
    tokens = strtrim(strsplit(panels, ','));
elseif iscellstr(panels) %#ok<ISCLSTR> R2019b compatible
    tokens = strtrim(panels);
else
    error('Figure6:PanelSelection', ...
        'PANELS must be a character vector or a cell array of character vectors.');
end
tokens = upper(tokens(~cellfun(@isempty, tokens)));
assert(~ismember('A', tokens), 'Figure6:PanelANotReproduced', ...
    ['Panel A is the spinning-disk confocal image. Its source TIFFs are not ', ...
     'in any NWB file and the author ruled it need not be reproduced, so it ', ...
     'is not a valid selector. Valid: %s.'], strjoin(valid, ', '));
unknown = setdiff(tokens, valid);
assert(isempty(unknown), 'Figure6:UnknownPanel', ...
    'Unknown panel selector(s): %s. Valid: %s.', strjoin(unknown, ', '), ...
    strjoin(valid, ', '));
selectedPanels = valid(ismember(valid, tokens));
end

function deviations = documentedDeviations()
deviations = { ...
    ['Motion-corrected spike times are not in the NWB (mcTrace_all was not ', ...
     'carried into the prepared MAT). The author waived the step on ', ...
     '2026-09-18, and the Figure 6 behavior scripts already have the call ', ...
     'commented out, so stored raw spike times are used.']; ...
    ['Panel B''s source cell selection adds a 2-sigma spike-height term that ', ...
     '/units/passed_qc does not carry. The author waived it on 2026-09-19; it ', ...
     'was measured to exclude no cell in M-YQ0201-29 sparseOpto. See ', ...
     'results.panelB.selectionCheck.']; ...
    ['UMAP and k-means are unseeded in the source and there is no rng call ', ...
     'anywhere in the Figure 6 scripts, so panels E, F and G are not ', ...
     'bit-reproducible. The k-means labels are aligned to the NPY subtype ', ...
     'firing-rate averages after clustering so panel G''s columns cannot ', ...
     'swap between runs. That fixes the orientation of the comparison only.']; ...
    ['The source''s throwaway 1-component UMAP is not computed. In the raster ', ...
     'script its sort order is immediately overwritten by the linear ', ...
     'discriminant; in the confusion-matrix script it orders only a label bar ', ...
     'that is not a Figure 6 panel.']; ...
    ['Third-party colour maps (colorcet) are replaced with built-in colours, ', ...
     'and the panel E raster is drawn inline rather than through the lab ', ...
     'raster_plot helper. Cosmetic only.']};
end

function summary = sessionSummary(session)
summary = struct( ...
    'nwbFile', session.nwbFile, ...
    'sessionType', session.sessionType, ...
    'numberOfCells', session.numberOfCells, ...
    'figureCellCount', numel(session.figureCellIndices), ...
    'npyPositiveInFigureCells', sum(session.gfpPositive(session.figureCellIndices)), ...
    'npyNegativeInFigureCells', sum(~session.gfpPositive(session.figureCellIndices)), ...
    'pixelSizeMicrometers', session.pixelSizeMicrometers);
if isfield(session, 'numberOfVoltageFrames')
    summary.numberOfVoltageFrames = session.numberOfVoltageFrames;
    summary.framesPerChunk = session.framesPerChunk;
    summary.numberOfChunks = session.numberOfChunks;
    summary.framePeriodSeconds = session.framePeriodSeconds;
end
end

function summary = streamSummary(stream)
summary = struct( ...
    'voltageEdges', stream.voltageEdges, ...
    'analysisFrameCount', stream.analysisFrameCount, ...
    'quietFrameCount', stream.quietFrameCount, ...
    'whiskOnsetFrameCount', stream.whiskOnsetFrameCount, ...
    'rejectedHighSynchronyFrames', stream.rejectedHighSynchronyFrames);
end

% =========================================================================
function stream = accumulateBehaviorStatistics(session)
% One streaming pass over the behavior recording that builds every full-length
% quantity panels H, I, J and K need. Blocks are aligned to whole recording
% chunks, which keeps the chunk-boundary guards exact.
%
% Accumulated, all restricted to the passed_qc cells:
%   populationSubthreshold  [1 x nFrames]        source vAllAvg      (panel I)
%   subtypeSubthreshold     [2 x nFrames]        source vAvg         (H, I, J)
%   voltageHistogram        [nCells x nV x 2]    source vHistAll     (panel K)
%   spikeVoltageHistogram   [nCells x nV x 2]    source vSpkHistAll  (panel K)
%
% The per-cell subthreshold traces are never held in memory. Panel I's source
% builds per-cell spike-triggered averages and then averages them by subtype;
% a triggered average is linear in the trace, so averaging by subtype first and
% triggering afterwards gives the identical result at a fraction of the cost.

dt = session.framePeriodSeconds;
nFrames = session.numberOfVoltageFrames;
qc = session.figureCellIndices;
nQC = numel(qc);
isNpy = session.gfpPositive(qc);

% ------------------------------------------------------------ panel K masks
% Source fig6_NPYpm_pSpkV_whisk.m lines 354-381.
tMask = false(1, nFrames);
tMask(1:min(round(2000 / dt), nFrames)) = true;
tMask(chunkBoundaryGuard(nFrames, session.framesPerChunk, 140)) = false;

runMask = runningMask(session);
whisk = selectWhiskingOnsets(session);
onsetKernel = zeros(1, 2 * 2500 + 1);
onsetKernel(2500 + 1 + (0:500)) = 1;
whiskOnsetMask = conv(double(whisk.onsetIndicator), onsetKernel, 'same') > 0;

stateMask = [tMask & session.quietMask; ...
             tMask & whiskOnsetMask & ~runMask];

voltageEdges = linspace(0, 1.5, 51);
nV = numel(voltageEdges) - 1;

stream = struct();
stream.voltageEdges = voltageEdges;
stream.voltageCenters = binCenters(voltageEdges);
stream.stateLabels = {'Quiet', 'Whisk'};
stream.analysisFrameCount = sum(tMask);
stream.quietFrameCount = sum(stateMask(1, :));
stream.whiskOnsetFrameCount = sum(stateMask(2, :));
stream.whisk = whisk;
stream.runMask = runMask;

populationSubthreshold = zeros(1, nFrames);
subtypeSubthreshold = zeros(2, nFrames);
voltageHistogram = zeros(nQC, nV, 2);
spikeVoltageHistogram = zeros(nQC, nV, 2);

spikeMaskByCell = spikeIndicatorRows(session, qc);

blockFrames = session.framesPerChunk;
nBlocks = ceil(nFrames / blockFrames);
for iBlock = 1:nBlocks
    firstFrame = (iBlock - 1) * blockFrames + 1;
    lastFrame = min(iBlock * blockFrames, nFrames);
    block = double(session.readVoltageRows('subthreshold_voltage', ...
        qc(1):qc(end), [firstFrame lastFrame]));
    block = block(qc - qc(1) + 1, :);
    columns = firstFrame:lastFrame;

    populationSubthreshold(columns) = mean(block, 1);
    subtypeSubthreshold(1, columns) = mean(block(~isNpy, :), 1);
    subtypeSubthreshold(2, columns) = mean(block(isNpy, :), 1);

    for iState = 1:2
        inState = stateMask(iState, columns);
        if ~any(inState)
            continue
        end
        for iCell = 1:nQC
            trace = block(iCell, inState);
            voltageHistogram(iCell, :, iState) = ...
                voltageHistogram(iCell, :, iState) + histcounts(trace, voltageEdges);
            spikes = spikeMaskByCell(iCell, columns) & inState;
            if any(spikes)
                spikeVoltageHistogram(iCell, :, iState) = ...
                    spikeVoltageHistogram(iCell, :, iState) + ...
                    histcounts(block(iCell, spikes), voltageEdges);
            end
        end
    end
    fprintf('  block %d/%d\n', iBlock, nBlocks);
end

stream.populationSubthreshold = populationSubthreshold;
stream.subtypeSubthreshold = subtypeSubthreshold;
stream.voltageHistogram = voltageHistogram;
stream.spikeVoltageHistogram = spikeVoltageHistogram;
stream.isNpy = isNpy;

% ---------------------------------------------------- spike-derived products
[rate, rejected] = subtypeSpikeRates(session);
stream.subtypeSpikeRate = rate.subtypeMean;
stream.subtypeSpikeRateRejected = rate.subtypeMeanRejected;
stream.rejectedHighSynchronyFrames = rejected;
end

% =========================================================================
function mask = spikeIndicatorRows(session, cellIndices)
% Logical [numel(cellIndices) x nFrames] spike indicator, the source
% spk_t_mat restricted to the figure cells.
nFrames = session.numberOfVoltageFrames;
mask = false(numel(cellIndices), nFrames);
for k = 1:numel(cellIndices)
    frames = session.spikeFrameIndex{cellIndices(k)};
    mask(k, frames) = true;
end
end

function [rate, rejectedFrames] = subtypeSpikeRates(session)
% Source spk_t_mat and spk_t_mat_c. The rejection zeroes whole frames in which
% more than 15 of the figure cells spike simultaneously, and it tests
% sum(spk_t_mat) while assigning into spk_t_mat_c. Both the raw and the
% rejected subtype means are returned: panel I uses the raw one, panels H and J
% the rejected one.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
isNpy = session.gfpPositive(qc);
mask = spikeIndicatorRows(session, qc);

perFrame = sum(mask, 1);
reject = perFrame > 15;
rejectedFrames = sum(reject);

rate = struct();
rate.subtypeMean = [mean(mask(~isNpy, :), 1); mean(mask(isNpy, :), 1)] / dt;
maskRejected = mask;
maskRejected(:, reject) = false;
rate.subtypeMeanRejected = ...
    [mean(maskRejected(~isNpy, :), 1); mean(maskRejected(isNpy, :), 1)] / dt;
end

% =========================================================================
function runMask = runningMask(session)
% Source tRun = conv(abs(loc) > 1e-3, ker, 'same') > 0 with a flat +/-1000
% frame kernel.
kernel = ones(1, 2001);
runMask = conv(double(abs(session.runningSpeed) > 1e-3), kernel, 'same') > 0;
end

function interpolated = interpolateWhisking(session, smoothFirst)
% Source lines 77-84. The face-camera trace is lifted onto the voltage timebase
% chunk by chunk with a bilinear imresize, optionally after an 8-sample moving
% median. mWhiskInterp is the raw version, mWhiskSmInterp the smoothed one.
trace = session.whiskingMotion;
if smoothFirst
    trace = movmedian(trace, 8);
end
byChunk = reshape(trace, session.faceFramesPerChunk, session.numberOfChunks);
resized = imresize(byChunk, [session.framesPerChunk session.numberOfChunks], ...
    'bilinear');
interpolated = reshape(resized, 1, []);
assert(numel(interpolated) == session.numberOfVoltageFrames, ...
    'Figure6:WhiskInterpolation', ...
    'Interpolated whisking trace does not match the voltage frame count.');
end

function quietMask = medianQuietMask(session)
% Source fig6_NPYpm_avgDynamics_raster.m line 235, the median rule. This is
% NOT the stored quiet_mask, which follows the Gaussian-mixture rule.
smoothed = interpolateWhisking(session, true);
quietMask = ~(movmean(smoothed > median(smoothed), 100) > 0);
end

function whiskMask = medianWhiskMask(session, runMask)
% Source line 236.
smoothed = interpolateWhisking(session, true);
whiskMask = smoothed > median(smoothed) & ~runMask;
end

% =========================================================================
function whisk = selectWhiskingOnsets(session)
% Source fig6_NPY_osc_whisk_TrigAvg_ranksumStats.m lines 343-392, shared by
% panels J and K.
%
% The regularisation strength is chosen as the one whose onset-triggered
% whisking average has the smallest pre-onset variance, then onsets are
% thresholded at 1e-3 and thinned so that consecutive onsets are at least 1.5 s
% apart. A different lambda is chosen per mouse (14 of 15 for M-YQ0201-27,
% 9 of 15 for M-YQ0201-29), so it must not be hard-coded.
dt = session.framePeriodSeconds;
nFrames = session.numberOfVoltageFrames;
reconstruction = session.whiskingMotionReconstruction;
nLambda = size(reconstruction, 2);

faceDt = 25e-3;
nBack = round(2 / faceDt);
nFront = round(3 / faceDt);
preOnsetStd = nan(nLambda, 1);
for iLambda = 1:nLambda
    active = reconstruction(:, iLambda)' > 1e-3;
    onsets = find([0 diff(active)] > 0);
    average = triggeredAverage(onsets, session.whiskingMotion, nBack, nFront);
    value = std(average(1:nBack));
    if value ~= 0
        preOnsetStd(iLambda) = value;
    end
end
[~, lambdaIndex] = min(preOnsetStd);

% The selected reconstruction is lifted onto the voltage timebase the same way
% the source does, with a single bilinear imresize of the whole 15-column array.
lifted = imresize(reconstruction', [nLambda nFrames], 'bilinear');
active = lifted(lambdaIndex, :) > 1e-3;
onsets = find([0 diff(active)] > 0);
nBeforeThinning = numel(onsets);
onsets = onsets([inf diff(onsets)] > 1.5 / dt);

whisk = struct();
whisk.lambdaIndex = lambdaIndex;
whisk.preOnsetStd = preOnsetStd;
whisk.onsetFrames = onsets;
whisk.onsetCountBeforeThinning = nBeforeThinning;
whisk.onsetCount = numel(onsets);
whisk.onsetIndicator = false(1, nFrames);
whisk.onsetIndicator(onsets) = true;
end

% =========================================================================
function average = triggeredAverage(triggerFrames, traces, nBack, nFront)
% Inline equivalent of get_sta_mat_single for one trigger set.
%
% traces  : [nRows x nFrames]
% average : [nRows x nBack+nFront+1]
%
% Triggers closer than nBack to the start or nFront to the end are dropped, and
% the average across triggers is a nanmean, so NaN-masked samples drop out
% exactly as they do in the source.
[nRows, nFrames] = size(traces);
average = zeros(nRows, nBack + nFront + 1);
triggerFrames = triggerFrames(:);
triggerFrames = triggerFrames(triggerFrames >= nBack + 1 & ...
    triggerFrames <= nFrames - nFront);
if isempty(triggerFrames)
    return
end
windows = triggerFrames + (-nBack:nFront);
for iRow = 1:nRows
    row = traces(iRow, :);
    average(iRow, :) = nanmean(row(windows), 1);
end
end

function guard = chunkBoundaryGuard(nFrames, framesPerChunk, halfWidth)
% Source pattern (1:nFrame:nFramesTotal+1) + (-halfWidth:halfWidth)', clipped
% into range. Returns linear frame indices to exclude.
boundaries = (1:framesPerChunk:nFrames + 1)';
guard = boundaries + (-halfWidth:halfWidth);
guard = max(min(guard(:), nFrames), 1);
end

function centers = binCenters(edges)
centers = mean([edges(1:end-1); edges(2:end)], 1);
end

% =========================================================================
% Panel B
% =========================================================================
function panel = plotSubtypeStimulusResponses(session, outputDirectory)
% Figure 6B. Source fig6_NPYpm_optoAvgTraces.m lines 227-287.
%
% The voltage trace is the source traces_all_sub_raw: the recording with a
% [-3:5] frame window around every spike replaced by linear interpolation, and
% no smoothing. normalized_voltage is already divided by the spike height, and
% linear interpolation commutes with that scaling, so interpolating the stored
% normalized trace reproduces traces_all_sub_raw ./ spkHgtNoBlue exactly.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
isNpy = session.gfpPositive(qc);

panel = struct();
panel.panel = 'B';
panel.selectionCheck = twoSigmaSelectionCheck(session);

trial = accumulateStimulusTrialAverages(session);
% The accumulator works on the full 320-cell units axis, so restrict to the
% figure cells before applying the NPY mask, which is on the figure-cell axis.
voltageAverage = trial.voltageAverage(qc, :);
rateAverage = trial.rateAverage(qc, :);
voltageByType = [mean(voltageAverage(~isNpy, :), 1); ...
                 mean(voltageAverage(isNpy, :), 1)];
rateByType = [mean(rateAverage(~isNpy, :), 1); ...
              mean(rateAverage(isNpy, :), 1)] / dt;

nBinStep = 4;
nBinRamp = 20;
nFrame = session.framesPerChunk;
timeSeconds = (1:nFrame) * dt;

panel.cellCount = numel(qc);
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.stimulatedTrialsPerCell = trial.trialCount(qc);
panel.timeSeconds = timeSeconds;
panel.blueCommand = trial.blueCommand;
panel.subtypeVoltage = voltageByType;
panel.subtypeRateHz = rateByType;
panel.stepWindowSeconds = [2.5 3] + [-0.05 0.05];
panel.rampWindowSeconds = [4.5 9.5] + [-0.25 0.25];
panel.stepRateTimeSeconds = imresize(timeSeconds, 1 / nBinStep, 'box');
panel.stepRateHz = imresize(rateByType, [2 nFrame / nBinStep], 'box');
panel.rampRateTimeSeconds = imresize(timeSeconds, 1 / nBinRamp, 'box');
panel.rampRateHz = imresize(rateByType, [2 nFrame / nBinRamp], 'box');

colours = npyColours();
figureHandle = figure('Name', 'Figure 6B', 'Color', 'w', ...
    'Position', [80 80 900 620]);

ax(1) = subplot(3, 2, 1);
plot(ax(1), timeSeconds, trial.blueCommand, 'Color', [0 0 0.85], 'LineWidth', 1.2);
ylabel(ax(1), 'Blue (a.u.)');
title(ax(1), 'Step stimulation');
xlim(ax(1), panel.stepWindowSeconds);

ax(2) = subplot(3, 2, 2);
plot(ax(2), timeSeconds, trial.blueCommand, 'Color', [0 0 0.85], 'LineWidth', 1.2);
title(ax(2), 'Ramp stimulation');
xlim(ax(2), panel.rampWindowSeconds);

ax(3) = subplot(3, 2, 3);
plotSubtypes(ax(3), timeSeconds, voltageByType, colours);
ylabel(ax(3), 'Voltage (norm)');
xlim(ax(3), panel.stepWindowSeconds);

ax(4) = subplot(3, 2, 4);
plotSubtypes(ax(4), timeSeconds, voltageByType, colours);
xlim(ax(4), panel.rampWindowSeconds);

ax(5) = subplot(3, 2, 5);
plotSubtypes(ax(5), panel.stepRateTimeSeconds, panel.stepRateHz, colours);
xlabel(ax(5), 'Time (s)');
ylabel(ax(5), 'Firing rate (Hz)');
xlim(ax(5), panel.stepWindowSeconds);

ax(6) = subplot(3, 2, 6);
plotSubtypes(ax(6), panel.rampRateTimeSeconds, panel.rampRateHz, colours);
xlabel(ax(6), 'Time (s)');
xlim(ax(6), panel.rampWindowSeconds);
legend(ax(6), session.npyLabels, 'Location', 'northwest', 'Box', 'off');

linkaxes(ax([1 3 5]), 'x');
linkaxes(ax([2 4 6]), 'x');
saveFigure(figureHandle, outputDirectory, 'figure6B_subtype_stimulus_responses_raw.png');
end

function trial = accumulateStimulusTrialAverages(session)
% One streaming pass over the sparseOpto recording. For every cell, average the
% spike-interpolated voltage and the spike train over the recording chunks in
% which that cell received blue light.
%
% Blocks carry a halo so the spike interpolation matches the source's single
% whole-recording interp1 across block boundaries. The chunk-boundary spike
% guard is the source's +/-140 frame t_noise.
nFrames = session.numberOfVoltageFrames;
nFrame = session.framesPerChunk;
nCells = session.numberOfCells;
qc = session.figureCellIndices;

spikeGuard = false(1, nFrames);
spikeGuard(chunkBoundaryGuard(nFrames, nFrame, 140)) = true;

voltageSum = zeros(nCells, nFrame);
rateSum = zeros(nCells, nFrame);
trialCount = zeros(nCells, 1);
blueCommand = [];

chunksPerBlock = max(1, floor(10000 * 8 / nFrame));
haloFrames = 2000;
spikeWindow = -3:5;

nBlocks = ceil(session.numberOfChunks / chunksPerBlock);
for iBlock = 1:nBlocks
    firstChunk = (iBlock - 1) * chunksPerBlock + 1;
    lastChunk = min(iBlock * chunksPerBlock, session.numberOfChunks);
    firstFrame = (firstChunk - 1) * nFrame + 1;
    lastFrame = lastChunk * nFrame;
    haloFirst = max(1, firstFrame - haloFrames);
    haloLast = min(nFrames, lastFrame + haloFrames);

    voltage = double(session.readVoltageRows('normalized_voltage', ...
        1:nCells, [haloFirst haloLast]));
    blue = double(session.readBlueRows(1:nCells, [firstFrame lastFrame]));
    nHalo = size(voltage, 2);

    for k = 1:numel(qc)
        iCell = qc(k);
        blueRow = reshape(blue(iCell, :), nFrame, []);
        stimulated = sum(blueRow, 1) > 0;
        if ~any(stimulated)
            continue
        end
        if isempty(blueCommand) && iCell == qc(1)
            firstStim = find(stimulated, 1);
            blueCommand = blueRow(:, firstStim)';
        end

        interpolated = interpolateOverSpikes(voltage(iCell, :), ...
            session.spikeFrameIndex{iCell} - haloFirst + 1, spikeWindow, nHalo);
        trimmed = interpolated(firstFrame - haloFirst + 1 : lastFrame - haloFirst + 1);
        voltageStack = reshape(trimmed, nFrame, []);

        spikes = false(1, lastFrame - firstFrame + 1);
        frames = session.spikeFrameIndex{iCell};
        frames = frames(frames >= firstFrame & frames <= lastFrame);
        spikes(frames - firstFrame + 1) = true;
        spikes(spikeGuard(firstFrame:lastFrame)) = false;
        spikeStack = reshape(spikes, nFrame, []);

        voltageSum(iCell, :) = voltageSum(iCell, :) + ...
            sum(voltageStack(:, stimulated), 2)';
        rateSum(iCell, :) = rateSum(iCell, :) + ...
            sum(spikeStack(:, stimulated), 2)';
        trialCount(iCell) = trialCount(iCell) + sum(stimulated);
    end
    fprintf('  chunks %d-%d of %d\n', firstChunk, lastChunk, session.numberOfChunks);
end

assert(~isempty(blueCommand), 'Figure6:NoStimulatedTrial', ...
    'No stimulated recording chunk was found for the first figure cell.');
safeCount = max(trialCount, 1);
trial = struct();
trial.voltageAverage = voltageSum ./ safeCount;
trial.rateAverage = rateSum ./ safeCount;
trial.trialCount = trialCount;
trial.blueCommand = blueCommand;
end

function interpolated = interpolateOverSpikes(trace, spikeFrames, spikeWindow, nFrames)
% Source traces_all_sub_raw: drop a window around every spike and fill it by
% linear interpolation from the surviving samples.
spikeFrames = spikeFrames(:)';
if isempty(spikeFrames)
    interpolated = trace;
    return
end
drop = spikeFrames + spikeWindow';
drop = drop(drop >= 1 & drop <= nFrames);
keep = true(1, nFrames);
keep(drop) = false;
if all(keep)
    interpolated = trace;
    return
end
assert(sum(keep) >= 2, 'Figure6:SpikeInterpolation', ...
    'Too few non-spike samples remain to interpolate the subthreshold trace.');
interpolated = interp1(find(keep), trace(keep), 1:nFrames, 'linear');
interpolated(isnan(interpolated)) = trace(isnan(interpolated));
end

function check = twoSigmaSelectionCheck(session)
% Panel B's source adds spkHgtNoBlue > mean - 2*std to the selection that
% /units/passed_qc carries. spkHgtNoBlue is not deposited, so the term cannot be
% re-applied. It was measured to be a strict no-op in M-YQ0201-29 sparseOpto and
% to drop 2 of 286 cells in M-YQ0201-27, so this records the waiver and warns if
% the file being read is not the one the waiver was measured on.
check = struct();
check.term = 'spkHgtNoBlue > mean(spkHgtNoBlue) - 2*std(spkHgtNoBlue)';
check.applied = false;
check.reason = ['Waived by the author on 2026-09-19. spkHgtNoBlue is not a ', ...
    '/units column, so the term cannot be re-applied from NWB. It was ', ...
    'measured to exclude no cell in M-YQ0201-29 sparseOpto, where panel B is ', ...
    'drawn, and to drop 2 of 286 cells in M-YQ0201-27 sparseOpto.'];
check.measuredNoOpForCellCount = 257;
check.cellCount = numel(session.figureCellIndices);
if check.cellCount ~= check.measuredNoOpForCellCount
    warning('Figure6:TwoSigmaWaiverUnverified', ...
        ['Panel B is being drawn on a session with %d figure cells, not the ', ...
         '%d of M-YQ0201-29 sparseOpto where the 2-sigma spike-height term ', ...
         'was measured to be a no-op. The waiver has not been verified here; ', ...
         'see FIGURE6_PHASE1_HANDOFF.md gap 2.'], check.cellCount, ...
        check.measuredNoOpForCellCount);
end
end

% =========================================================================
% Panel C
% =========================================================================
function panel = plotIntrinsicPropertyScatter(session, outputDirectory)
% Figure 6C. Source fig6_NPYpm_ephys.m lines 139-159, figure(20).
% props_all(:,1) is the after-depolarization and (:,2) the optical rheobase,
% matching the published V_ADP and V_Rh axes.
qc = session.figureCellIndices;
adp = session.intrinsicProperties(qc, 1);
rheobase = session.intrinsicProperties(qc, 2);
isNpy = session.gfpPositive(qc);

xLimits = [-0.05 0.4];
yLimits = [-0.2 0.7];

panel = struct();
panel.panel = 'C';
panel.cellCount = numel(qc);
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.afterDepolarization = adp;
panel.opticalRheobase = rheobase;
panel.npyPositive = isNpy;
panel.axisLimits = struct('x', xLimits, 'y', yLimits);
panel.insideAxisLimits = sum(adp >= xLimits(1) & adp <= xLimits(2) & ...
    rheobase >= yLimits(1) & rheobase <= yLimits(2));
panel.medianByType = [median(adp(~isNpy)) median(rheobase(~isNpy)); ...
                      median(adp(isNpy))  median(rheobase(isNpy))];

colours = npyColours();
figureHandle = figure('Name', 'Figure 6C', 'Color', 'w', ...
    'Position', [100 100 520 470]);
ax = axes(figureHandle);
hold(ax, 'on');
scatter(ax, adp(~isNpy), rheobase(~isNpy), 36, colours(1, :), 'filled');
scatter(ax, adp(isNpy), rheobase(isNpy), 36, colours(2, :), 'filled');
hold(ax, 'off');
xlim(ax, xLimits);
ylim(ax, yLimits);
xlabel(ax, 'V_{ADP} (norm)');
ylabel(ax, 'V_{Rh} (norm)');
title(ax, sprintf('n = %d cells (%d NPY-, %d NPY+)', panel.cellCount, ...
    panel.npyNegativeCount, panel.npyPositiveCount));
legend(ax, session.npyLabels, 'Location', 'northeast', 'Box', 'off');
box(ax, 'on');
set(ax, 'FontSize', 11);

saveFigure(figureHandle, outputDirectory, 'figure6C_intrinsic_property_scatter_raw.png');
end

% =========================================================================
% Panel D
% =========================================================================
function panel = plotSubtypeIpspWaveforms(session, outputDirectory)
% Figure 6D. Source fig6_NPYpm_conn.m lines 160-271, figure(1).
%
% Since 2026-09-24 the source takes pMat from the unwhitened IPSP detection
% (ipspThetaFitsUnwhitened, lines 57-116), the same mask dataPrepsForMAT.m
% stored as synaptic_connection, so the stored column is used as is. The older
% permutation pMat from synapticConn_expFit.mat (line 35) is archival only.
%
% Three stages, all of which matter:
%   1. fit a pw_doubleExp template to the average IPSP waveform of the nearest
%      unconnected pairs, 60-120 um;
%   2. regress that template out of every unconnected pair in the 60-400 um
%      band, over the first 18 post-stimulus frames only, keeping the slope
%      coefficient and discarding the fitted intercept;
%   3. average the corrected waveforms by presynaptic and postsynaptic NPY
%      status over all pairs in the 60-400 um band, connected and unconnected
%      alike, with a plain nanmean.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
nQC = numel(qc);
isNpy = session.gfpPositive(qc);
dx = session.pixelSizeMicrometers;

distance = session.cellDistancePixels(qc, qc);
connected = session.metrics.synapticConnection(qc, qc);

fprintf('  reading the pairwise IPSP waveform ...\n');
waveform = session.readPairwiseSeries('blue_pulse_ipsp', ...
    1:session.numberOfCells, 1:session.numberOfCells, []);
waveform = double(waveform(qc, qc, :));
waveform(isinf(waveform)) = NaN;
nTime = size(waveform, 3);
zeroLagIndex = (nTime + 1) / 2;
assert(zeroLagIndex == round(zeroLagIndex), 'Figure6:IpspWaveformLength', ...
    'The pairwise IPSP waveform must have an odd sample count.');
nBack = zeroLagIndex - 1;
lagFrames = (-nBack:nBack)';

% ------------------------------------------------ stage 1, template fit
nearUnconnected = distance > 60 / dx & distance < 120 / dx & ~connected;
controlAverage = squeeze(nanmean(waveform .* nanMask(nearUnconnected), [1 2]));
startAmplitude = -max(controlAverage(zeroLagIndex:end));
templateFit = fitDoubleExponential(lagFrames, controlAverage, startAmplitude);
template = evaluateDoubleExponential(lagFrames, templateFit.coefficients);

% -------------------------------------- stage 2, remove the blue-light artifact
distanceMask = distance > 60 / dx & distance < 400 / dx;
[preIndex, postIndex] = find(~connected & distanceMask);
fitIndex = zeroLagIndex + (0:17);
corrected = waveform;
coefficients = zeros(numel(preIndex), 1);
for iPair = 1:numel(preIndex)
    trace = squeeze(waveform(preIndex(iPair), postIndex(iPair), :));
    coefficient = regressionSlope(trace(fitIndex)', template(fitIndex)');
    coefficients(iPair) = coefficient;
    corrected(preIndex(iPair), postIndex(iPair), :) = trace - template * coefficient;
end

% --------------------------------------------- stage 3, cell-type averages
subtypeAverage = zeros(2, 2, nTime);
pairCount = zeros(2, 2);
for iPre = 1:2
    for iPost = 1:2
        mask = distanceMask & (isNpy == (iPre - 1)) & (isNpy == (iPost - 1))';
        pairCount(iPre, iPost) = sum(mask(:));
        subtypeAverage(iPre, iPost, :) = ...
            squeeze(nanmean(corrected .* nanMask(mask), [1 2]));
    end
end

blueWaveform = bluePulseCommand(session, nBack, nBack);

panel = struct();
panel.panel = 'D';
panel.cellCount = nQC;
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.lagSeconds = lagFrames * dt;
panel.subtypeAverage = subtypeAverage;
panel.pairCountMatrix = pairCount;
panel.pairCount = sum(pairCount(:));
panel.pairCountLabels = {'NPY- -> NPY-', 'NPY- -> NPY+', 'NPY+ -> NPY-', 'NPY+ -> NPY+'};
panel.captionPairCounts = [pairCount(1, 1) pairCount(1, 2) pairCount(2, 2)];
panel.templateCoefficients = templateFit.coefficients;
panel.templateCoefficientNames = templateFit.names;
panel.templateGoodnessOfFit = templateFit.gof;
panel.controlAverage = controlAverage;
panel.template = template;
panel.correctedPairCount = numel(preIndex);
panel.correctionCoefficientSummary = [min(coefficients) median(coefficients) ...
    max(coefficients)];
panel.blueCommand = blueWaveform;
panel.distanceBandMicrometers = [60 400];

% The manuscript arranges these as a 2 x 2 grid, presynaptic NPY status down
% the rows and postsynaptic across the columns, with the stimulation pulse
% shaded on each trace. The source script overlays all four on one axis; the
% grid is easier to compare against the published panel and changes nothing.
colours = npyColours();
pulseOn = panel.lagSeconds(find(blueWaveform > 0.5 * max(blueWaveform), 1, 'first'));
pulseOff = panel.lagSeconds(find(blueWaveform > 0.5 * max(blueWaveform), 1, 'last'));
yLimits = [min(subtypeAverage(:)) max(subtypeAverage(:))] + ...
    [-0.1 0.1] * range(subtypeAverage(:));

figureHandle = figure('Name', 'Figure 6D', 'Color', 'w', ...
    'Position', [120 120 760 620]);
ax = gobjects(2, 2);
for iPre = 1:2
    for iPost = 1:2
        ax(iPre, iPost) = subplot(2, 2, (iPre - 1) * 2 + iPost);
        hold(ax(iPre, iPost), 'on');
        patch(ax(iPre, iPost), [pulseOn pulseOff pulseOff pulseOn], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.6 0.6 0.95], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5);
        plot(ax(iPre, iPost), panel.lagSeconds, ...
            squeeze(subtypeAverage(iPre, iPost, :)), ...
            'Color', colours(iPre, :), 'LineWidth', 1.5);
        hold(ax(iPre, iPost), 'off');
        box(ax(iPre, iPost), 'on');
        xlim(ax(iPre, iPost), [-0.05 0.2]);
        ylim(ax(iPre, iPost), yLimits);
        title(ax(iPre, iPost), sprintf('%s -> %s  (n = %d)', ...
            session.npyLabels{iPre}, session.npyLabels{iPost}, ...
            pairCount(iPre, iPost)));
        set(ax(iPre, iPost), 'FontSize', 10);
        if iPre == 2
            xlabel(ax(iPre, iPost), 'Time (s)');
        end
        if iPost == 1
            ylabel(ax(iPre, iPost), 'Voltage (norm)');
        end
    end
end
linkaxes(ax(:), 'xy');

saveFigure(figureHandle, outputDirectory, 'figure6D_subtype_ipsp_waveforms_raw.png');
end

function mask = nanMask(logicalMask)
mask = nan(size(logicalMask));
mask(logicalMask) = 1;
end

function fitResult = fitDoubleExponential(lagFrames, trace, startAmplitude)
% Source fig6_NPYpm_conn.m lines 177-185. pw_doubleExp is a lab helper that
% lives outside the main MATLAB tree, so it is reproduced inline below and the
% fittype is built from a function handle rather than from the name. The two
% forms were verified to give identical coefficient ordering and identical fits.
names = {'a', 'delta1', 'delta2', 'tau1', 'tau2'};
model = fittype(@(a, delta1, delta2, tau1, tau2, x) ...
    doubleExponentialPulse(x, a, delta1, delta2, tau1, tau2), ...
    'independent', 'x', 'coefficients', names);
[fitObject, gof] = fit(double(lagFrames), double(trace), model, ...
    'startPoint', [startAmplitude 0 20 5 50], ...
    'lower', [-inf 0 16 1 1], ...
    'upper', [0 16 100 500 500]);
fitResult = struct();
fitResult.names = names;
fitResult.coefficients = coeffvalues(fitObject);
fitResult.gof = gof;
end

function value = evaluateDoubleExponential(lagFrames, coefficients)
value = doubleExponentialPulse(double(lagFrames), coefficients(1), ...
    coefficients(2), coefficients(3), coefficients(4), coefficients(5));
end

function y = doubleExponentialPulse(t, A, Del1, Del2, tau1, tau2)
% Inline copy of the lab helper pw_doubleExp: a rising exponential from Del1,
% then a decaying exponential from Del2.
y = zeros(size(t));
rising = t > Del1 & t < Del2;
decaying = t >= Del2;
y(rising) = -A * (1 - exp(-(t(rising) - Del1) / tau1));
y(decaying) = -A * (1 - exp(-(Del2 - Del1) / tau1)) .* ...
    exp(-(t(decaying) - Del2) / tau2);
end

function slope = regressionSlope(trace, reference)
% Inline equivalent of SeeResiduals_vec(trace, reference, 1, 1), returning the
% second coefficient: the reference is mean-subtracted and fitted alongside a
% constant term, and the caller keeps only the slope.
reference = reference(:)';
trace = trace(:)';
design = [ones(1, numel(reference)); reference - mean(reference)];
coefficients = trace * design' / (design * design');
slope = coefficients(2);
end

function command = bluePulseCommand(session, nBack, nFront)
% Source fig6_NPYpm_conn.m lines 250-253: the blue command of cell 1, triggered
% on the rising edges of its highest level.
nFrames = min(session.numberOfVoltageFrames, session.framesPerChunk * 4);
blue = double(session.readBlueRows(1:1, [1 nFrames]));
levels = unique(blue);
onsets = find(diff(blue == levels(end)) > 0);
if isempty(onsets)
    command = zeros(1, nBack + nFront + 1);
    return
end
command = triggeredAverage(onsets, blue, nBack, nFront);
end

% =========================================================================
% Panels E, F, G: the UMAP grouping
% =========================================================================
function grouping = umapGrouping(session, guardHalfWidth)
% Source fig6_NPYpm_avgDynamics_raster.m lines 120-160 and
% fig6_NPY_umap_confusionmat.m lines 74-120. The two differ only in the
% chunk-boundary guard: +/-100 frames in the raster script, +/-5000 in the
% confusion-matrix script. Everything else is identical.
%
% Returns [] when no MATLAB UMAP implementation is on the path, so the rest of
% the figure still builds, as reproduce_figure2_from_nwb.m does for its UMAP
% scatter.
if isempty(which('UMAP'))
    warning('Figure6:NoUmap', ...
        ['No MATLAB UMAP implementation found on the path, so panels E, F ', ...
         'and G are skipped. Every other panel is unaffected.']);
    grouping = [];
    return
end

dt = session.framePeriodSeconds;
nFrames = session.numberOfVoltageFrames;
qc = session.figureCellIndices;
nQC = numel(qc);
isNpy = session.gfpPositive(qc);

binFrames = round(500e-3 / dt);
nBinTotal = ceil(nFrames / binFrames);

% Source fr_bin_mat = imresize(spk_t_mat_c/dt, [nCells nBinTotal], 'box').
% imresize with a unit scale on the cell axis is exactly row-separable, which
% was verified numerically, so the rows are built one at a time and the full
% [nCells x nFrames] double matrix is never allocated.
mask = spikeIndicatorRows(session, qc);
reject = sum(mask, 1) > 15;
binnedRate = zeros(nQC, nBinTotal);
for iCell = 1:nQC
    row = double(mask(iCell, :));
    row(reject) = 0;
    binnedRate(iCell, :) = imresize(row / dt, [1 nBinTotal], 'box');
end

tMask = false(1, nFrames);
tMask(1:min(round(2000 / dt), nFrames)) = true;
tMask(chunkBoundaryGuard(nFrames, session.framesPerChunk, guardHalfWidth)) = false;
binMask = imresize(double(tMask), [1 nBinTotal], 'box') > 0.8;

[u, s] = svds(zscore(binnedRate(:, binMask), [], 2), nQC);
embeddingInput = u * s;

umap = UMAP('n_neighbors', 20, 'n_components', 2);
coordinates = umap.fit_transform(embeddingInput);
assert(size(coordinates, 1) == nQC, 'Figure6:UmapOutput', ...
    'UMAP returned %d rows for %d cells.', size(coordinates, 1), nQC);

groupIndex = kmeans(coordinates, 2);
groupIndex = 3 - groupIndex;   % source relabel, kept for provenance

alignment = alignGroupsToCellTypes(groupIndex, binnedRate, isNpy);
groupIndex = alignment.groupIndex;

sortOrder = lineardiscriminantSort(coordinates, groupIndex);

grouping = struct();
grouping.sessionFile = session.nwbFile;
grouping.guardHalfWidth = guardHalfWidth;
grouping.binFrames = binFrames;
grouping.numberOfBins = nBinTotal;
grouping.analysisBinCount = sum(binMask);
grouping.coordinates = coordinates;
grouping.groupIndex = groupIndex;
grouping.sortOrder = sortOrder;
grouping.npyPositive = isNpy;
grouping.groupSizes = [sum(groupIndex == 1) sum(groupIndex == 2)];
grouping.alignment = alignment;
end

function alignment = alignGroupsToCellTypes(groupIndex, binnedRate, isNpy)
% The k-means labels are arbitrary and the source's fixed idxUmap = 3 - idxUmap
% does not stabilise them, so a rerun can swap the two columns of panel G and
% invert panel E's group colour bar.
%
% Compare each cluster's average binned firing-rate trajectory with each cell
% type's and keep whichever of the two labellings puts group 1 with NPY(-) and
% group 2 with NPY(+), the published orientation (157/203 Group 1 cells NPY(-),
% 202/302 Group 2 cells NPY(+)).
%
% This fixes the ORIENTATION of the comparison only. Both labellings carry the
% same off-diagonal mass, which is what panel G reports, so nothing is being
% assumed into the result. A near-tie in the two trace sums would mean the
% clusters do not separate along the cell-type axis at all and the panel
% deserves a second look, so the margin is returned.
groupMean = [mean(binnedRate(groupIndex == 1, :), 1); ...
             mean(binnedRate(groupIndex == 2, :), 1)];
typeMean = [mean(binnedRate(~isNpy, :), 1); ...
            mean(binnedRate(isNpy, :), 1)];
correlation = corr(groupMean', typeMean');
aligned = correlation(1, 1) + correlation(2, 2);
swapped = correlation(1, 2) + correlation(2, 1);

alignment = struct();
alignment.correlation = correlation;
alignment.alignedTrace = aligned;
alignment.swappedTrace = swapped;
alignment.margin = aligned - swapped;
alignment.swapApplied = swapped > aligned;
if alignment.swapApplied
    groupIndex = 3 - groupIndex;
end
alignment.groupIndex = groupIndex;
if abs(alignment.margin) < 0.05
    warning('Figure6:WeakGroupAlignment', ...
        ['The UMAP clusters barely separate along the NPY axis (margin ', ...
         '%.3g). Panel G''s column orientation is then close to arbitrary; ', ...
         'inspect the embedding before trusting it.'], alignment.margin);
end
end

function sortOrder = lineardiscriminantSort(coordinates, groupIndex)
% Source fig6_NPYpm_avgDynamics_raster.m lines 147-160: a Fisher linear
% discriminant on the 2-D UMAP coordinates, used to order the raster. The
% 1-component UMAP the source computes first is immediately overwritten by this
% and is therefore not reproduced.
scatterWithin = zeros(2);
groupMean = zeros(2);
for iGroup = 1:2
    points = coordinates(groupIndex == iGroup, :);
    centre = mean(points, 1);
    scatterWithin = scatterWithin + (points - centre)' * (points - centre);
    groupMean(iGroup, :) = centre;
end
direction = scatterWithin \ diff(groupMean, 1, 1)';
[~, sortOrder] = sort(coordinates * direction, 'ascend');
end

% =========================================================================
% Panel E
% =========================================================================
function panel = plotUmapSortedDynamics(session, grouping, outputDirectory)
% Figure 6E. Source fig6_NPYpm_avgDynamics_raster.m lines 193-230, figures 1
% and 3. The heatmap and raster are drawn only over the displayed window; the
% source computes them over a longer stretch and then sets xlim, which is a
% display choice with no scientific content.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
order = grouping.sortOrder;

displaySeconds = [290 360];
firstFrame = max(1, floor(displaySeconds(1) / dt));
lastFrame = min(session.numberOfVoltageFrames, ceil(displaySeconds(2) / dt));
frames = firstFrame:lastFrame;
timeSeconds = frames * dt;

voltage = double(session.readVoltageRows('subthreshold_voltage', ...
    qc(1):qc(end), [firstFrame lastFrame]));
voltage = voltage(qc - qc(1) + 1, :);
voltage = voltage(order, :);

whiskTrace = interpolateWhisking(session, false);

panel = struct();
panel.panel = 'E';
panel.cellCount = numel(qc);
panel.displaySeconds = displaySeconds;
panel.sortOrder = order;
panel.groupIndex = grouping.groupIndex;
panel.npyPositive = grouping.npyPositive;
panel.groupSizes = grouping.groupSizes;
panel.alignment = grouping.alignment;
panel.colourScaleLimits = [-1.25 1.25];

groupColours = umapGroupColours();
npyColour = npyColours();
figureHandle = figure('Name', 'Figure 6E', 'Color', 'w', ...
    'Position', [60 60 980 760]);

% Explicit axes positions, because the manuscript puts the two identity colour
% bars vertically at the right edge of the raster rather than below it.
ax1 = axes(figureHandle, 'Position', [0.09 0.865 0.77 0.095]);
plot(ax1, timeSeconds, whiskTrace(frames), 'k', 'LineWidth', 0.8);
ylabel(ax1, 'Whisking');
title(ax1, sprintf('UMAP-sorted spontaneous dynamics (n = %d cells): group 1 n = %d, group 2 n = %d', ...
    numel(qc), grouping.groupSizes(1), grouping.groupSizes(2)));
set(ax1, 'XTickLabel', [], 'Box', 'on');

ax2 = axes(figureHandle, 'Position', [0.09 0.495 0.77 0.355]);
imagesc(ax2, timeSeconds, 1:numel(qc), voltage);
caxis(ax2, panel.colourScaleLimits);
colormap(ax2, divergingColormap());
ylabel(ax2, 'Cell #');
set(ax2, 'YDir', 'normal', 'XTickLabel', []);
colorbar(ax2, 'Position', [0.935 0.495 0.015 0.355]);

ax3 = axes(figureHandle, 'Position', [0.09 0.09 0.77 0.39]);
hold(ax3, 'on');
for k = 1:numel(order)
    spikes = session.spikeFrameIndex{qc(order(k))};
    spikes = spikes(spikes >= firstFrame & spikes <= lastFrame);
    if isempty(spikes)
        continue
    end
    x = repmat(spikes(:)' * dt, 3, 1);
    y = repmat([k - 0.45; k + 0.45; NaN], 1, numel(spikes));
    plot(ax3, x(:), y(:), 'k-', 'LineWidth', 0.4);
end
hold(ax3, 'off');
ylim(ax3, [0 numel(order) + 1]);
ylabel(ax3, 'Cell #');
xlabel(ax3, 'Time (s)');
box(ax3, 'on');

% Two vertical identity strips beside the raster, cell order matching it:
% UMAP group on the left, NPY status on the right.
bars = zeros(numel(order), 2, 3);
for k = 1:numel(order)
    bars(k, 1, :) = groupColours(grouping.groupIndex(order(k)), :);
    bars(k, 2, :) = npyColour(grouping.npyPositive(order(k)) + 1, :);
end
ax4 = axes(figureHandle, 'Position', [0.875 0.09 0.045 0.39]);
image(ax4, 1:2, 1:numel(order), bars);
set(ax4, 'YDir', 'normal', 'YTick', [], 'XTick', [1 2], ...
    'XTickLabel', {'grp', 'NPY'}, 'FontSize', 9);

linkaxes([ax1 ax2 ax3], 'x');
xlim(ax1, displaySeconds);
ylim(ax4, [0.5 numel(order) + 0.5]);

saveFigure(figureHandle, outputDirectory, 'figure6E_umap_sorted_dynamics_raw.png');
end

% =========================================================================
% Panel F
% =========================================================================
function panel = plotUmapEmbedding(grouping, rasterGrouping, outputDirectory)
% Figure 6F. Source fig6_NPY_umap_confusionmat.m lines 96-116, figure(2).
%
% Which script produced the displayed panel F is not settled: the raster script
% draws the same scatter from an embedding built with a different
% chunk-boundary guard. Panel G definitely comes from the confusion-matrix
% script, so that grouping is plotted here, and the agreement between the two
% is reported so the ambiguity can be closed with data.
panel = struct();
panel.panel = 'F';
panel.guardHalfWidth = grouping.guardHalfWidth;
panel.coordinates = grouping.coordinates;
panel.groupIndex = grouping.groupIndex;
panel.groupSizes = grouping.groupSizes;
panel.alignment = grouping.alignment;
panel.cellCount = numel(grouping.groupIndex);

if ~isempty(rasterGrouping)
    agreement = mean(rasterGrouping.groupIndex == grouping.groupIndex);
    panel.guardComparison = struct( ...
        'rasterGuardHalfWidth', rasterGrouping.guardHalfWidth, ...
        'confusionGuardHalfWidth', grouping.guardHalfWidth, ...
        'groupAgreementFraction', agreement, ...
        'note', ['Uncertainty 2 in FIGURE_CODE_MAP.md: the raster script and ', ...
                 'the confusion-matrix script build the UMAP embedding with ', ...
                 'different chunk-boundary guards, so they can disagree. Both ', ...
                 'labellings are aligned to the NPY averages first, so this ', ...
                 'compares the clusterings, not the labelling convention.']);
end

colours = umapGroupColours();
figureHandle = figure('Name', 'Figure 6F', 'Color', 'w', ...
    'Position', [140 140 480 460]);
ax = axes(figureHandle);
hold(ax, 'on');
for iGroup = 1:2
    rows = grouping.groupIndex == iGroup;
    scatter(ax, grouping.coordinates(rows, 1), grouping.coordinates(rows, 2), ...
        30, colours(iGroup, :), 'filled');
end
hold(ax, 'off');
xlabel(ax, 'UMAP 1');
ylabel(ax, 'UMAP 2');
title(ax, 'Spontaneous spiking dynamics');
legend(ax, {sprintf('Group 1 (n = %d)', grouping.groupSizes(1)), ...
            sprintf('Group 2 (n = %d)', grouping.groupSizes(2))}, ...
    'Location', 'best', 'Box', 'off');
set(ax, 'XTick', [], 'YTick', [], 'FontSize', 11);
box(ax, 'on');

saveFigure(figureHandle, outputDirectory, 'figure6F_umap_embedding_raw.png');
end

% =========================================================================
% Panel G
% =========================================================================
function panel = plotConfusionMatrix(behaviorNwbFiles, firstSession, ...
    firstGrouping, outputDirectory)
% Figure 6G. Source fig6_NPY_umap_confusionmat.m lines 121-161, figure(99).
% Rows are NPY status, columns are UMAP group, pooled over both mice.
if isempty(which('UMAP'))
    panel = [];
    return
end

npyLabels = [];
groupLabels = [];
perSession = struct('nwbFile', {}, 'cellCount', {}, 'matrix', {}, 'alignment', {});

for iFile = 1:numel(behaviorNwbFiles)
    if iFile == 1
        session = firstSession;
        grouping = firstGrouping;
    else
        fprintf('  loading behavior session %s ...\n', behaviorNwbFiles{iFile});
        session = loadFigure6SessionFromNWB(char(behaviorNwbFiles{iFile}));
        assert(strcmp(session.sessionType, 'behavior'), 'Figure6:WrongBehaviorFile', ...
            '%s is not a behavior file.', behaviorNwbFiles{iFile});
        grouping = [];
    end
    if isempty(grouping)
        fprintf('  UMAP grouping for %s ...\n', session.nwbFile);
        grouping = umapGrouping(session, 5000);
    end
    if isempty(grouping)
        panel = [];
        return
    end
    npy = grouping.npyPositive(:);
    group = grouping.groupIndex(:);
    npyLabels = [npyLabels; npy]; %#ok<AGROW>
    groupLabels = [groupLabels; group]; %#ok<AGROW>
    perSession(end + 1) = struct( ...
        'nwbFile', session.nwbFile, ...
        'cellCount', numel(npy), ...
        'matrix', confusionCounts(npy, group), ...
        'alignment', grouping.alignment); %#ok<AGROW>
end

matrix = confusionCounts(npyLabels, groupLabels);

panel = struct();
panel.panel = 'G';
panel.cellCount = numel(npyLabels);
panel.mouseCount = numel(behaviorNwbFiles);
panel.matrix = matrix;
panel.rowLabels = {'NPY-', 'NPY+'};
panel.columnLabels = {'Group 1', 'Group 2'};
panel.rowTotals = sum(matrix, 2)';
panel.columnTotals = sum(matrix, 1);
panel.perSession = perSession;

figureHandle = figure('Name', 'Figure 6G', 'Color', 'w', ...
    'Position', [160 160 460 420]);
ax = axes(figureHandle);
imagesc(ax, matrix);
colormap(ax, flipud(gray(64)));
for iRow = 1:2
    for iColumn = 1:2
        text(ax, iColumn, iRow, sprintf('%d', matrix(iRow, iColumn)), ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold', ...
            'Color', textContrast(matrix(iRow, iColumn), max(matrix(:))));
    end
end
set(ax, 'XTick', 1:2, 'XTickLabel', panel.columnLabels, ...
    'YTick', 1:2, 'YTickLabel', panel.rowLabels, 'FontSize', 12);
xlabel(ax, 'UMAP group');
ylabel(ax, 'NPY');
title(ax, sprintf('n = %d cells, N = %d mice', panel.cellCount, panel.mouseCount));

saveFigure(figureHandle, outputDirectory, 'figure6G_npy_umap_confusion_raw.png');
end

function counts = confusionCounts(npyPositive, groupIndex)
counts = zeros(2);
for iRow = 1:2
    for iColumn = 1:2
        counts(iRow, iColumn) = ...
            sum(npyPositive == (iRow - 1) & groupIndex == iColumn);
    end
end
end

function colour = textContrast(value, maximum)
if value > 0.6 * maximum
    colour = [1 1 1];
else
    colour = [0 0 0];
end
end

% =========================================================================
% Panel H
% =========================================================================
function panel = plotSubtypeSpontaneousDynamics(session, stream, outputDirectory)
% Figure 6H. Source fig6_NPYpm_avgDynamics_raster.m lines 231-284, figure(2).
%
% This block re-defines the quiet and whisking masks with the source's MEDIAN
% rule, which is not the Gaussian-mixture rule stored in the NWB and used by
% panels I, J and K. Both are reported.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
isNpy = session.gfpPositive(qc);

runMask = stream.runMask;
quietMask = medianQuietMask(session);
whiskMask = medianWhiskMask(session, runMask);

voltage = stream.subtypeSubthreshold;
rate = smoothdata(stream.subtypeSpikeRateRejected, 2, 'gaussian', 32);
whiskTrace = interpolateWhisking(session, false);

displaySeconds = [290 320];
frames = max(1, floor(displaySeconds(1) / dt)):min(session.numberOfVoltageFrames, ...
    ceil(displaySeconds(2) / dt));
timeSeconds = frames * dt;

panel = struct();
panel.panel = 'H';
panel.cellCount = numel(qc);
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.displaySeconds = displaySeconds;
panel.subtypeVoltage = voltage(:, frames);
panel.subtypeRateHz = rate(:, frames);
panel.timeSeconds = timeSeconds;
panel.subtypeMeanRateHz = mean(stream.subtypeSpikeRateRejected, 2)';
panel.subtypeMeanVoltage = mean(voltage, 2)';
panel.quietMaskComparison = struct( ...
    'medianRuleQuietFrames', sum(quietMask), ...
    'storedMixtureQuietFrames', sum(session.quietMask), ...
    'agreementFraction', mean(quietMask == session.quietMask), ...
    'note', ['This panel''s state bar uses the source median-threshold rule ', ...
             'from fig6_NPYpm_avgDynamics_raster.m line 235. The stored ', ...
             'quiet_mask is the Gaussian-mixture rule used by panels I, J ', ...
             'and K. Both are reported; the plotted bar uses the median rule.']);

colours = npyColours();
figureHandle = figure('Name', 'Figure 6H', 'Color', 'w', ...
    'Position', [80 80 960 620]);

ax1 = subplot(3, 1, 1);
hold(ax1, 'on');
plot(ax1, timeSeconds, whiskTrace(frames), 'k', 'LineWidth', 0.8);
yBar = max(whiskTrace(frames)) * 1.1;
plotStateBar(ax1, timeSeconds, quietMask(frames), yBar, [0.5 0.5 0.5]);
plotStateBar(ax1, timeSeconds, whiskMask(frames), yBar, [0 0 0]);
hold(ax1, 'off');
ylabel(ax1, 'Whisking');
title(ax1, sprintf('Subtype average spontaneous dynamics (n = %d NPY-, %d NPY+)', ...
    panel.npyNegativeCount, panel.npyPositiveCount));
box(ax1, 'on');

ax2 = subplot(3, 1, 2);
plotSubtypes(ax2, timeSeconds, panel.subtypeVoltage, colours);
ylabel(ax2, 'Voltage (norm)');
ylim(ax2, [-1.5 0.5]);

ax3 = subplot(3, 1, 3);
plotSubtypes(ax3, timeSeconds, panel.subtypeRateHz, colours);
ylabel(ax3, 'Spike rate (Hz)');
xlabel(ax3, 'Time (s)');
legend(ax3, session.npyLabels, 'Location', 'northeast', 'Box', 'off');

linkaxes([ax1 ax2 ax3], 'x');
xlim(ax1, displaySeconds);

saveFigure(figureHandle, outputDirectory, 'figure6H_subtype_spontaneous_dynamics_raw.png');
end

function plotStateBar(ax, timeSeconds, mask, yValue, colour)
y = nan(size(timeSeconds));
y(mask) = yValue;
plot(ax, timeSeconds, y, 'Color', colour, 'LineWidth', 4);
end

% =========================================================================
% Panel I
% =========================================================================
function panel = plotOscillationTriggeredAverage(session, stream, outputDirectory)
% Figure 6I. Source fig6_NPY_osc_whisk_TrigAvg_ranksumStats.m lines 192-227.
%
% Triggers are the peaks of the 1.5-10 Hz band-passed population-average
% subthreshold voltage, restricted to the stored quiet mask and to the first
% 1.3e6 frames. The source builds per-cell spike-triggered averages and then
% averages them by subtype; a triggered average is linear, so the subtype means
% are triggered directly here, which is identical and far cheaper.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
isNpy = session.gfpPositive(qc);

filtered = bandpassPopulation(stream.populationSubthreshold, dt, [1.5 10], 4);
[~, peakFrames] = findpeaks(filtered, 'MinPeakHeight', 0);
peakCountRaw = numel(peakFrames);
peakFrames = intersect(peakFrames, find(session.quietMask));
peakCountQuiet = numel(peakFrames);
peakFrames(peakFrames > 1.3e6) = [];

nBack = 300;
nFront = 300;
lagSeconds = (-nBack:nFront) * dt;

voltageAverage = triggeredAverage(peakFrames, stream.subtypeSubthreshold, ...
    nBack, nFront);
rateAverage = triggeredAverage(peakFrames, stream.subtypeSpikeRate, nBack, nFront);

panel = struct();
panel.panel = 'I';
panel.cellCount = numel(qc);
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.eventCount = numel(peakFrames);
panel.eventCountStages = [peakCountRaw peakCountQuiet numel(peakFrames)];
panel.eventCountStageLabels = {'findpeaks', 'after quiet mask', 'after t <= 1.3e6'};
panel.lagSeconds = lagSeconds;
panel.subtypeVoltage = voltageAverage;
panel.subtypeRateHz = rateAverage;
panel.bandHz = [1.5 10];
panel.rateMatrixNote = ['The rate STA uses the raw spike matrix, not the ', ...
    '>15-rejected one, exactly as the source does. The author confirmed on ', ...
    '2026-09-19 that this is immaterial for an event-triggered average.'];

colours = npyColours();
figureHandle = figure('Name', 'Figure 6I', 'Color', 'w', ...
    'Position', [100 100 520 620]);

ax1 = subplot(2, 1, 1);
plotSubtypes(ax1, lagSeconds, voltageAverage, colours);
ylabel(ax1, 'Voltage (norm)');
title(ax1, sprintf('Osc triggered average, quiet (n = %d events)', panel.eventCount));

ax2 = subplot(2, 1, 2);
plotSubtypes(ax2, lagSeconds, rateAverage, colours);
ylabel(ax2, 'Spike rate (Hz)');
xlabel(ax2, 'Time (s)');
legend(ax2, session.npyLabels, 'Location', 'northeast', 'Box', 'off');

linkaxes([ax1 ax2], 'x');
xlim(ax1, [-0.4 0.4]);

saveFigure(figureHandle, outputDirectory, 'figure6I_oscillation_triggered_average_raw.png');
end

function filtered = bandpassPopulation(trace, dt, band, order)
% Inline equivalent of the lab helper butterworth_filt: a zero-phase
% Butterworth band-pass of the given order.
%
% Two details here are the source's, not conveniences, and both move the
% reported event count of panel I:
%
%   * The normalised band is computed as band*2/sampleRate with
%     sampleRate = 1/dt, exactly as butterworth_filt does. Writing the
%     algebraically equal band*2*dt differs in the last unit in the last place,
%     which perturbs the filter coefficients enough to add or drop one shallow
%     peak at MinPeakHeight 0.
%   * findpeaks is run on a single-precision trace, because line 194 of the
%     source reads
%       vAllAvg_lo = single(butterworth_filt(double(vAllAvg)',4,[1.5 10],1/dt)');
sampleRateHz = 1 / dt;
effectiveOrder = max(floor(order / 2), 1) * 2;
[b, a] = butter(effectiveOrder, band * 2 / sampleRateHz);
filtered = single(filtfilt(b, a, double(trace(:)))');
end

% =========================================================================
% Panel J
% =========================================================================
function panel = plotWhiskingTriggeredAverage(session, stream, outputDirectory)
% Figure 6J. Source fig6_NPY_osc_whisk_TrigAvg_ranksumStats.m lines 343-430.
%
% Note the deliberate asymmetry in the source: the voltage average is NaN-masked
% +/-140 frames around every recording-chunk boundary, the rate average is not.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
isNpy = session.gfpPositive(qc);
whisk = stream.whisk;

nBack = round(1.5 / dt);
nFront = round(1.5 / dt);
lagSeconds = (-nBack:nFront) * dt;

noiseMask = ones(1, session.numberOfVoltageFrames);
noiseMask(chunkBoundaryGuard(session.numberOfVoltageFrames, ...
    session.framesPerChunk, 140)) = NaN;

voltageAverage = triggeredAverage(whisk.onsetFrames, ...
    stream.subtypeSubthreshold .* noiseMask, nBack, nFront);
rateAverage = triggeredAverage(whisk.onsetFrames, ...
    stream.subtypeSpikeRateRejected, nBack, nFront);
whiskTrace = interpolateWhisking(session, false);
whiskAverage = triggeredAverage(whisk.onsetFrames, whiskTrace, nBack, nFront);

panel = struct();
panel.panel = 'J';
panel.cellCount = numel(qc);
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.eventCount = whisk.onsetCount;
panel.eventCountBeforeThinning = whisk.onsetCountBeforeThinning;
panel.selectedLambdaIndex = whisk.lambdaIndex;
panel.preOnsetStd = whisk.preOnsetStd;
panel.lagSeconds = lagSeconds;
panel.whiskingAverage = whiskAverage;
panel.subtypeVoltage = voltageAverage;
panel.subtypeRateHz = rateAverage;
panel.maskingNote = ['The voltage average is NaN-masked +/-140 frames around ', ...
    'every chunk boundary and the rate average is not, exactly as in the ', ...
    'source. Uncertainty 4 in FIGURE_CODE_MAP.md.'];

colours = npyColours();
figureHandle = figure('Name', 'Figure 6J', 'Color', 'w', ...
    'Position', [120 120 520 700]);

ax1 = subplot(3, 1, 1);
plot(ax1, lagSeconds, whiskAverage, 'k', 'LineWidth', 1.2);
ylabel(ax1, 'Whisking');
title(ax1, sprintf('Whisk triggered average (n = %d events, lambda %d of %d)', ...
    panel.eventCount, panel.selectedLambdaIndex, numel(whisk.preOnsetStd)));

ax2 = subplot(3, 1, 2);
plotSubtypes(ax2, lagSeconds, voltageAverage, colours);
ylabel(ax2, 'Voltage (norm)');

ax3 = subplot(3, 1, 3);
plotSubtypes(ax3, lagSeconds, rateAverage, colours);
ylabel(ax3, 'Spike rate (Hz)');
xlabel(ax3, 'Peri-whisk time (s)');
legend(ax3, session.npyLabels, 'Location', 'northeast', 'Box', 'off');

linkaxes([ax1 ax2 ax3], 'x');
xlim(ax1, [-1.5 1.5]);

saveFigure(figureHandle, outputDirectory, 'figure6J_whisking_triggered_average_raw.png');
end

% =========================================================================
% Panel K
% =========================================================================
function panel = plotStateDependentExcitability(session, stream, outputDirectory)
% Figure 6K. Source fig6_NPYpm_pSpkV_whisk.m lines 354-426, figure(11).
%
% "Whisking" here is the 0-635 ms window after a selected whisking onset, not
% the whole whisking epoch. The four curves pool counts across cells before
% dividing; they are not means of per-cell ratios.
dt = session.framePeriodSeconds;
qc = session.figureCellIndices;
isNpy = session.gfpPositive(qc);

curves = zeros(4, numel(stream.voltageCenters));
labels = cell(4, 1);
stateLabels = stream.stateLabels;
typeLabels = session.npyLabels;
for iState = 1:2
    for iType = 1:2
        rows = (isNpy == (iType - 1));
        row = iType + (iState - 1) * 2;
        spikeCounts = sum(stream.spikeVoltageHistogram(rows, :, iState), 1);
        allCounts = sum(stream.voltageHistogram(rows, :, iState), 1);
        curves(row, :) = spikeCounts ./ allCounts / dt;
        labels{row} = sprintf('%s, %s', typeLabels{iType}, stateLabels{iState});
    end
end

panel = struct();
panel.panel = 'K';
panel.cellCount = numel(qc);
panel.npyPositiveCount = sum(isNpy);
panel.npyNegativeCount = sum(~isNpy);
panel.voltageCenters = stream.voltageCenters;
panel.voltageEdges = stream.voltageEdges;
panel.spikeProbabilityHz = curves;
panel.curveLabels = labels;
panel.quietFrameCount = stream.quietFrameCount;
panel.whiskOnsetFrameCount = stream.whiskOnsetFrameCount;
panel.selectedLambdaIndex = stream.whisk.lambdaIndex;
panel.whiskOnsetEventCount = stream.whisk.onsetCount;

colours = npyColours();
figureHandle = figure('Name', 'Figure 6K', 'Color', 'w', ...
    'Position', [140 140 520 470]);
ax = axes(figureHandle);
hold(ax, 'on');
for iType = 1:2
    plot(ax, panel.voltageCenters, curves(iType, :), '-', ...
        'Color', colours(iType, :), 'LineWidth', 2);
end
for iType = 1:2
    plot(ax, panel.voltageCenters, curves(iType + 2, :), '--', ...
        'Color', colours(iType, :), 'LineWidth', 2);
end
hold(ax, 'off');
xlim(ax, [0 1]);
xlabel(ax, 'Voltage (norm)');
ylabel(ax, 'P(spk | V) (Hz)');
title(ax, sprintf('State-dependent excitability (n = %d cells)', panel.cellCount));
legend(ax, labels, 'Location', 'northwest', 'Box', 'off');
box(ax, 'on');
set(ax, 'FontSize', 11);

saveFigure(figureHandle, outputDirectory, 'figure6K_state_dependent_excitability_raw.png');
end

% =========================================================================
% Shared plotting helpers
% =========================================================================
function colours = npyColours()
% The source builds its NPY pair from colorcet('cbtd1'), a third-party map:
% lc = max([1 0 0; cmCB(1,:)] - .2, 0), a dark red and a dark cyan. These are
% fixed equivalents so the reproduction has no third-party colour dependency.
colours = [0.8 0 0; 0 0.55 0.65];
end

function colours = umapGroupColours()
% The source uses the two ends of colorcet('cbd2'), rendered yellow and blue in
% the manuscript. Cosmetic.
colours = [0.93 0.78 0.13; 0.14 0.42 0.80];
end

function map = divergingColormap()
% Stand-in for colorcet('d1'), blue-white-red, for the panel E voltage heatmap.
n = 128;
half = linspace(0, 1, n)';
map = [ [half * 0.9 + 0.1, half, ones(n, 1)]; ...
        flipud([ones(n, 1), half, half * 0.9 + 0.1]) ];
map(n, :) = [1 1 1];
map(n + 1, :) = [1 1 1];
end

function plotSubtypes(ax, x, y, colours)
% y is [2 x nSample] or [nSample x 2]; orient to match x.
if size(y, 2) ~= numel(x) && size(y, 1) == numel(x)
    y = y';
end
hold(ax, 'on');
for iType = 1:size(y, 1)
    plot(ax, x, y(iType, :), 'LineWidth', 1.5, 'Color', colours(iType, :));
end
hold(ax, 'off');
box(ax, 'on');
set(ax, 'FontSize', 10);
end

function saveFigure(figureHandle, outputDirectory, fileName)
% Matches reproduce_figure3/4/5_from_nwb.m: print at -r300, then leave the
% figure open and visible so the panels can be inspected after the run.
set(figureHandle, 'PaperPositionMode', 'auto');
target = fullfile(outputDirectory, fileName);
writeWithRetry(@() print(figureHandle, target, '-dpng', '-r300'), target);
set(figureHandle, 'Visible', 'on');
drawnow;
end

function saveResults(resultsPath, results)
% save() resolves variable names in the calling workspace, so it must be invoked
% here rather than inside an anonymous function, which has its own scope.
for attempt = 1:3
    try
        save(resultsPath, 'results', '-v7.3', '-nocompression');
        return
    catch err
        if attempt == 3
            rethrow(err);
        end
        fprintf('  write retry %d for %s (%s)\n', attempt, resultsPath, err.identifier);
        pause(1);
    end
end
end

function writeWithRetry(writeFcn, target)
% Dropbox-synced folders occasionally hold a brief lock on a new file.
for attempt = 1:3
    try
        writeFcn();
        return
    catch err
        if attempt == 3
            rethrow(err);
        end
        fprintf('  write retry %d for %s (%s)\n', attempt, target, err.identifier);
        pause(1);
    end
end
end
