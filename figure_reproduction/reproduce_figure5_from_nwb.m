function results = reproduce_figure5_from_nwb(optoNwbFile, connNwbFile, ...
    outputDirectory, panels)
% REPRODUCE_FIGURE5_FROM_NWB Recreate the data panels of manuscript Figure 5.
%
%   RESULTS = REPRODUCE_FIGURE5_FROM_NWB(OPTONWBFILE, CONNNWBFILE, OUTPUTDIRECTORY, PANELS)
%   reads the two M-YQ0201-27 NWB files Figure 5 needs, writes raw panels and a
%   reference MAT file to OUTPUTDIRECTORY, and returns the numerical results. It
%   does not edit or execute the legacy scripts.
%
%   Every Figure 5 panel uses M-YQ0201-27 (confirmed by the author on
%   2026-09-18). Two files are needed because the figure spans two protocols:
%
%       A, B, C      M_YQ0201_27_sparseOpto.nwb        intrinsic properties
%       D, E, F, G   M_YQ0201_27_sparsePulseHad.nwb    connectivity
%
%   The two-mouse numbers quoted in the Figure 5 narrative come from
%   sup_fig5_connMap_stats_norm.m, not from the panel scripts, and are out of
%   scope here.
%
%   All four arguments are optional and follow the same convention as
%   reproduce_figure2/3/4_from_nwb.m: an omitted or empty argument falls back to
%   its default, so REPRODUCE_FIGURE5_FROM_NWB() regenerates the whole figure
%   from the documented session paths.
%
%     optoNwbFile      default M_YQ0201_27_sparseOpto.nwb
%     connNwbFile      default M_YQ0201_27_sparsePulseHad.nwb
%     outputDirectory  default figure5_nwb_reproduction/ beside this file
%     panels           default all
%
%   PANELS accepts a panel letter, a comma-separated character vector, or a cell
%   array, for example 'C', 'D,E', or {'F','G'}. Only the NWB sessions the
%   request needs are opened, so 'D,E,F,G' never touches the sparseOpto file and
%   'A,B,C' never touches the connectivity file.
%
%   Each panel is drawn in its own named, visible figure window, printed to PNG
%   at 300 dpi, and left open for inspection, as in the other reproduction
%   scripts.
%
%   Examples
%   --------
%     reproduce_figure5_from_nwb();                       % everything, 70-95 s
%     reproduce_figure5_from_nwb([], [], [], 'C');        % one panel, ~60 s
%     reproduce_figure5_from_nwb([], [], [], 'D,E,F,G');  % connectivity, ~10 s
%     reproduce_figure5_from_nwb(optoFile, connFile, 'C:\tmp\fig5');
%
%   Any of A, B or C runs the one streaming pass over the sparseOpto recording,
%   so 'C' alone costs about as much as 'A,B,C'. Panel C needs it because its
%   autocorrelogram triggers use the per-chunk blue-neighbour mask.
%
%   A subset run writes a reference MAT holding only the requested panels, so
%   point OUTPUTDIRECTORY elsewhere to keep the full reference file.
%
%   Raw panel arrangement and styling are diagnostic, not a replacement for the
%   externally assembled manuscript figure. In particular the source script puts
%   the spontaneous curve above the evoked one in panel A; the manuscript prints
%   evoked on top. Only the arrangement differs.
%
%   Source provenance
%   -----------------
%   Line numbers refer to the legacy scripts as of 2026-09-26.
%   Panels A, B, C   fig5_ephys_vs_spiking.m
%                      A: lines 367-438 with idxPlt = 2 (optical rheobase)
%                      B: lines 735-880 with idxPlt = 4 (spike-rate adaptation);
%                         dvdt_all line 737, median tQuiet line 758, plot 798-880
%                      C: lines 541-602 with idxPlt = 1 (after-depolarization);
%                         selfStaRate_q lines 472-473
%   Panels D, E      fig5_synapConn_gjConn_spontAct_ctrlRes.m lines 330-413
%                    (compute) and 452-490 (plot)
%   Panels F, G      fig5_synapConn_gjConn_spontAct_ctrlRes.m lines 520-705
%   Full panel-level evidence is in FIGURE_CODE_MAP.md; execution notes and the
%   preserved-behavior list are in FIGURE5_PHASE1_HANDOFF.md.
%
%   Preserved source behavior
%   -------------------------
%   These are copied deliberately and must not be "cleaned up":
%     * Quintile edges are linspace(prctile(p,1), prctile(p,99), 6), so cells
%       outside the 1st-99th percentile land in bin 0 and are dropped. That is
%       why n = 280 rather than the 286 cells that passed QC.
%     * Panel A pools counts before dividing, sum(spk)/sum(all) per quintile.
%       It is not a mean of per-cell ratios.
%     * Panel A's ramp mask keeps only the lower half of the ramp
%       (rampAlpha = 0.5) and excludes step transitions.
%     * Panel B's dV/dt is an 8-frame centred difference padded with 4 zeros at
%       each end, then circularly shifted by nShift = 4, and its marginal is
%       taken only over 0.2 < V < 0.8.
%     * Panels D-G average with spike-count weights, except panel D's IPSP
%       waveform which is a plain nanmean.
%     * Panels E and G are "measured relative to distance-matched unconnected
%       controls" (caption wording): each quintile average is regressed onto
%       its distance-matched unconnected-pair average, fitted on the pre-spike
%       samples only (idxFit = 1:nBack-1) and then applied to the whole trace.
%     * Control pairs are rebuilt per quintile and per presynaptic cell.
%     * Panels D-G all restrict the connected pairs to the 60-400 um band, in
%       the quintile edges, the mask and the reported n. For D and E this is the
%       source's [distfix 2026-09-24] change; before it, D and E used every
%       connected pair and applied the band only to the controls.
%     * Spontaneous spikes for panels A and C, including panel C's
%       autocorrelogram triggers, are taken only in chunks where no QC cell
%       within 80 um, itself included, received blue.
%     * The six baseline conventions across D-G all differ; see each panel.
%
%   Documented deviations from the legacy scripts
%   ---------------------------------------------
%   1. Motion-corrected spike times. The source scripts drop spikes inside
%      high-motion windows via spike_times_motion_correction_aggressive, which
%      needs mcTrace_all. That variable is not in the NWB. The author confirmed
%      on 2026-09-18 that this step can be ignored, so this script uses the
%      stored raw spike times. results.deviations records it. In
%      fig5_synapConn_gjConn_spontAct_ctrlRes.m the no-blue spikes are already
%      built from the raw spk_t_mat (line 309), so there the correction reaches
%      only selfHadStaFN, which deviation 3 replaces anyway.
%   2. Panel B's quiet mask. The source re-defines tQuiet with a median
%      threshold at line 758, while the stored quiet_mask is the Gaussian
%      mixture boundary used by panels A and C. Both are computed here and both
%      are returned; the plotted curve uses the source's median rule. See
%      results.panelB.quietMaskComparison.
%   3. Panel F's top trace. The source plots mean(selfHadStaFN,1), the
%      self-STA of all Hadamard-evoked spikes (line 660). The NWB stores, on the
%      diagonal of hadamard_cross_spike_sta, a self-STA built from the
%      Hadamard-orthogonal subset instead. This script plots the plain mean of
%      that stored diagonal over the QC cells and records the difference rather
%      than re-deriving the source variant.
%   4. The >20-simultaneous-spike rejection. fig5_ephys_vs_spiking.m line 95
%      zeroes every frame in which more than 20 cells spike. The stored
%      spike_times and spike_times_no_blue come from the raw spike matrix, and
%      the matching dataPrepsForMAT.m line 565 is commented out, so the
%      rejection is not applied. It travels with deviation 1 but is a separate
%      step.
%   5. The neighbour pool behind spike_times_no_blue. The source builds its
%      frame-wise no-blue spikes from the 80 um neighbours among the QC cells;
%      dataPrepsForMAT.m built the stored spike_times_no_blue from all 320
%      cells. The chunk-level mask computed here uses the QC pool, as the
%      source does, so only the frame-wise exclusion inside the stored spike
%      times differs. Measured on 2026-09-26 by re-deriving that exclusion from
%      the raw spike times with the QC pool: 27 of 286 QC cells change, their
%      no-blue spikes rise by 1.4% in total, the panel A spontaneous curves
%      move by at most 4.6% of their maximum, panel B by 0.6% and the panel C
%      autocorrelogram peaks by at most 4.3%. The stored spike times are used
%      as they are. For panels E and G the analogous question, blueOn pooled
%      over QC cells (source line 300) or all 320 cells, makes no difference:
%      the two pools mark identical frames in this session.
%
%   Outputs
%   -------
%   One PNG per panel plus figure5_nwb_reference_results.mat in
%   OUTPUTDIRECTORY, and a RESULTS struct with one field per panel, a
%   sampleSizes field carrying the caption counts, and a deviations field.
%
%   See also LOADFIGURE5SESSIONFROMNWB, REPRODUCE_FIGURE4_FROM_NWB.

% Placeholders: pass the paths to your local copies of these files instead.
if nargin < 1 || isempty(optoNwbFile)
    optoNwbFile = 'path/to/M_YQ0201_27_sparseOpto.nwb';
end
if nargin < 2 || isempty(connNwbFile)
    connNwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';
end
if nargin < 3 || isempty(outputDirectory)
    outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
        'figure5_nwb_reproduction');
end
if nargin < 4
    panels = [];
end
selected = normalizePanelSelection(panels);
if exist(outputDirectory, 'dir') ~= 7
    mkdir(outputDirectory);
end

needOpto = any(ismember({'A', 'B', 'C'}, selected));
needConn = any(ismember({'D', 'E', 'F', 'G'}, selected));
if needOpto
    assert(exist(optoNwbFile, 'file') == 2, 'Figure5:MissingOptoFile', ...
        'Panels A, B and C need the sparseOpto NWB file, not found: %s', optoNwbFile);
end
if needConn
    assert(exist(connNwbFile, 'file') == 2, 'Figure5:MissingConnFile', ...
        'Panels D, E, F and G need the sparsePulseHad NWB file, not found: %s', ...
        connNwbFile);
end

results = struct();
results.generatedOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');
results.panelsRequested = selected;
results.deviations = documentedDeviations();
results.sampleSizes = struct();

nBin = 5;   % quintiles throughout Figure 5

% ------------------------------------------------------------ panels A, B, C
if needOpto
    fprintf('loading sparseOpto session ...\n');
    opto = loadFigure5SessionFromNWB(char(optoNwbFile));
    assert(strcmp(opto.sessionType, 'sparseOpto'), 'Figure5:WrongOptoFile', ...
        '%s is not a sparseOpto file.', optoNwbFile);
    results.optoSession = sessionSummary(opto);

    % Panel C needs the streaming pass too: its autocorrelogram triggers are
    % the source spk_t_quiet, which depend on the per-chunk blue-neighbour mask
    % and therefore on the stimulation of every cell.
    fprintf('streaming voltage and stimulation for panels A, B and C ...\n');
    stream = accumulateOptoStatistics(opto);
    results.optoStream = streamSummary(stream);

    if ismember('A', selected)
        fprintf('panel A ...\n');
        results.panelA = plotExcitabilityByRheobase(opto, stream, nBin, outputDirectory);
        results.sampleSizes.panelA_cells = results.panelA.cellCount;
    end
    if ismember('B', selected)
        fprintf('panel B ...\n');
        results.panelB = plotSpikeProbabilityByAdaptation(opto, stream, nBin, outputDirectory);
        results.sampleSizes.panelB_cells = results.panelB.cellCount;
    end
    if ismember('C', selected)
        fprintf('panel C ...\n');
        results.panelC = plotAutocorrelogramByADP(opto, stream, nBin, outputDirectory);
        results.sampleSizes.panelC_cells = results.panelC.cellCount;
    end
end

% --------------------------------------------------------- panels D, E, F, G
if needConn
    fprintf('loading sparsePulseHad session ...\n');
    conn = loadFigure5SessionFromNWB(char(connNwbFile));
    assert(strcmp(conn.sessionType, 'sparsePulseHad'), 'Figure5:WrongConnFile', ...
        '%s is not a sparsePulseHad file.', connNwbFile);
    results.connSession = sessionSummary(conn);

    if any(ismember({'D', 'E'}, selected))
        fprintf('panels D and E: IPSP-amplitude quintiles ...\n');
        inhibition = quintileConnectivityAverages(conn, 'synaptic', nBin);
        results.sampleSizes.panelDE_pairs = inhibition.pairCount;
        if ismember('D', selected)
            results.panelD = plotIpspWaveforms(conn, inhibition, outputDirectory);
        end
        if ismember('E', selected)
            results.panelE = plotSpontaneousResiduals(inhibition, 'E', ...
                'ipsp_amplitude', outputDirectory);
        end
    end

    if any(ismember({'F', 'G'}, selected))
        fprintf('panels F and G: spikelet-amplitude quintiles ...\n');
        junction = quintileConnectivityAverages(conn, 'gap_junction', nBin);
        results.sampleSizes.panelFG_pairs = junction.pairCount;
        if ismember('F', selected)
            results.panelF = plotSpikeletWaveforms(junction, outputDirectory);
        end
        if ismember('G', selected)
            results.panelG = plotSpontaneousResiduals(junction, 'G', ...
                'spikelet_amplitude', outputDirectory);
        end
    end
end

resultsPath = fullfile(outputDirectory, 'figure5_nwb_reference_results.mat');
saveResults(resultsPath, results);
fprintf('Figure 5 panel(s) %s written to %s\n', strjoin(selected, ','), outputDirectory);
end

% =========================================================================
function selectedPanels = normalizePanelSelection(panels)
valid = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};
if isempty(panels) || (ischar(panels) && strcmpi(strtrim(panels), 'all'))
    selectedPanels = valid;
    return
end
if ischar(panels)
    tokens = strtrim(strsplit(panels, ','));
elseif iscellstr(panels) %#ok<ISCLSTR> R2019b compatible
    tokens = strtrim(panels);
else
    error('Figure5:PanelSelection', ...
        'PANELS must be a character vector or a cell array of character vectors.');
end
tokens = upper(tokens(~cellfun(@isempty, tokens)));
unknown = setdiff(tokens, valid);
assert(isempty(unknown), 'Figure5:UnknownPanel', ...
    'Unknown panel selector(s): %s. Valid: %s.', strjoin(unknown, ', '), ...
    strjoin(valid, ', '));
selectedPanels = valid(ismember(valid, tokens));
end

function deviations = documentedDeviations()
deviations = { ...
    ['Motion-corrected spike times are not in the NWB (mcTrace_all was not ', ...
     'carried into the prepared MAT). The source applies ', ...
     'spike_times_motion_correction_aggressive; the author waived this on ', ...
     '2026-09-18, so stored raw spike times are used.']; ...
    ['Panel B uses the source median-threshold quiet rule, recomputed here ', ...
     'from whisking_motion. The stored quiet_mask is the Gaussian-mixture ', ...
     'rule used by panels A and C. Both are reported.']; ...
    ['Panel F top trace is the mean over QC cells of the stored diagonal of ', ...
     'hadamard_cross_spike_sta, which dataPrepsForMAT.m built from the ', ...
     'Hadamard-orthogonal spike subset, not from the source selfHadStaFN ', ...
     'which used all Hadamard-evoked spikes.']; ...
    ['The >20-simultaneous-spike rejection (fig5_ephys_vs_spiking.m line 95) ', ...
     'is not applied: the stored spike times come from the raw spike matrix ', ...
     'and dataPrepsForMAT.m line 565 is commented out.']; ...
    ['Stored spike_times_no_blue excludes frames with blue on any of all 320 ', ...
     'cells within 80 um; the source uses only the QC cells. The chunk-level ', ...
     'blue-neighbour mask computed here uses the QC pool, as the source does. ', ...
     'Re-deriving the frame-wise exclusion with the QC pool changes 27 of 286 ', ...
     'QC cells (+1.4% no-blue spikes) and moves panels A-C by at most 4.6%. ', ...
     'For panels E and G, blueOn over QC cells and over all 320 cells mark ', ...
     'identical frames in this session, so the nNoBlueSpk weights are exact.']};
end

function summary = sessionSummary(session)
summary = struct( ...
    'nwbFile', session.nwbFile, ...
    'sessionType', session.sessionType, ...
    'numberOfCells', session.numberOfCells, ...
    'figureCellCount', numel(session.figureCellIndices), ...
    'numberOfVoltageFrames', session.numberOfVoltageFrames, ...
    'framesPerChunk', session.framesPerChunk, ...
    'numberOfChunks', session.numberOfChunks, ...
    'framePeriodSeconds', session.framePeriodSeconds, ...
    'pixelSizeMicrometers', session.pixelSizeMicrometers);
end

function summary = streamSummary(stream)
summary = struct( ...
    'voltageEdges', stream.voltageEdges, ...
    'dvdtEdges', stream.dvdtEdges, ...
    'blueStepLevels', stream.blueStepLevels, ...
    'analysisFrameCount', stream.analysisFrameCount, ...
    'rampFrameCount', stream.rampFrameCount, ...
    'stimulatedChunksPerCell', stream.stimulatedChunkCount);
end

% =========================================================================
function stream = accumulateOptoStatistics(session)
% One streaming pass over the sparseOpto recording that builds every quantity
% panels A and B need, plus panel C's autocorrelogram triggers. Blocks are
% aligned to whole recording chunks, which is required because the
% blue-stimulation masks and the trial averages are both defined per chunk.
%
% Accumulated, all on the full 320-cell axis:
%   voltageHistogram        [nCells x nV]      source vHistAll
%   voltageSpikeHistogram   [nCells x nV]      source vSpkHistAll
%   rampHistogram           [nCells x nV]      source vBlueHistAll
%   rampSpikeHistogram      [nCells x nV]      source vBlueSpkHistAll
%   jointHistogram          [nCells x nVj x nD]    source vHistAll  (panel B)
%   jointSpikeHistogram     [nCells x nVj x nD]    source vSpkHistAll (panel B)
%   trialVoltageSum         [nCells x framesPerChunk]
%   trialSpikeSum           [nCells x framesPerChunk]
%   quietSpikeFrames        {nCells x 1}       source spk_t_quiet (panel C)

nCells = session.numberOfCells;
nFrames = session.numberOfVoltageFrames;
framesPerChunk = session.framesPerChunk;
dt = session.framePeriodSeconds;

% ---------------------------------------------------------------- binning
stream.voltageEdges = linspace(-1, 1.5, 51);          % source vEdge, panels A
stream.voltageCenters = binCenters(stream.voltageEdges);
nV = numel(stream.voltageEdges) - 1;

stream.jointVoltageEdges = linspace(-1, 1.5, 11);     % source vEdge, panel B
stream.jointVoltageCenters = binCenters(stream.jointVoltageEdges);
nVj = numel(stream.jointVoltageEdges) - 1;

stream.dvdtEdges = linspace(-0.5, 1, 51);             % source dvdtEdge
stream.dvdtCenters = binCenters(stream.dvdtEdges);
nD = numel(stream.dvdtEdges) - 1;

% ------------------------------------------------------- run and time masks
% Source lines 255-258: tRun dilates any nonzero running speed by +/-1000
% frames. conv with a 2001-frame box, then > 0.
runKernel = ones(1, 2001);
runMask = conv(double(abs(session.runningSpeed) > 1e-3), runKernel, 'same') > 0;

% Source lines 268-274 (panels A and C): keep frames 5000:2e6, drop a +/-100
% frame guard at every chunk boundary, then require quiet and not running.
guardA = chunkBoundaryGuard(nFrames, framesPerChunk, 100);
maskA = false(1, nFrames);
maskA(5000:min(2e6, nFrames)) = true;
maskA(guardA) = false;
maskA = maskA & session.quietMask & ~runMask;

% Source lines 753-765 (panel B): a different window, a +/-140 frame guard, and
% the median whisking rule instead of the stored mixture mask.
medianQuiet = medianQuietMask(session);
guardB = chunkBoundaryGuard(nFrames, framesPerChunk, 140);
maskB = false(1, nFrames);
maskB(5000:min(round(2000 / dt), nFrames)) = true;
maskB(guardB) = false;
maskB = maskB & ~runMask & medianQuiet;

stream.medianQuietMask = medianQuiet;
stream.storedQuietMask = session.quietMask;
stream.runMask = runMask;

% --------------------------------------------- blue step levels, source bVStep
stream.blueStepLevels = blueStepLevels(session);
smallestStep = stream.blueStepLevels(1);
largestStep = stream.blueStepLevels(end);
rampAlpha = 0.5;   % source line 294

% ------------------------------------------------------ neighbour index, 80 um
% Source lines 176-179 build distMat from cell_coord after the idxUse subset,
% so the neighbour pool is the QC cells only. Each cell is kept in its own pool
% so the rows of cells outside QC, which no panel uses, stay defined.
distanceMicrometers = session.cellDistancePixels * session.pixelSizeMicrometers;
inPool = session.passedQC(:)';
neighbourWithin80 = distanceMicrometers < 80 & ...
    (repmat(inPool, session.numberOfCells, 1) | eye(session.numberOfCells));

% ------------------------------------------------------------- accumulators
stream.voltageHistogram = zeros(nCells, nV);
stream.voltageSpikeHistogram = zeros(nCells, nV);
stream.rampHistogram = zeros(nCells, nV);
stream.rampSpikeHistogram = zeros(nCells, nV);
stream.jointHistogram = zeros(nCells, nVj, nD);
stream.jointSpikeHistogram = zeros(nCells, nVj, nD);
stream.trialVoltageSum = zeros(nCells, framesPerChunk);
stream.trialSpikeSum = zeros(nCells, framesPerChunk);
stream.stimulatedChunkCount = zeros(nCells, 1);
stream.analysisFrameCount = zeros(nCells, 1);
stream.rampFrameCount = zeros(nCells, 1);
quietSpikeFrames = cell(nCells, 1);
for iCell = 1:nCells
    quietSpikeFrames{iCell} = zeros(1, 0);
end

% Per-cell logical spike maps are rebuilt per block from the frame indices.
spikeFrames = session.spikeFrameIndex;
noBlueFrames = session.spikeFrameIndexNoBlue;

% --------------------------------------------------------------- block loop
chunksPerBlock = max(1, floor(session.voltageChunkFrames * 4 / framesPerChunk));
halo = 8;   % enough for the 8-frame centred difference and the nShift = 4 shift
firstChunk = 1;
while firstChunk <= session.numberOfChunks
    lastChunk = min(firstChunk + chunksPerBlock - 1, session.numberOfChunks);
    firstFrame = (firstChunk - 1) * framesPerChunk + 1;
    lastFrame = lastChunk * framesPerChunk;

    haloStart = max(1, firstFrame - halo);
    haloStop = min(nFrames, lastFrame + halo);
    leadIn = firstFrame - haloStart;

    voltage = double(session.readVoltageRows('subthreshold_voltage', ...
        1:nCells, [haloStart haloStop]));
    blue = double(session.readBlueRows(1:nCells, [haloStart haloStop]));

    % dV/dt on the haloed block, then trim. Source line 737 pads with 4 zeros
    % at each end of the whole recording; inside the recording the halo makes
    % the centred difference exact.
    dvdt = centredDifference8(voltage);
    dvdt = circshift(dvdt, 4, 2);          % source nShift = 4

    coreCols = leadIn + (1:(lastFrame - firstFrame + 1));
    voltage = voltage(:, coreCols);
    blueCore = blue(:, coreCols);
    dvdt = dvdt(:, coreCols);
    blockFrames = firstFrame:lastFrame;
    nBlockFrames = numel(blockFrames);

    % Ramp mask needs the first difference of blue, which the halo supplies.
    blueDiff = [zeros(nCells, 1) diff(blue, 1, 2)];
    blueDiff = blueDiff(:, coreCols);

    % Per-chunk blue presence inside this block.
    nBlockChunks = lastChunk - firstChunk + 1;
    blueByChunk = reshape(blueCore, nCells, framesPerChunk, nBlockChunks);
    cellHadBlue = squeeze(any(blueByChunk > 0, 2));
    if nBlockChunks == 1
        cellHadBlue = cellHadBlue(:);
    end

    blockMaskA = maskA(blockFrames);
    blockMaskB = maskB(blockFrames);

    for iCell = 1:nCells
        % ------------------------------------------- blueStimMask_all, lines 190-192
        % A chunk is usable for cell iCell only if no QC cell within 80 um,
        % itself included, received blue anywhere in that chunk.
        neighbours = neighbourWithin80(iCell, :);
        chunkUsable = ~any(cellHadBlue(neighbours, :), 1);
        usable = repelem(chunkUsable, 1, framesPerChunk);

        % ------------------------------------------------ panel A spontaneous
        cellMaskA = blockMaskA & usable;
        if any(cellMaskA)
            v = voltage(iCell, cellMaskA);
            stream.voltageHistogram(iCell, :) = stream.voltageHistogram(iCell, :) + ...
                histcounts(v, stream.voltageEdges);
            stream.analysisFrameCount(iCell) = stream.analysisFrameCount(iCell) + ...
                sum(cellMaskA);

            spikeHere = blockSpikeMask(noBlueFrames{iCell}, firstFrame, nBlockFrames);
            spikeInMask = spikeHere & cellMaskA;
            if any(spikeInMask)
                stream.voltageSpikeHistogram(iCell, :) = ...
                    stream.voltageSpikeHistogram(iCell, :) + ...
                    histcounts(voltage(iCell, spikeInMask), stream.voltageEdges);
                quietSpikeFrames{iCell} = [quietSpikeFrames{iCell} ...
                    blockFrames(spikeInMask)];
            end
        end

        % ------------------------------------------------------ panel A ramp
        % Source line 302: rising, but slower than the smallest step jump, and
        % only over the lower rampAlpha of the largest step.
        rampMask = blueDiff(iCell, :) > 0 & blueDiff(iCell, :) < smallestStep & ...
            blueCore(iCell, :) < rampAlpha * largestStep;
        if any(rampMask)
            stream.rampHistogram(iCell, :) = stream.rampHistogram(iCell, :) + ...
                histcounts(voltage(iCell, rampMask), stream.voltageEdges);
            stream.rampFrameCount(iCell) = stream.rampFrameCount(iCell) + sum(rampMask);

            blueSpike = blockSpikeMask(spikeFrames{iCell}, firstFrame, nBlockFrames) & ...
                blueCore(iCell, :) > 0;
            rampSpike = blueSpike & rampMask;
            if any(rampSpike)
                stream.rampSpikeHistogram(iCell, :) = ...
                    stream.rampSpikeHistogram(iCell, :) + ...
                    histcounts(voltage(iCell, rampSpike), stream.voltageEdges);
            end
        end

        % ---------------------------------------------------------- panel B
        cellMaskB = blockMaskB;
        if any(cellMaskB)
            jointAll = histcounts2(voltage(iCell, cellMaskB), dvdt(iCell, cellMaskB), ...
                stream.jointVoltageEdges, stream.dvdtEdges);
            stream.jointHistogram(iCell, :, :) = ...
                squeeze(stream.jointHistogram(iCell, :, :)) + jointAll;

            spikeB = blockSpikeMask(noBlueFrames{iCell}, firstFrame, nBlockFrames) & cellMaskB;
            if any(spikeB)
                jointSpike = histcounts2(voltage(iCell, spikeB), dvdt(iCell, spikeB), ...
                    stream.jointVoltageEdges, stream.dvdtEdges);
                stream.jointSpikeHistogram(iCell, :, :) = ...
                    squeeze(stream.jointSpikeHistogram(iCell, :, :)) + jointSpike;
            end
        end

        % --------------------------------- trial averages, source lines 237-250
        % Averaged over the chunks in which this cell itself was stimulated.
        stimChunks = find(cellHadBlue(iCell, :));
        if ~isempty(stimChunks)
            voltageStack = reshape(voltage(iCell, :), framesPerChunk, nBlockChunks);
            spikeStack = reshape( ...
                blockSpikeMask(spikeFrames{iCell}, firstFrame, nBlockFrames), ...
                framesPerChunk, nBlockChunks);
            stream.trialVoltageSum(iCell, :) = stream.trialVoltageSum(iCell, :) + ...
                sum(voltageStack(:, stimChunks), 2)';
            stream.trialSpikeSum(iCell, :) = stream.trialSpikeSum(iCell, :) + ...
                sum(double(spikeStack(:, stimChunks)), 2)';
            stream.stimulatedChunkCount(iCell) = ...
                stream.stimulatedChunkCount(iCell) + numel(stimChunks);
        end
    end

    fprintf('  chunks %d-%d of %d\n', firstChunk, lastChunk, session.numberOfChunks);
    firstChunk = lastChunk + 1;
end

stream.quietSpikeFrames = quietSpikeFrames;
stream.trialVoltageAverage = stream.trialVoltageSum ./ ...
    max(stream.stimulatedChunkCount, 1);
stream.trialSpikeAverage = stream.trialSpikeSum ./ ...
    max(stream.stimulatedChunkCount, 1);
stream.blueTrialWaveform = representativeBlueTrial(session);
end

% =========================================================================
function centers = binCenters(edges)
centers = mean([edges(1:end-1); edges(2:end)], 1);
end

function guard = chunkBoundaryGuard(nFrames, framesPerChunk, halfWidth)
% Source pattern (1:nFrame:nFramesTotal+1) + (-halfWidth:halfWidth)', clipped
% into range. Returns linear frame indices to exclude.
boundaries = (1:framesPerChunk:nFrames + 1)';
guard = boundaries + (-halfWidth:halfWidth);
guard = max(min(guard(:), nFrames), 1);
end

function dvdt = centredDifference8(voltage)
% Source line 737: [zeros(n,4) v(:,9:end)-v(:,1:end-8) zeros(n,4)].
[nRows, nCols] = size(voltage);
dvdt = zeros(nRows, nCols);
if nCols > 8
    dvdt(:, 5:nCols - 4) = voltage(:, 9:end) - voltage(:, 1:end - 8);
end
end

function mask = blockSpikeMask(frameIndices, firstFrame, nBlockFrames)
mask = false(1, nBlockFrames);
if isempty(frameIndices)
    return
end
local = frameIndices - firstFrame + 1;
local = local(local >= 1 & local <= nBlockFrames);
mask(local) = true;
end

function quietMask = medianQuietMask(session)
% Source lines 98-104 and 758. The face-camera whisking trace is smoothed with
% an 8-sample moving median, interpolated onto the voltage timebase chunk by
% chunk with bilinear imresize, then thresholded at its own median and dilated
% by a 100-frame moving mean.
nChunks = session.numberOfChunks;
faceFramesPerChunk = session.numberOfFaceFrames / nChunks;
assert(faceFramesPerChunk == round(faceFramesPerChunk), ...
    'Figure5:FaceChunkStructure', ...
    'The face-video frame count is not a whole multiple of the chunk count.');

smoothed = movmedian(session.whiskingMotion, 8);
byChunk = reshape(smoothed, faceFramesPerChunk, nChunks);
interpolated = imresize(byChunk, [session.framesPerChunk nChunks], 'bilinear');
interpolated = reshape(interpolated, 1, []);
assert(numel(interpolated) == session.numberOfVoltageFrames, ...
    'Figure5:WhiskInterpolation', ...
    'Interpolated whisking trace does not match the voltage frame count.');
quietMask = ~(movmean(interpolated > median(interpolated), 100) > 0);
end

function levels = blueStepLevels(session)
% Source line 296: the plateau levels of the first QC cell's blue command
% (blue_all(1,:) after the idxUse subset), that is the values where blue is
% nonzero and not changing. The source scans the whole recording; the first 20
% chunks already contain every level (all four, checked against the whole
% recording on 2026-09-26), so only those are read.
firstCell = session.figureCellIndices(1);
nProbeChunks = min(20, session.numberOfChunks);
lastFrame = nProbeChunks * session.framesPerChunk;
blue = double(session.readBlueRows(firstCell:firstCell, [1 lastFrame]));
held = [0 diff(blue)] == 0;
levels = unique(blue(blue > 0 & held));
assert(~isempty(levels), 'Figure5:NoBlueSteps', ...
    ['Could not find held blue levels for cell %d in the first %d chunks; ', ...
     'the source bVStep cannot be reproduced.'], firstCell, nProbeChunks);
levels = sort(levels(:))';
end

function waveform = representativeBlueTrial(session)
% The blue command of the first chunk in which the first QC cell was
% stimulated, used as the panel A and B stimulation trace (source
% squeeze(blue_all_stim(1,:,1)), after the idxUse subset).
framesPerChunk = session.framesPerChunk;
firstCell = session.figureCellIndices(1);
for iChunk = 1:min(40, session.numberOfChunks)
    range = [(iChunk - 1) * framesPerChunk + 1, iChunk * framesPerChunk];
    blue = double(session.readBlueRows(firstCell:firstCell, range));
    if any(blue > 0)
        waveform = blue;
        return
    end
end
waveform = zeros(1, framesPerChunk);
end

% =========================================================================
function [bin, edges, cellIndices] = quintileBins(values, nBin)
% Source pattern: percentile-spanned equal-width edges, not equal-count bins.
% Cells outside the 1st-99th percentile fall in bin 0 and are excluded.
edges = linspace(prctile(values, 1), prctile(values, 99), nBin + 1);
[~, ~, bin] = histcounts(values, edges);
cellIndices = find(bin > 0);
end

% =========================================================================
function panel = plotExcitabilityByRheobase(session, stream, nBin, outputDirectory)
% Figure 5A. Source fig5_ephys_vs_spiking.m lines 367-438, idxPlt = 2.
propertyColumn = 2;
qc = session.figureCellIndices;
property = session.intrinsicProperties(qc, propertyColumn);
[bin, edges] = quintileBins(property, nBin);

dt = session.framePeriodSeconds;
spontaneous = zeros(nBin, numel(stream.voltageCenters));
ramp = zeros(nBin, numel(stream.voltageCenters));
trialVoltage = zeros(nBin, session.framesPerChunk);
propertyMean = zeros(nBin, 1);
binCount = zeros(nBin, 1);

for iBin = 1:nBin
    rows = qc(bin == iBin);
    binCount(iBin) = numel(rows);
    % Counts are pooled before dividing; this is not a mean of ratios.
    spontaneous(iBin, :) = sum(stream.voltageSpikeHistogram(rows, :), 1) ./ ...
        sum(stream.voltageHistogram(rows, :), 1);
    ramp(iBin, :) = sum(stream.rampSpikeHistogram(rows, :), 1) ./ ...
        sum(stream.rampHistogram(rows, :), 1);
    trialVoltage(iBin, :) = mean(stream.trialVoltageAverage(rows, :), 1);
    propertyMean(iBin) = mean(session.intrinsicProperties(rows, propertyColumn));
end

panel = struct();
panel.panel = 'A';
panel.quintileVariable = session.intrinsicPropertyLabels{propertyColumn};
panel.quintileEdges = edges;
panel.quintileCellCount = binCount;
panel.quintilePropertyMean = propertyMean;
panel.cellCount = sum(bin > 0);
panel.voltageCenters = stream.voltageCenters;
panel.spontaneousRateHz = spontaneous / dt;
panel.rampRateHz = ramp / dt;
panel.trialVoltageAverage = trialVoltage;
panel.trialTimeSeconds = (1:session.framesPerChunk) * dt;

figureHandle = figure('Name', 'Figure 5A', 'Color', 'w', 'Position', [100 100 640 760]);
colours = quintileColours(nBin);

ax1 = subplot(3, 1, 1);
plot(ax1, panel.trialTimeSeconds, stream.blueTrialWaveform, 'b');
ylabel(ax1, 'Blue (a.u.)');
title(ax1, 'Figure 5A stimulation');
xlim(ax1, [4.2 10]);

ax2 = subplot(3, 1, 2);
plotQuintiles(ax2, panel.voltageCenters, panel.rampRateHz, colours);
xlabel(ax2, 'Voltage (norm)');
ylabel(ax2, 'Rate (Hz)');
title(ax2, 'Evoked: ramp stimulation');
xlim(ax2, [0 1]);

ax3 = subplot(3, 1, 3);
plotQuintiles(ax3, panel.voltageCenters, panel.spontaneousRateHz, colours);
xlabel(ax3, 'Voltage (norm)');
ylabel(ax3, 'Rate (Hz)');
title(ax3, sprintf('Spontaneous, quiet state (n = %d cells)', panel.cellCount));
xlim(ax3, [0 1]);

saveFigure(figureHandle, outputDirectory, 'figure5A_excitability_by_rheobase_raw.png');
end

% =========================================================================
function panel = plotSpikeProbabilityByAdaptation(session, stream, nBin, outputDirectory)
% Figure 5B. Source fig5_ephys_vs_spiking.m lines 735-880, idxPlt = 4.
propertyColumn = 4;
qc = session.figureCellIndices;
property = session.intrinsicProperties(qc, propertyColumn);
[bin, edges] = quintileBins(property, nBin);

dt = session.framePeriodSeconds;
voltageWindow = [0.2 0.8];        % source vLim
inWindow = stream.jointVoltageCenters > voltageWindow(1) & ...
    stream.jointVoltageCenters < voltageWindow(2);

marginal = zeros(nBin, numel(stream.dvdtCenters));
trialRate = zeros(nBin, session.framesPerChunk);
propertyMean = zeros(nBin, 1);
binCount = zeros(nBin, 1);

for iBin = 1:nBin
    rows = qc(bin == iBin);
    binCount(iBin) = numel(rows);
    spikeCounts = sum(sum(stream.jointSpikeHistogram(rows, inWindow, :), 1), 2);
    allCounts = sum(sum(stream.jointHistogram(rows, inWindow, :), 1), 2);
    marginal(iBin, :) = squeeze(spikeCounts) ./ squeeze(allCounts);
    trialRate(iBin, :) = mean(stream.trialSpikeAverage(rows, :), 1);
    propertyMean(iBin) = mean(session.intrinsicProperties(rows, propertyColumn));
end

panel = struct();
panel.panel = 'B';
panel.quintileVariable = session.intrinsicPropertyLabels{propertyColumn};
panel.quintileEdges = edges;
panel.quintileCellCount = binCount;
panel.quintilePropertyMean = propertyMean;
panel.cellCount = sum(bin > 0);
panel.dvdtCenters = stream.dvdtCenters;
panel.spikeRateHz = marginal / dt;
panel.voltageWindow = voltageWindow;
panel.trialRateHz = trialRate / dt;
panel.trialTimeSeconds = (1:session.framesPerChunk) * dt;
panel.quietMaskComparison = struct( ...
    'medianRuleFrames', sum(stream.medianQuietMask), ...
    'storedMixtureFrames', sum(stream.storedQuietMask), ...
    'agreementFraction', mean(stream.medianQuietMask == stream.storedQuietMask), ...
    'note', ['The plotted curve uses the source median rule. The stored ', ...
             'quiet_mask is the Gaussian-mixture rule used by panels A and C.']);

figureHandle = figure('Name', 'Figure 5B', 'Color', 'w', 'Position', [100 100 640 760]);
colours = quintileColours(nBin);

ax1 = subplot(3, 1, 1);
plot(ax1, panel.trialTimeSeconds, stream.blueTrialWaveform, 'b');
ylabel(ax1, 'Blue (a.u.)');
title(ax1, 'Figure 5B stimulation');
xlim(ax1, [2.25 3.25]);

ax2 = subplot(3, 1, 2);
plotQuintiles(ax2, panel.trialTimeSeconds, movmean(panel.trialRateHz, 40, 2), colours);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Rate (Hz)');
title(ax2, 'Evoked: step stimulation');
xlim(ax2, [2.25 3.25]);

ax3 = subplot(3, 1, 3);
plotQuintiles(ax3, panel.dvdtCenters, panel.spikeRateHz, colours);
xlabel(ax3, 'dV/dt');
ylabel(ax3, 'Rate (Hz)');
title(ax3, sprintf('Spontaneous, quiet state (n = %d cells)', panel.cellCount));
xlim(ax3, [-0.5 0.5]);

saveFigure(figureHandle, outputDirectory, 'figure5B_spike_probability_by_adaptation_raw.png');
end

% =========================================================================
function panel = plotAutocorrelogramByADP(session, stream, nBin, outputDirectory)
% Figure 5C. Source fig5_ephys_vs_spiking.m lines 541-602, idxPlt = 1, with
% selfStaRate_q from lines 472-473. The top trace is the stored blue-evoked
% STA. The autocorrelogram triggers come from the streaming pass, because the
% source spk_t_quiet (lines 276-287) requires tMask & blueStimMask_all(ii,:),
% the same per-chunk blue-neighbour mask as panel A's spontaneous histogram.
propertyColumn = 1;
qc = session.figureCellIndices;
property = session.intrinsicProperties(qc, propertyColumn);
[bin, edges] = quintileBins(property, nBin);

dt = session.framePeriodSeconds;
nBack = 140;
nFront = 140;
lagFrames = -nBack:nFront;

% Quiet-state autocorrelogram. The source triggers on quiet-state spikes in
% blue-free chunks and counts all spontaneous spikes of the same cell, with
% the zero-lag bin forced to zero (source line 473).
quietSpikes = stream.quietSpikeFrames;
autocorrelogram = zeros(session.numberOfCells, numel(lagFrames));
for iCell = 1:session.numberOfCells
    autocorrelogram(iCell, :) = selfCorrelogram(quietSpikes{iCell}, ...
        session.spikeFrameIndexNoBlue{iCell}, nBack, nFront, ...
        session.numberOfVoltageFrames);
end
autocorrelogram(:, nBack + 1) = 0;

spikeWaveform = session.blueEvokedSpikeSTA;
waveformLag = session.staLagFrames;

waveformByBin = zeros(nBin, numel(waveformLag));
correlogramByBin = zeros(nBin, numel(lagFrames));
propertyMean = zeros(nBin, 1);
binCount = zeros(nBin, 1);
for iBin = 1:nBin
    rows = qc(bin == iBin);
    binCount(iBin) = numel(rows);
    waveformByBin(iBin, :) = mean(spikeWaveform(rows, :), 1);
    correlogramByBin(iBin, :) = mean(autocorrelogram(rows, :), 1);
    propertyMean(iBin) = mean(session.intrinsicProperties(rows, propertyColumn));
end
% Source line 576 baseline-subtracts the single sample at nBack+1-2.
baselineIndex = session.staZeroLagIndex - 2;
waveformByBin = waveformByBin - waveformByBin(:, baselineIndex);

panel = struct();
panel.panel = 'C';
panel.quintileVariable = session.intrinsicPropertyLabels{propertyColumn};
panel.quintileEdges = edges;
panel.quintileCellCount = binCount;
panel.quintilePropertyMean = propertyMean;
panel.cellCount = sum(bin > 0);
panel.waveformLagSeconds = waveformLag * dt;
panel.spikeWaveform = waveformByBin;
panel.correlogramLagSeconds = lagFrames * dt;
panel.autocorrelogramHz = correlogramByBin / dt;
panel.triggerSpikeCount = cellfun(@numel, quietSpikes);

figureHandle = figure('Name', 'Figure 5C', 'Color', 'w', 'Position', [100 100 560 700]);
colours = quintileColours(nBin);

ax1 = subplot(2, 1, 1);
plotQuintiles(ax1, panel.waveformLagSeconds, panel.spikeWaveform, colours);
xlabel(ax1, 'Peri-spike time (s)');
ylabel(ax1, 'Voltage (norm)');
title(ax1, 'Blue-evoked spike waveform');
xlim(ax1, [-1 1] * 20e-3);

ax2 = subplot(2, 1, 2);
plotQuintiles(ax2, panel.correlogramLagSeconds, panel.autocorrelogramHz, colours);
xlabel(ax2, 'Peri-spike time (s)');
ylabel(ax2, 'Rate (Hz)');
title(ax2, sprintf('Spontaneous autocorrelogram (n = %d cells)', panel.cellCount));
xlim(ax2, [-1 1] * 20e-3);

saveFigure(figureHandle, outputDirectory, 'figure5C_autocorrelogram_by_adp_raw.png');
end

function correlogram = selfCorrelogram(triggerFrames, targetFrames, ...
    nBack, nFront, nFrames)
% Probability that the cell spikes at each lag around one of its own trigger
% spikes, matching get_sta_mat_self on a binary spike train.
correlogram = zeros(1, nBack + nFront + 1);
if isempty(triggerFrames) || isempty(targetFrames)
    return
end
keep = triggerFrames > nBack & triggerFrames <= nFrames - nFront;
triggerFrames = triggerFrames(keep);
if isempty(triggerFrames)
    return
end
targetMask = false(1, nFrames);
targetMask(targetFrames) = true;
windows = triggerFrames(:) + (-nBack:nFront);
correlogram = mean(targetMask(windows), 1);
end

% =========================================================================
function quintile = quintileConnectivityAverages(session, mode, nBin)
% Panels D-G share one structure: bin connected pairs into quintiles of a
% connectivity amplitude, average the evoked and spontaneous pairwise series
% with spike-count weights, build a distance-matched unconnected control set per
% quintile, and regress the control out of the spontaneous averages.
%
% mode 'synaptic'     panels D and E, IPSP amplitude quintiles
% mode 'gap_junction' panels F and G, spikelet amplitude quintiles

qc = session.figureCellIndices;
nQC = numel(qc);
dx = session.pixelSizeMicrometers;
distance = session.metrics.distancePixels(qc, qc);
distanceMask = distance > 60 / dx & distance < 400 / dx;
neighbourRadius = 400 / dx;

switch mode
    case 'synaptic'
        connected = session.metrics.synapticConnection(qc, qc);
        amplitude = session.metrics.ipspAmplitude(qc, qc);
        % Source lines 332, 337 and 357-359 [distfix 2026-09-24] restrict the
        % percentiles, the mask and the reported n to the 60-400 um band, the
        % same form as the gap-junction case. Two in-band pairs have a
        % non-finite amplitude, fall in bin 0 and are not counted.
        selectionMask = connected & distanceMask;
        edges = prctile(amplitude(selectionMask), linspace(0, 100, nBin + 1));
        evokedSeries = 'blue_pulse_ipsp';
    case 'gap_junction'
        connected = session.metrics.gapJunctionConnection(qc, qc);
        amplitude = session.metrics.spikeletAmplitude(qc, qc);
        % Source lines 528-529 restrict both the percentiles and the mask to
        % the 60-400 um band.
        selectionMask = connected & distanceMask;
        edges = prctile(amplitude(selectionMask), linspace(0, 100, nBin + 1));
        evokedSeries = 'hadamard_cross_spike_sta';
    otherwise
        error('Figure5:UnknownQuintileMode', 'Unknown mode: %s', mode);
end

[~, ~, bin] = histcounts(amplitude, edges);
quintile = struct();
quintile.mode = mode;
quintile.amplitudeEdges = edges;
quintile.pairCount = sum(bin(:) > 0 & selectionMask(:));
quintile.nBin = nBin;

fprintf('  %s quintiles: n = %d pairs\n', mode, quintile.pairCount);

% ------------------------------------------------------------------ inputs
nBack = 140;
nFront = 140;
dt = session.framePeriodSeconds;
quintile.lagSeconds = (-nBack:nFront) * dt;
fitIndex = 1:(nBack - 1);          % source idxFit = 1:(nBack+1-2)
quintile.fitIndexNote = 'Regression fitted on pre-spike samples 1:nBack-1 only.';

fprintf('    reading pairwise series ...\n');
evoked = session.readPairwiseSeries(evokedSeries, 1:session.numberOfCells, ...
    1:session.numberOfCells, []);
evoked = double(evoked(qc, qc, :));
if strcmp(mode, 'gap_junction')
    % Panel F's top trace, the source mean(selfHadStaFN,1) at line 660: a plain
    % mean over the QC cells of each cell's own spike-triggered voltage. The
    % stored self-STAs sit on the diagonal of hadamard_cross_spike_sta; see
    % documented deviation 3 for how they differ from selfHadStaFN.
    selfSta = zeros(nQC, size(evoked, 3));
    for iCell = 1:nQC
        selfSta(iCell, :) = squeeze(evoked(iCell, iCell, :))';
    end
    quintile.selfSpikeWaveform = mean(selfSta, 1)';
end
if strcmp(mode, 'synaptic')
    crosstalk = session.readPairwiseSeries('blue_light_crosstalk', ...
        1:session.numberOfCells, 1:session.numberOfCells, []);
    evoked = evoked - double(crosstalk(qc, qc, :));
end
spontaneousVoltage = double(session.readPairwiseSeries( ...
    'spontaneous_cross_spike_sta', 1:session.numberOfCells, ...
    1:session.numberOfCells, []));
spontaneousVoltage = spontaneousVoltage(qc, qc, :);
spontaneousRate = double(session.readPairwiseSeries( ...
    'spontaneous_spike_cross_correlogram', 1:session.numberOfCells, ...
    1:session.numberOfCells, []));
spontaneousRate = spontaneousRate(qc, qc, :);

spontaneousWeight = repmat(session.noBlueSpikeCount(qc), 1, nQC);
switch mode
    case 'synaptic'
        evokedWeight = [];      % panel D uses a plain nanmean
    case 'gap_junction'
        evokedWeight = session.hadamardSpikeCount(qc, qc);
end

nLag = size(spontaneousVoltage, 3);
quintile.evokedAverage = zeros(nLag, nBin);
quintile.evokedControl = zeros(nLag, nBin);
quintile.spontaneousVoltage = zeros(nLag, nBin);
quintile.spontaneousVoltageControl = zeros(nLag, nBin);
quintile.spontaneousVoltageResidual = zeros(nLag, nBin);
quintile.spontaneousRate = zeros(nLag, nBin);
quintile.spontaneousRateControl = zeros(nLag, nBin);
quintile.spontaneousRateResidual = zeros(nLag, nBin);
quintile.pairsPerBin = zeros(nBin, 1);
quintile.controlPairsPerBin = zeros(nBin, 1);
quintile.voltageRegressionCoefficient = zeros(nBin, 1);
quintile.rateRegressionCoefficient = zeros(nBin, 1);

for iBin = 1:nBin
    mask = bin == iBin & selectionMask;
    control = controlPairMask(mask, connected, distanceMask, distance, ...
        neighbourRadius, nQC);
    quintile.pairsPerBin(iBin) = sum(mask(:));
    quintile.controlPairsPerBin(iBin) = sum(control(:));

    if isempty(evokedWeight)
        quintile.evokedAverage(:, iBin) = maskedMean(evoked, mask);
        quintile.evokedControl(:, iBin) = maskedMean(evoked, control);
    else
        quintile.evokedAverage(:, iBin) = weightedMean(evoked, mask, evokedWeight);
        quintile.evokedControl(:, iBin) = weightedMean(evoked, control, evokedWeight);
    end

    voltageAverage = weightedMean(spontaneousVoltage, mask, spontaneousWeight);
    voltageControl = weightedMean(spontaneousVoltage, control, spontaneousWeight);
    rateAverage = weightedMean(spontaneousRate, mask, spontaneousWeight);
    rateControl = weightedMean(spontaneousRate, control, spontaneousWeight);

    quintile.spontaneousVoltage(:, iBin) = voltageAverage;
    quintile.spontaneousVoltageControl(:, iBin) = voltageControl;
    quintile.spontaneousRate(:, iBin) = rateAverage;
    quintile.spontaneousRateControl(:, iBin) = rateControl;

    [quintile.spontaneousVoltageResidual(:, iBin), ...
        quintile.voltageRegressionCoefficient(iBin)] = ...
        removeSharedDynamics(voltageAverage, voltageControl, fitIndex);
    [quintile.spontaneousRateResidual(:, iBin), ...
        quintile.rateRegressionCoefficient(iBin)] = ...
        removeSharedDynamics(rateAverage, rateControl, fitIndex);
end

quintile.zeroLagIndex = nBack + 1;
quintile.framePeriodSeconds = dt;
end

function control = controlPairMask(mask, connected, distanceMask, distance, ...
    neighbourRadius, nQC)
% Source lines 363-370 and 560-567. For every presynaptic cell that contributes
% a pair to this quintile, the controls are cells within neighbourRadius of any
% of that cell's quintile partners, that are not connected to it, and that lie
% in the distance band.
control = false(nQC);
preCells = find(any(mask, 2));
for iPre = 1:numel(preCells)
    row = preCells(iPre);
    nearPartner = any(distance(mask(row, :), :) < neighbourRadius, 1);
    control(row, :) = nearPartner & ~connected(row, :) & distanceMask(row, :);
end
end

function average = maskedMean(series, mask)
% Plain nanmean over the masked pairs, matching source line 374.
maskNan = double(mask);
maskNan(~mask) = NaN;
masked = series .* maskNan;
average = squeeze(nanmean(nanmean(masked, 1), 2));
average = average(:);
end

function average = weightedMean(series, mask, weight)
% Spike-count-weighted average, matching the source nansum(X.*maskNan.*n)/
% nansum(n.*maskNan) pattern.
maskNan = double(mask);
maskNan(~mask) = NaN;
weighted = maskNan .* weight;
numerator = squeeze(nansum(nansum(series .* weighted, 1), 2));
denominator = nansum(nansum(weighted, 1), 2);
average = numerator(:) / denominator;
end

function [residual, coefficient] = removeSharedDynamics(trace, reference, fitIndex)
% Source lines 392-412 and 596-618. Both traces are mean-subtracted over the
% full window, the scale factor is least-squares fitted on the pre-spike
% samples only, and the residual is formed over the whole window.
% SeeResiduals_vec solves the same single-regressor least-squares problem,
% reproduced here so the reproduction does not depend on a lab helper.
trace = trace(:) - mean(trace(:));
reference = reference(:) - mean(reference(:));
x = reference(fitIndex);
y = trace(fitIndex);
if all(x == 0) || any(~isfinite(x)) || any(~isfinite(y))
    coefficient = 0;
else
    coefficient = (x' * y) / (x' * x);
end
residual = trace - coefficient * reference;
end

% =========================================================================
function panel = plotIpspWaveforms(session, quintile, outputDirectory)
% Figure 5D. Source lines 374 (average) and 452-465 (plot): blue pulse
% waveform above the quintile average IPSP response. No baseline subtraction;
% the crosstalk estimate was already removed inside the average. The current
% caption no longer lists a blue waveform on top of panel D; the blue-pulse
% subplot is kept here as a diagnostic, as the source still draws it.
panel = struct();
panel.panel = 'D';
panel.lagSeconds = quintile.lagSeconds;
panel.ipspWaveform = quintile.evokedAverage;
panel.pairCount = quintile.pairCount;
panel.pairsPerBin = quintile.pairsPerBin;
panel.amplitudeEdges = quintile.amplitudeEdges;
panel.bluePulseWaveform = bluePulseTriggeredBlue(session, 140, 140);

figureHandle = figure('Name', 'Figure 5D', 'Color', 'w', 'Position', [100 100 560 620]);
colours = quintileColours(quintile.nBin);

ax1 = subplot(2, 1, 1);
plot(ax1, panel.lagSeconds, panel.bluePulseWaveform, 'b');
ylabel(ax1, 'Blue (a.u.)');
title(ax1, 'Figure 5D blue pulse');

ax2 = subplot(2, 1, 2);
plotQuintiles(ax2, panel.lagSeconds, panel.ipspWaveform', colours);
xlabel(ax2, 'Peri-spike time (s)');
ylabel(ax2, 'Voltage (SH)');
title(ax2, sprintf('IPSP by amplitude quintile (n = %d pairs)', panel.pairCount));

linkaxes([ax1 ax2], 'x');
xlim(ax1, [panel.lagSeconds(1) panel.lagSeconds(end)]);

saveFigure(figureHandle, outputDirectory, 'figure5D_ipsp_waveforms_raw.png');
end

function waveform = bluePulseTriggeredBlue(session, nBack, nFront)
% Source bluePulseStaBlue: the blue command averaged around its own sparse-pulse
% rising edges. The sparse-pulse protocol occupies the first three of every four
% recording chunks, so one early block is enough to recover the pulse shape.
framesPerChunk = session.framesPerChunk;
nProbeChunks = min(3, session.numberOfChunks);
blue = double(session.readBlueRows(1:session.numberOfCells, ...
    [1 nProbeChunks * framesPerChunk]));
rising = find(any(diff(blue > 0, 1, 2) > 0, 1)) + 1;
rising = rising(rising > nBack & rising <= size(blue, 2) - nFront);
if isempty(rising)
    waveform = zeros(1, nBack + nFront + 1);
    return
end
pooled = any(blue > 0, 1);
windows = rising(:) + (-nBack:nFront);
waveform = mean(double(pooled(windows)), 1);
end

% =========================================================================
function panel = plotSpikeletWaveforms(quintile, outputDirectory)
% Figure 5F. Source lines 670-671: the control trace is subtracted from the
% quintile average and the result is baseline-subtracted at the single sample
% nBack+1-2. The top trace is the spike waveform, source line 660, plotted
% without a baseline as in the source.
baselineIndex = quintile.zeroLagIndex - 2;
waveform = quintile.evokedAverage - quintile.evokedControl;
waveform = waveform - waveform(baselineIndex, :);

panel = struct();
panel.panel = 'F';
panel.lagSeconds = quintile.lagSeconds;
panel.spikeletWaveform = waveform;
% Mean over QC cells of the stored self-STA diagonal (deviation 3). Before
% 2026-09-26 this field wrongly averaged the five quintile spikelet traces.
panel.preJunctionWaveform = quintile.selfSpikeWaveform;
panel.pairCount = quintile.pairCount;
panel.pairsPerBin = quintile.pairsPerBin;
panel.amplitudeEdges = quintile.amplitudeEdges;
panel.baselineNote = ['Control-subtracted, then the single sample at ', ...
    'nBack+1-2 removed. Source lines 670-671.'];

figureHandle = figure('Name', 'Figure 5F', 'Color', 'w', 'Position', [100 100 560 620]);
colours = quintileColours(quintile.nBin);

ax1 = subplot(2, 1, 1);
plot(ax1, panel.lagSeconds, panel.preJunctionWaveform, 'k');
ylabel(ax1, 'Voltage (SH)');
title(ax1, 'Spike waveform (mean of stored self-STA diagonal)');
xlim(ax1, [-1 1] * 10e-3);

ax2 = subplot(2, 1, 2);
plotQuintiles(ax2, panel.lagSeconds, waveform', colours);
xlabel(ax2, 'Peri-spike time (s)');
ylabel(ax2, 'Voltage (SH)');
title(ax2, sprintf('Spikelet by amplitude quintile (n = %d pairs)', panel.pairCount));
xlim(ax2, [-1 1] * 10e-3);

saveFigure(figureHandle, outputDirectory, 'figure5F_spikelet_waveforms_raw.png');
end

% =========================================================================
function panel = plotSpontaneousResiduals(quintile, panelLetter, quintileName, ...
    outputDirectory)
% Figures 5E and 5G. Both plot the cross STA voltage above the cross STA
% firing rate, each measured relative to distance-matched unconnected controls
% (the regression in removeSharedDynamics). Their baseline conventions differ,
% and panel G's rate trace deliberately has none.
dt = quintile.framePeriodSeconds;
zeroLag = quintile.zeroLagIndex;

voltage = quintile.spontaneousVoltageResidual;
rate = quintile.spontaneousRateResidual / dt;
switch panelLetter
    case 'E'
        voltage = voltage - mean(voltage(zeroLag - 2, :), 1);
        rate = rate - mean(rate(zeroLag + (-32:-2), :), 1);
        baselineNote = ['Voltage: single sample nBack+1-2. Rate: mean of ', ...
            'nBack+1+(-32:-2). Source lines 468 and 479.'];
        xLimit = [quintile.lagSeconds(1) quintile.lagSeconds(end)];
    case 'G'
        voltage = voltage - mean(voltage(zeroLag + (-10:-2), :), 1);
        baselineNote = ['Voltage: mean of nBack+1+(-10:-2). Rate: none, the ', ...
            'source baseline line 690 is commented out. Voltage baseline ', ...
            'source line 679.'];
        xLimit = [-1 1] * 10e-3;
    otherwise
        error('Figure5:UnknownResidualPanel', 'Unknown panel: %s', panelLetter);
end

panel = struct();
panel.panel = panelLetter;
panel.quintileVariable = quintileName;
panel.lagSeconds = quintile.lagSeconds;
panel.crossStaVoltage = voltage;
panel.crossStaRateHz = rate;
panel.pairCount = quintile.pairCount;
panel.pairsPerBin = quintile.pairsPerBin;
panel.controlPairsPerBin = quintile.controlPairsPerBin;
panel.voltageRegressionCoefficient = quintile.voltageRegressionCoefficient;
panel.rateRegressionCoefficient = quintile.rateRegressionCoefficient;
panel.baselineNote = baselineNote;

figureHandle = figure('Name', ['Figure 5' panelLetter], 'Color', 'w', 'Position', [100 100 560 620]);
colours = quintileColours(quintile.nBin);

ax1 = subplot(2, 1, 1);
plotQuintiles(ax1, panel.lagSeconds, voltage', colours);
ylabel(ax1, 'Voltage (SH)');
title(ax1, sprintf('Figure 5%s cross STA voltage (n = %d pairs)', ...
    panelLetter, panel.pairCount));
xlim(ax1, xLimit);

ax2 = subplot(2, 1, 2);
plotQuintiles(ax2, panel.lagSeconds, rate', colours);
xlabel(ax2, 'Peri-spike time (s)');
ylabel(ax2, 'P(post spk | pre spk) (Hz)');
title(ax2, 'Cross STA rate, relative to distance-matched controls');
xlim(ax2, xLimit);

saveFigure(figureHandle, outputDirectory, ...
    sprintf('figure5%s_cross_sta_%s_raw.png', panelLetter, quintileName));
end

% =========================================================================
function colours = quintileColours(nBin)
% The source uses colorcet('cbd2') for panels A-C and winter/autumn ramps for
% D-G. Those are cosmetic, and colorcet is a third-party helper, so this
% reproduction uses a built-in perceptual ramp of the same length.
colours = parula(nBin);
end

function plotQuintiles(ax, x, y, colours)
% y is [nSample x nBin] or [nBin x nSample]; orient to match x.
if size(y, 1) ~= numel(x) && size(y, 2) == numel(x)
    y = y';
end
hold(ax, 'on');
for iBin = 1:size(y, 2)
    plot(ax, x, y(:, iBin), 'LineWidth', 1.5, 'Color', colours(iBin, :));
end
hold(ax, 'off');
box(ax, 'on');
set(ax, 'FontSize', 10);
end

function saveFigure(figureHandle, outputDirectory, fileName)
% Matches reproduce_figure3/4_from_nwb.m: print at -r300, then leave the figure
% open and visible so the panels can be inspected interactively after the run.
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
