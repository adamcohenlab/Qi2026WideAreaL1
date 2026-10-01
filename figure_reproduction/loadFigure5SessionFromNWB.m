function session = loadFigure5SessionFromNWB(nwbFile)
% LOADFIGURE5SESSIONFROMNWB Load Figure 5 inputs from either session type.
%
%   SESSION = LOADFIGURE5SESSIONFROMNWB(NWBFILE) reads one of the two NWB files
%   Figure 5 needs and returns ordinary MATLAB arrays plus lazy readers for the
%   large datasets. The session type is detected from the file contents, so the
%   same call works for both:
%
%     'sparseOpto'      panels A, B, C. Intrinsic properties measured under
%                       optogenetic stimulation, compared with spontaneous
%                       activity in the same cells.
%     'sparsePulseHad'  panels D, E, F, G. Synaptic and gap-junctional
%                       connectivity, compared with spontaneous activity.
%
%   Figure 5 uses **M-YQ0201-27 for every panel** (confirmed by the author on
%   2026-09-18, and independently reproduced: the caption's n = 280 cells,
%   n = 800 pairs and n = 3553 pairs all come out of that mouse, while
%   M-YQ0201-29 gives 251, 667 and 2199). The D/E count of 800 is on the
%   sparsePulseHad deposit regenerated on 2026-09-25, with the 60-400 um band
%   the panel script applies since its [distfix 2026-09-24] change; it was 814
%   on the earlier deposit without the band. The two-mouse numbers quoted in
%   the Figure 5 narrative come from sup_fig5_connMap_stats_norm.m, not from
%   the panel scripts.
%
%   The loader uses built-in HDF5 functions only, so it is compatible with
%   MATLAB R2019b and does not require MatNWB at reproduction time. Dataset
%   paths follow the embedded MATLAB-to-NWB lookup tables written by
%   convert_sparseOpto_to_nwb.m and convert_sparsePulseHad_to_nwb.m.
%
%   Relationship to the other loaders
%   ---------------------------------
%   For a sparsePulseHad file this function **delegates to
%   LOADFIGURE3SESSIONFROMNWB** for the units axis, the pairwise connectivity
%   metrics and the pairwise timebase, then adds the fields Figure 5 needs on
%   top. Figure 3 and Figure 5 read the same connectivity products, so there is
%   one implementation of that layer, not two. Everything the Figure 3 loader
%   returns is still present on the returned struct.
%
%   Cell axis
%   ---------
%   Both files store all 320 detected cells and mark the figure scripts'
%   selection in /units/passed_qc. Unlike the Figure 4 behavior files there is
%   no second, pre-subset cell axis: every array here is on the 320-cell units
%   axis, and the Figure 5 panels restrict to session.figureCellIndices.
%
%   The two selection rules are deliberately different and must not be unified:
%     sparseOpto      spkHgtNoBlue > 0 & fpRates < 0.5 & idxUseMerge
%                     -> 286 of 320 cells
%     sparsePulseHad  the same, plus spkHgtNoBlue > mean - 2*std
%                     -> 290 of 320 cells
%   dataPrepsForMAT.m applies each file's own rule, so passed_qc is exactly the
%   corresponding figure script's selection.
%
%   Units and conventions
%   ---------------------
%   Voltage      : dimensionless, already divided by each cell's mean no-blue
%                  spike height ("SH" in the source scripts). Both
%                  normalized_voltage and subthreshold_voltage are normalized.
%   Coordinates  : cellCoordinatePixels is in camera pixels, as in the source
%                  cell_coord. Multiply by pixelSizeMicrometers (6.5) for
%                  microns. Distances are returned in pixels, the unit the
%                  source distance limits are expressed in.
%   Time         : NWB stores seconds. framePeriodSeconds is 1.27e-3 s, the
%                  voltage-camera exposure, which is the source scripts' dt.
%   Spike times  : /units/spike_times and /units/spike_times_no_blue are in
%                  seconds. This loader converts them back to the 1-based
%                  voltage-frame indices the source scripts use, in
%                  spikeFrameIndex and spikeFrameIndexNoBlue.
%   Masks        : quiet_mask and whisking_mask are stored as uint8 and are
%                  returned as logical.
%
%   sparseOpto fields (panels A, B, C)
%   ----------------------------------
%   intrinsicProperties   [nCells x 5] in the source props_all column order:
%                         [adp_blue vRheobase frMax frPulseAdp membraneC].
%                         Panel A uses column 2, panel B column 4, panel C
%                         column 1. intrinsicPropertyLabels names them.
%   blueEvokedSpikeSTA    [nCells x 281] source selfBlueStaFN_ADP, the panel C
%                         top trace. staLagFrames and staZeroLagIndex give its
%                         lag axis (zero lag at index 141).
%   framesPerChunk        8000 for sparseOpto. The recording is 300 such chunks.
%   readVoltageRows       Lazy reader; see the performance note below.
%   readBlueRows          Lazy reader for the per-cell blue stimulation.
%
%   sparsePulseHad extra fields (panels D, E, F, G)
%   ----------------------------------------------
%   readPairwiseSeries    Extended beyond the Figure 3 loader with the two
%                         spontaneous series Figure 5 panels E and G need:
%                         'spontaneous_cross_spike_sta' and
%                         'spontaneous_spike_cross_correlogram'.
%   noBlueSpikeCount      [nCells x 1] per-cell spontaneous spike counts, the
%                         source nNoBlueSpk weights for panels E and G.
%   hadamardSpikeCount    [nCells x nCells] per-pair Hadamard-orthogonal
%                         presynaptic spike counts, the source nHadSpk weights
%                         for panel F. Recovered from the ragged
%                         hadamard_orthogonal_spike_times column added on
%                         2026-09-18; see HADAMARD_ORTHOGONAL_SPIKES_PATCH.md.
%                         A zero means "not computed" for pairs beyond the
%                         source rLim of 400 um, not "no spikes".
%   readHadamardSpikeTimes(i, j)
%                         The per-pair trigger spike times themselves, in
%                         1-based voltage-frame indices.
%
%   Reading the voltage datasets efficiently
%   ----------------------------------------
%   normalized_voltage, subthreshold_voltage and blue_stimulation_per_cell are
%   [nCells x nFrames] in MATLAB order and HDF5-chunked as [nCells x 10000], so
%   every chunk spans *all* cells for a 10000-frame window. Reading a block of
%   frames for all cells is contiguous and fast; reading one cell's full row
%   costs the same as reading all cells, because it must still touch every
%   chunk. Iterate over frame blocks, not cells. Every Figure 5 quantity that
%   needs the full recording is a histogram or triggered-average accumulation
%   over time, so block-wise accumulation is exact.
%   session.voltageChunkFrames reports the chunk extent along the frame axis.
%
%   Example
%   -------
%     opto = loadFigure5SessionFromNWB('...\M_YQ0201_27_sparseOpto.nwb');
%     rheobase = opto.intrinsicProperties(opto.figureCellIndices, 2);
%     block = opto.readVoltageRows('subthreshold_voltage', 1:opto.numberOfCells, [1 10000]);
%
%     conn = loadFigure5SessionFromNWB('...\M_YQ0201_27_sparsePulseHad.nwb');
%     ipsp = conn.readPairwiseSeries('blue_pulse_ipsp', 1:conn.numberOfCells, ...
%         1:conn.numberOfCells, []);
%
%   See also REPRODUCE_FIGURE5_FROM_NWB, LOADFIGURE3SESSIONFROMNWB,
%   LOADFIGURE4SESSIONFROMNWB.

if nargin < 1 || isempty(nwbFile)
    error('Figure5NWB:MissingFile', ...
        'Provide a sparseOpto or sparsePulseHad NWB file.');
end
validateattributes(nwbFile, {'char'}, {'row', 'nonempty'}, mfilename, 'nwbFile');
assert(exist(nwbFile, 'file') == 2, 'Figure5NWB:FileNotFound', ...
    'NWB file not found: %s', nwbFile);

sessionType = detectSessionType(nwbFile);
switch sessionType
    case 'sparseOpto'
        session = loadSparseOptoSession(nwbFile);
    case 'sparsePulseHad'
        session = loadSparsePulseHadSession(nwbFile);
    otherwise
        error('Figure5NWB:UnknownSessionType', ...
            'Unrecognized session type: %s', sessionType);
end
session.sessionType = sessionType;
session.nwbFile = nwbFile;
end

% =========================================================================
function sessionType = detectSessionType(nwbFile)
% The five intrinsic-property columns exist only in a sparseOpto file, and the
% pairwise connectivity table only in a sparsePulseHad file. Probe for both
% rather than guessing from the file name.
hasProperties = datasetExists(nwbFile, '/units/optical_rheobase_norm');
hasPairwise = datasetExists(nwbFile, ...
    '/analysis/pairwise_connectivity_metrics/synaptic_connection');
assert(hasProperties || hasPairwise, 'Figure5NWB:UnrecognizedFile', ...
    ['%s is neither a sparseOpto file (no /units/optical_rheobase_norm) nor a ', ...
     'sparsePulseHad file (no pairwise connectivity table).'], nwbFile);
assert(~(hasProperties && hasPairwise), 'Figure5NWB:AmbiguousFile', ...
    '%s looks like both session types at once.', nwbFile);
if hasProperties
    sessionType = 'sparseOpto';
else
    sessionType = 'sparsePulseHad';
end
end

% =========================================================================
function session = loadSparseOptoSession(nwbFile)
paths = sparseOptoPaths();
required = struct2cell(paths);
for iPath = 1:numel(required)
    assertDatasetExists(nwbFile, required{iPath});
end

session = struct();
session.paths = paths;

% ------------------------------------------------------------- units axis
session.cellId = double(h5read(nwbFile, paths.cellId));
session.cellId = session.cellId(:);
session.numberOfCells = numel(session.cellId);
assert(session.numberOfCells > 0, 'Figure5NWB:NoCells', 'The units table is empty.');

session.cellCoordinatePixels = [ ...
    readColumn(nwbFile, paths.xPixel, session.numberOfCells, 'x_pixel'), ...
    readColumn(nwbFile, paths.yPixel, session.numberOfCells, 'y_pixel')];
session.passedQC = logical(readColumn(nwbFile, paths.passedQC, ...
    session.numberOfCells, 'passed_qc'));
session.gfpPositive = logical(readColumn(nwbFile, paths.gfpPositive, ...
    session.numberOfCells, 'gfp_positive'));
session.figureCellIndices = find(session.passedQC);
session.pixelSizeMicrometers = pixelSizeFromUnits(nwbFile, paths, ...
    session.cellCoordinatePixels);
session.cellDistancePixels = squareform(pdist(session.cellCoordinatePixels));

% ------------------------------------------- intrinsic properties, props_all
% Column order is the source props_all order in fig5_ephys_vs_spiking.m line 80.
session.intrinsicPropertyLabels = {'ADP (SH)', 'Rheobase (SH)', 'frMax (Hz)', ...
    'Adaptation', 'membraneC (ms)'};
session.intrinsicProperties = [ ...
    readColumn(nwbFile, paths.adpBlue, session.numberOfCells, 'adp_blue_norm'), ...
    readColumn(nwbFile, paths.rheobase, session.numberOfCells, 'optical_rheobase_norm'), ...
    readColumn(nwbFile, paths.maxFiringRate, session.numberOfCells, 'max_firing_rate_hz'), ...
    readColumn(nwbFile, paths.pulseAdaptation, session.numberOfCells, 'pulse_adaptation'), ...
    readColumn(nwbFile, paths.membraneConstant, session.numberOfCells, 'membrane_time_constant_ms')];

% --------------------------------------------------------------- timebase
session.voltageTimestampsSeconds = double(h5read(nwbFile, paths.voltageTimestamps));
session.voltageTimestampsSeconds = session.voltageTimestampsSeconds(:);
session.numberOfVoltageFrames = numel(session.voltageTimestampsSeconds);
session.framePeriodSeconds = 1.27e-3;

voltageInfo = h5info(nwbFile, paths.subthresholdVoltage);
assert(isequal(voltageInfo.Dataspace.Size, ...
    [session.numberOfCells session.numberOfVoltageFrames]), ...
    'Figure5NWB:VoltageLayout', ...
    'Expected MATLAB HDF5 layout cell x frame at %s.', paths.subthresholdVoltage);

session.recordingChunkStartSeconds = double(h5read(nwbFile, paths.chunkStart));
session.recordingChunkStopSeconds = double(h5read(nwbFile, paths.chunkStop));
session.numberOfChunks = numel(session.recordingChunkStartSeconds);
assert(session.numberOfChunks > 0 && ...
    mod(session.numberOfVoltageFrames, session.numberOfChunks) == 0, ...
    'Figure5NWB:ChunkStructure', ...
    'The voltage frame count is not a whole multiple of the chunk count.');
session.framesPerChunk = session.numberOfVoltageFrames / session.numberOfChunks;

% ------------------------------------------------ blue-evoked spike-triggered
% Source selfBlueStaFN_ADP, the panel C top trace, stored as [nCells x nLag].
staInfo = h5info(nwbFile, paths.blueEvokedSTA);
assert(numel(staInfo.Dataspace.Size) == 2 && ...
    staInfo.Dataspace.Size(1) == session.numberOfCells, ...
    'Figure5NWB:BlueStaLayout', ...
    'Expected MATLAB HDF5 layout cell x lag at %s.', paths.blueEvokedSTA);
session.blueEvokedSpikeSTA = double(h5read(nwbFile, paths.blueEvokedSTA));
nStaSamples = size(session.blueEvokedSpikeSTA, 2);
assert(mod(nStaSamples, 2) == 1, 'Figure5NWB:BlueStaLength', ...
    'The blue-evoked STA must have an odd sample count so zero lag is centered.');
session.staZeroLagIndex = (nStaSamples + 1) / 2;
session.staLagFrames = (1:nStaSamples) - session.staZeroLagIndex;
session.staLagSeconds = session.staLagFrames * session.framePeriodSeconds;

% -------------------------------------------------------------- behavior
session.runningSpeed = readRowVector(nwbFile, paths.runningSpeed, ...
    session.numberOfVoltageFrames, 'running_speed');
session.quietMask = logical(readRowVector(nwbFile, paths.quietMask, ...
    session.numberOfVoltageFrames, 'quiet_mask'));
session.whiskingMask = logical(readRowVector(nwbFile, paths.whiskingMask, ...
    session.numberOfVoltageFrames, 'whisking_mask'));
session.whiskingMotion = double(h5read(nwbFile, paths.whiskingMotion));
session.whiskingMotion = session.whiskingMotion(:)';
session.whiskingMotionTimestampsSeconds = ...
    double(h5read(nwbFile, paths.whiskingMotionTimestamps));
session.numberOfFaceFrames = numel(session.whiskingMotion);
session.maskProvenanceNote = [ ...
    'quiet_mask and whisking_mask were produced by an independent, unseeded ', ...
    'fitgmdist draw in dataPrepsForMAT.m, so they are the same definition as, ', ...
    'but not bit-identical to, the masks used when the legacy figure script ran. ', ...
    'Panel B''s source block uses a different, median-threshold rule; see ', ...
    'reproduce_figure5_from_nwb.m.'];

% ----------------------------------------------------------- spike times
% spike_times_no_blue was built by dataPrepsForMAT.m lines 646-654 from the raw
% spike matrix, excluding frames with blue on any cell within 80 um among all
% 320 cells. fig5_ephys_vs_spiking.m pools only the QC cells; see deviation 5
% in reproduce_figure5_from_nwb.m for the measured difference.
session.spikeFrameIndex = readSpikeFrameIndices(nwbFile, paths.spikeTimes, ...
    paths.spikeTimesIndex, session.numberOfCells, ...
    session.voltageTimestampsSeconds, session.framePeriodSeconds, 'spike_times');
session.spikeFrameIndexNoBlue = readSpikeFrameIndices(nwbFile, ...
    paths.spikeTimesNoBlue, paths.spikeTimesNoBlueIndex, session.numberOfCells, ...
    session.voltageTimestampsSeconds, session.framePeriodSeconds, ...
    'spike_times_no_blue');
session.spikeCount = cellfun(@numel, session.spikeFrameIndex);
session.noBlueSpikeCount = cellfun(@numel, session.spikeFrameIndexNoBlue);

% -------------------------------------------------------------- readers
session.readVoltageRows = @(seriesName, cellIndex, frameRange) ...
    readCellFrameDataset(nwbFile, sparseOptoSeriesPath(paths, seriesName), ...
    session.numberOfCells, session.numberOfVoltageFrames, cellIndex, frameRange);
session.readBlueRows = @(cellIndex, frameRange) ...
    readCellFrameDataset(nwbFile, paths.blueStimulation, ...
    session.numberOfCells, session.numberOfVoltageFrames, cellIndex, frameRange);
session.voltageChunkFrames = chunkFrames(nwbFile, paths.subthresholdVoltage);
end

% =========================================================================
function session = loadSparsePulseHadSession(nwbFile)
% Reuse the Figure 3 loader for the units axis, pairwise metrics and timebase,
% then add what Figure 5 panels D-G need on top.
assert(~isempty(which('loadFigure3SessionFromNWB')), 'Figure5NWB:MissingFigure3Loader', ...
    ['loadFigure3SessionFromNWB.m must be on the path; Figure 5 panels D-G ', ...
     'reuse it for the pairwise connectivity layer.']);
session = loadFigure3SessionFromNWB(nwbFile);

paths = sparsePulseHadExtraPaths();
extra = struct2cell(paths);
for iPath = 1:numel(extra)
    assertDatasetExists(nwbFile, extra{iPath});
end
session.extraPaths = paths;

nCells = session.numberOfCells;
session.figureCellIndices = find(session.passedQC);
session.framePeriodSeconds = 1.27e-3;

session.voltageTimestampsSeconds = double(h5read(nwbFile, paths.voltageTimestamps));
session.voltageTimestampsSeconds = session.voltageTimestampsSeconds(:);
session.numberOfVoltageFrames = numel(session.voltageTimestampsSeconds);

session.recordingChunkStartSeconds = double(h5read(nwbFile, paths.chunkStart));
session.recordingChunkStopSeconds = double(h5read(nwbFile, paths.chunkStop));
session.numberOfChunks = numel(session.recordingChunkStartSeconds);
assert(session.numberOfChunks > 0 && ...
    mod(session.numberOfVoltageFrames, session.numberOfChunks) == 0, ...
    'Figure5NWB:ChunkStructure', ...
    'The voltage frame count is not a whole multiple of the chunk count.');
session.framesPerChunk = session.numberOfVoltageFrames / session.numberOfChunks;

% -------------------------------------------- spontaneous spike-time weights
session.spikeFrameIndex = readSpikeFrameIndices(nwbFile, paths.spikeTimes, ...
    paths.spikeTimesIndex, nCells, session.voltageTimestampsSeconds, ...
    session.framePeriodSeconds, 'spike_times');
session.spikeFrameIndexNoBlue = readSpikeFrameIndices(nwbFile, ...
    paths.spikeTimesNoBlue, paths.spikeTimesNoBlueIndex, nCells, ...
    session.voltageTimestampsSeconds, session.framePeriodSeconds, ...
    'spike_times_no_blue');
session.spikeCount = cellfun(@numel, session.spikeFrameIndex);
% Source nNoBlueSpk, the panel E and G averaging weights.
session.noBlueSpikeCount = cellfun(@numel, session.spikeFrameIndexNoBlue);

% ------------------------------- per-pair Hadamard trigger counts, nHadSpk
[session.hadamardSpikeCount, hadamardEndOffset] = ...
    readHadamardSpikeCounts(nwbFile, paths, nCells);
session.hadamardSpikeCountNote = [ ...
    'Source nHadSpk = cellfun(@length, spk_t_had_ortho). Zero means the pair ', ...
    'was beyond the source rLim of 400 um and was never computed, not that it ', ...
    'had no spikes. See HADAMARD_ORTHOGONAL_SPIKES_PATCH.md.'];
session.readHadamardSpikeTimes = @(preIndex, postIndex) ...
    readHadamardPairSpikeFrames(nwbFile, paths, nCells, hadamardEndOffset, ...
    session.voltageTimestampsSeconds, session.framePeriodSeconds, ...
    preIndex, postIndex);

% ------------------------------------------- extended pairwise series reader
% The Figure 3 loader knows three pairwise series. Figure 5 panels E and G also
% need the two spontaneous ones, so wrap its reader and fall through.
figure3Reader = session.readPairwiseSeries;
session.readPairwiseSeries = @(seriesName, preIndex, postIndex, timeRange) ...
    readFigure5PairwiseSeries(nwbFile, paths, nCells, figure3Reader, ...
    seriesName, preIndex, postIndex, timeRange);

% --------------------------------------------------- cell x frame readers
session.readVoltageRows = @(cellIndex, frameRange) ...
    readCellFrameDataset(nwbFile, paths.normalizedVoltage, nCells, ...
    session.numberOfVoltageFrames, cellIndex, frameRange);
session.readBlueRows = @(cellIndex, frameRange) ...
    readCellFrameDataset(nwbFile, paths.blueStimulation, nCells, ...
    session.numberOfVoltageFrames, cellIndex, frameRange);
session.voltageChunkFrames = chunkFrames(nwbFile, paths.normalizedVoltage);
end

% =========================================================================
function paths = sparseOptoPaths()
paths = struct( ...
    'cellId', '/units/id', ...
    'xPixel', '/units/x_pixel', ...
    'yPixel', '/units/y_pixel', ...
    'xMicrometers', '/units/x_um', ...
    'yMicrometers', '/units/y_um', ...
    'passedQC', '/units/passed_qc', ...
    'gfpPositive', '/units/gfp_positive', ...
    'adpBlue', '/units/adp_blue_norm', ...
    'rheobase', '/units/optical_rheobase_norm', ...
    'maxFiringRate', '/units/max_firing_rate_hz', ...
    'pulseAdaptation', '/units/pulse_adaptation', ...
    'membraneConstant', '/units/membrane_time_constant_ms', ...
    'spikeTimes', '/units/spike_times', ...
    'spikeTimesIndex', '/units/spike_times_index', ...
    'spikeTimesNoBlue', '/units/spike_times_no_blue', ...
    'spikeTimesNoBlueIndex', '/units/spike_times_no_blue_index', ...
    'subthresholdVoltage', '/processing/optical_voltage/subthreshold_voltage/data', ...
    'normalizedVoltage', '/processing/optical_voltage/normalized_voltage/data', ...
    'voltageTimestamps', '/processing/optical_voltage/normalized_voltage/timestamps', ...
    'blueEvokedSTA', '/processing/optical_voltage/blue_evoked_spike_triggered_average/data', ...
    'blueStimulation', '/stimulus/presentation/blue_stimulation_per_cell/data', ...
    'runningSpeed', '/processing/behavior/behavioral_time_series/running_speed/data', ...
    'quietMask', '/processing/behavior/behavioral_time_series/quiet_mask/data', ...
    'whiskingMask', '/processing/behavior/behavioral_time_series/whisking_mask/data', ...
    'whiskingMotion', '/processing/behavior/behavioral_time_series/whisking_motion/data', ...
    'whiskingMotionTimestamps', '/processing/behavior/behavioral_time_series/whisking_motion/timestamps', ...
    'chunkStart', '/intervals/recording_chunks/start_time', ...
    'chunkStop', '/intervals/recording_chunks/stop_time');
end

function paths = sparsePulseHadExtraPaths()
base = '/processing/connectivity_mapping/';
paths = struct( ...
    'spontaneousCrossSpikeSTA', [base 'spontaneous_cross_spike_sta/data'], ...
    'spontaneousSpikeCrossCorrelogram', [base 'spontaneous_spike_cross_correlogram/data'], ...
    'normalizedVoltage', '/processing/optical_voltage/normalized_voltage/data', ...
    'voltageTimestamps', '/processing/optical_voltage/normalized_voltage/timestamps', ...
    'blueStimulation', '/stimulus/presentation/blue_stimulation_per_cell/data', ...
    'spikeTimes', '/units/spike_times', ...
    'spikeTimesIndex', '/units/spike_times_index', ...
    'spikeTimesNoBlue', '/units/spike_times_no_blue', ...
    'spikeTimesNoBlueIndex', '/units/spike_times_no_blue_index', ...
    'hadamardSpikeTimes', ...
        '/analysis/pairwise_connectivity_metrics/hadamard_orthogonal_spike_times', ...
    'hadamardSpikeTimesIndex', ...
        '/analysis/pairwise_connectivity_metrics/hadamard_orthogonal_spike_times_index', ...
    'chunkStart', '/intervals/recording_chunks/start_time', ...
    'chunkStop', '/intervals/recording_chunks/stop_time');
end

% =========================================================================
function [counts, endOffset] = readHadamardSpikeCounts(nwbFile, paths, nCells)
% The ragged column stores cumulative end offsets, so the per-pair counts are
% the successive differences. Rows follow the source MATLAB serialization, so
% the presynaptic index varies fastest and a plain reshape restores the matrix.
endOffset = double(h5read(nwbFile, paths.hadamardSpikeTimesIndex));
endOffset = endOffset(:);
assert(numel(endOffset) == nCells * nCells, 'Figure5NWB:HadamardIndexSize', ...
    ['hadamard_orthogonal_spike_times_index must have one row per directed ', ...
     'pair (%d), found %d.'], nCells * nCells, numel(endOffset));
assert(all(diff(endOffset) >= 0), 'Figure5NWB:HadamardIndexOrder', ...
    'The ragged index must be non-decreasing.');
counts = reshape(diff([0; endOffset]), nCells, nCells);
end

function frames = readHadamardPairSpikeFrames(nwbFile, paths, nCells, endOffset, ...
    voltageSeconds, framePeriodSeconds, preIndex, postIndex)
validateattributes(preIndex, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', nCells});
validateattributes(postIndex, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', nCells});
row = sub2ind([nCells nCells], preIndex, postIndex);
stopRow = endOffset(row);
if row == 1
    startRow = 0;
else
    startRow = endOffset(row - 1);
end
if stopRow <= startRow
    frames = zeros(1, 0);
    return
end
seconds = double(h5read(nwbFile, paths.hadamardSpikeTimes, ...
    startRow + 1, stopRow - startRow));
frames = secondsToFrameIndices(seconds, voltageSeconds, framePeriodSeconds, ...
    'hadamard_orthogonal_spike_times');
frames = frames(:)';
end

function values = readFigure5PairwiseSeries(nwbFile, paths, nCells, figure3Reader, ...
    seriesName, preIndex, postIndex, timeRange)
switch lower(strrep(strtrim(seriesName), ' ', '_'))
    case {'spontaneous_cross_spike_sta', 'spontaneouscrossspikesta', 'crossnobluestasubfn'}
        datasetPath = paths.spontaneousCrossSpikeSTA;
    case {'spontaneous_spike_cross_correlogram', 'spontaneousspikecrosscorrelogram', ...
            'crossnobluestarate'}
        datasetPath = paths.spontaneousSpikeCrossCorrelogram;
    otherwise
        % Not one of the Figure 5 additions; let the Figure 3 loader resolve it.
        values = figure3Reader(seriesName, preIndex, postIndex, timeRange);
        return
end
values = readPairwiseDataset(nwbFile, datasetPath, nCells, preIndex, postIndex, timeRange);
end

function values = readPairwiseDataset(nwbFile, datasetPath, nCells, ...
    preIndex, postIndex, timeRange)
% Mirrors the Figure 3 loader's reader, including its orientation contract:
% MATLAB's HDF5 interface restores the source pre x post x time order even
% though a generic reader reports the physical layout with time first.
validateattributes(preIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
validateattributes(postIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
preIndex = preIndex(:)';
postIndex = postIndex(:)';
assert(isequal(preIndex, preIndex(1):preIndex(end)) && ...
    isequal(postIndex, postIndex(1):postIndex(end)), 'Figure5NWB:PairRange', ...
    'HDF5 pairwise reads require contiguous pre- and postsynaptic index ranges.');
info = h5info(nwbFile, datasetPath);
datasetSize = info.Dataspace.Size;
assert(numel(datasetSize) == 3 && datasetSize(1) == nCells && datasetSize(2) == nCells, ...
    'Figure5NWB:PairDatasetLayout', ...
    'Expected MATLAB HDF5 layout pre x post x time at %s.', datasetPath);
nTime = datasetSize(3);
if isempty(timeRange)
    timeRange = [1 nTime];
end
validateattributes(timeRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', nTime});
assert(timeRange(1) <= timeRange(2), 'Figure5NWB:PairTimeRange', ...
    'timeRange must be ascending.');
start = [preIndex(1) postIndex(1) timeRange(1)];
count = [numel(preIndex) numel(postIndex) timeRange(2) - timeRange(1) + 1];
values = h5read(nwbFile, datasetPath, start, count);
end

% =========================================================================
function values = readCellFrameDataset(nwbFile, datasetPath, nCells, nFrames, ...
    cellIndex, frameRange)
validateattributes(cellIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
cellIndex = cellIndex(:)';
assert(isequal(cellIndex, cellIndex(1):cellIndex(end)), 'Figure5NWB:CellRange', ...
    'HDF5 cell-time reads require contiguous cell indices.');
if isempty(frameRange)
    frameRange = [1 nFrames];
end
validateattributes(frameRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', nFrames});
assert(frameRange(1) <= frameRange(2), 'Figure5NWB:FrameRange', ...
    'frameRange must be ascending.');
start = [cellIndex(1) frameRange(1)];
count = [numel(cellIndex) diff(frameRange) + 1];
values = h5read(nwbFile, datasetPath, start, count);
end

function datasetPath = sparseOptoSeriesPath(paths, seriesName)
switch lower(strrep(strtrim(seriesName), ' ', '_'))
    case {'subthreshold_voltage', 'subthresholdvoltage', 'traces_all_sub_n', 'sub'}
        datasetPath = paths.subthresholdVoltage;
    case {'normalized_voltage', 'normalizedvoltage', 'traces_all_n', 'v'}
        datasetPath = paths.normalizedVoltage;
    otherwise
        error('Figure5NWB:UnknownVoltageSeries', ...
            ['Unknown voltage series "%s". Use subthreshold_voltage or ', ...
             'normalized_voltage.'], seriesName);
end
end

function nFramesPerChunk = chunkFrames(nwbFile, datasetPath)
% Reports the HDF5 chunk extent along the frame axis, the natural read-block
% size. Chunks span all cells, so frame blocks are the only cheap slice.
info = h5info(nwbFile, datasetPath);
if isempty(info.ChunkSize)
    nFramesPerChunk = 10000;
    return
end
nFramesPerChunk = info.ChunkSize(2);
end

% =========================================================================
function spikeFrameIndex = readSpikeFrameIndices(nwbFile, timesPath, indexPath, ...
    nCells, voltageSeconds, framePeriodSeconds, label)
% The converters stored seconds as voltageSeconds(sourceFrameIndex). Invert that
% lookup so the returned indices are the 1-based voltage-frame indices the
% source scripts use to index traces.
spikeSeconds = double(h5read(nwbFile, timesPath));
spikeSeconds = spikeSeconds(:);
endOffset = double(h5read(nwbFile, indexPath));
endOffset = endOffset(:);
assert(numel(endOffset) == nCells, 'Figure5NWB:SpikeIndexSize', ...
    '%s_index must have one entry per unit.', label);
assert(endOffset(end) == numel(spikeSeconds), 'Figure5NWB:SpikeIndexTotal', ...
    'The last %s_index entry must equal the total spike count.', label);

frames = secondsToFrameIndices(spikeSeconds, voltageSeconds, ...
    framePeriodSeconds, label);

spikeFrameIndex = cell(nCells, 1);
startRow = 1;
for iCell = 1:nCells
    stopRow = endOffset(iCell);
    spikeFrameIndex{iCell} = frames(startRow:stopRow)';
    startRow = stopRow + 1;
end
end

function frames = secondsToFrameIndices(seconds, voltageSeconds, ...
    framePeriodSeconds, label)
nFrames = numel(voltageSeconds);
frames = interp1(voltageSeconds, (1:nFrames)', seconds(:), 'nearest', 'extrap');
residual = abs(voltageSeconds(frames) - seconds(:));
assert(isempty(residual) || max(residual) < framePeriodSeconds / 2, ...
    'Figure5NWB:FrameLookup', ...
    ['%s could not be matched to voltage frames: largest residual %.3g s ', ...
     'exceeds half a frame (%.3g s).'], label, max(residual), framePeriodSeconds / 2);
end

% =========================================================================
function values = readColumn(nwbFile, datasetPath, expectedLength, label)
values = double(h5read(nwbFile, datasetPath));
values = values(:);
assert(numel(values) == expectedLength, 'Figure5NWB:ColumnLength', ...
    '%s must have %d entries, found %d.', label, expectedLength, numel(values));
end

function values = readRowVector(nwbFile, datasetPath, expectedLength, label)
values = double(h5read(nwbFile, datasetPath));
values = reshape(values, 1, []);
assert(numel(values) == expectedLength, 'Figure5NWB:RowLength', ...
    '%s must have %d samples, found %d.', label, expectedLength, numel(values));
end

function pixelSize = pixelSizeFromUnits(nwbFile, paths, coordinatePixels)
% The converters wrote both pixel and micrometer coordinates, so the recorded
% calibration is recoverable rather than assumed.
xMicrometers = double(h5read(nwbFile, paths.xMicrometers));
xMicrometers = xMicrometers(:);
scale = xMicrometers ./ coordinatePixels(:, 1);
scale = scale(isfinite(scale) & coordinatePixels(:, 1) ~= 0);
assert(~isempty(scale), 'Figure5NWB:PixelSize', ...
    'Could not recover the pixel size from /units/x_um and /units/x_pixel.');
pixelSize = median(scale);
assert(max(abs(scale - pixelSize)) < 1e-6 * max(1, pixelSize), ...
    'Figure5NWB:PixelSizeInconsistent', ...
    'The recorded micrometer and pixel coordinates imply inconsistent scaling.');
end

function tf = datasetExists(nwbFile, datasetPath)
try
    h5info(nwbFile, datasetPath);
    tf = true;
catch
    tf = false;
end
end

function assertDatasetExists(nwbFile, datasetPath)
assert(datasetExists(nwbFile, datasetPath), 'Figure5NWB:MissingDataset', ...
    'Required dataset is missing from %s: %s', nwbFile, datasetPath);
end
