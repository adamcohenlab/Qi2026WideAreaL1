function session = loadFigure6SessionFromNWB(nwbFile)
% LOADFIGURE6SESSIONFROMNWB Load Figure 6 inputs from any of three session types.
%
%   SESSION = LOADFIGURE6SESSIONFROMNWB(NWBFILE) reads one of the NWB files
%   Figure 6 needs and returns ordinary MATLAB arrays plus lazy readers for the
%   large datasets. The session type is detected from the file contents, so the
%   same call works for all three:
%
%     'behavior'        panels E, F, G, H, I, J, K. Spontaneous activity with
%                       no optogenetic stimulation.
%     'sparseOpto'      panels B, C. Intrinsic properties measured under
%                       optogenetic stimulation.
%     'sparsePulseHad'  panel D. Synaptic connectivity.
%
%   Session assignment
%   ------------------
%   Figure 6 spans two mice and three protocols, all confirmed numerically
%   during the Phase 1 audit (see FIGURE6_PHASE1_HANDOFF.md):
%
%     B, C              M-YQ0201-29 sparseOpto
%     D                 M-YQ0201-27 sparsePulseHad   (n = 1420/1818/2224 pairs)
%     E, F, H, I, J, K  M-YQ0201-29 behavior         (n = 141/117 cells,
%                                                     5888 and 547 events)
%     G                 both behavior sessions       (505 cells, 2 mice)
%
%   Relationship to the other loaders
%   ---------------------------------
%   For a sparseOpto or sparsePulseHad file this function **delegates to
%   LOADFIGURE5SESSIONFROMNWB**, which in turn delegates the connectivity layer
%   to LOADFIGURE3SESSIONFROMNWB, and then adds only what Figure 6 needs on top.
%   Everything those loaders return is still present on the returned struct.
%   The behavior mode is new: no earlier figure needed the spontaneous-activity
%   sessions through this interface.
%
%   The loader uses built-in HDF5 functions only, so it is compatible with
%   MATLAB R2019b and does not require MatNWB at reproduction time. Dataset
%   paths follow the embedded MATLAB-to-NWB lookup tables written by
%   convert_behavior_to_nwb.m, convert_sparseOpto_to_nwb.m and
%   convert_sparsePulseHad_to_nwb.m.
%
%   Cell axis and cell selection
%   ----------------------------
%   All three file types store all 320 detected cells and mark the figure
%   scripts' selection in /units/passed_qc. session.figureCellIndices is
%   find(passedQC).
%
%   The behavior rule is  spkHgtNoBlue > 0 & fpRates < 0.5 & idxUseMerge, which
%   is exactly what every Figure 6 behavior-session script applies. The
%   sparseOpto rule is the same. fig6_NPYpm_optoAvgTraces.m (panel B) adds a
%   2-sigma spike-height term that passed_qc does not carry; the author waived
%   it on 2026-09-19 and it was measured to exclude no cell in M-YQ0201-29
%   sparseOpto, so passed_qc is exact there. It is *not* inert in M-YQ0201-27
%   sparseOpto, where it would drop 2 of 286 cells. REPRODUCE_FIGURE6_FROM_NWB
%   re-runs that check and warns if it ever bites.
%
%   NPY labelling
%   -------------
%   /units/gfp_positive is the manuscript's NPY indicator for every session.
%   The legacy scripts override gfp_idx with gfpMaskSD from
%   SD_to_FF_registration.mat when that file exists, which it does for all three
%   M-YQ0201-29 sessions. gfpMaskSD was verified element-for-element identical
%   to the deposited gfp_positive in all three, so there is nothing to override
%   here. See FIGURE6_PHASE1_HANDOFF.md.
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
%   Spike times  : /units/spike_times is in seconds. This loader converts them
%                  back to the 1-based voltage-frame indices the source scripts
%                  index traces with, in spikeFrameIndex. Behavior sessions have
%                  no blue light and therefore no spike_times_no_blue.
%   Masks        : quiet_mask and whisking_mask are stored as uint8 and are
%                  returned as logical. They are the Gaussian-mixture rule,
%                  which is what panels I, J and K use. Panel H's state bar uses
%                  a different, median-threshold rule that the reproduction
%                  recomputes from whiskingMotion.
%
%   behavior fields (panels E to K)
%   -------------------------------
%   whiskingMotion                  [1 x nFace] source mWhisk, face timebase.
%   whiskingMotionReconstruction    [nFace x 15] source mWhisk_rec, the sparse
%                                   deconvolution at 15 regularisation
%                                   strengths. Added to both behavior NWB files
%                                   on 2026-09-19. **h5read returns it as
%                                   [15 x nFace]; this loader transposes it back
%                                   to the source time-by-repeat order.** The
%                                   regularisation constants themselves are not
%                                   stored, so a lambda can only be named by its
%                                   column index, which is all the analysis uses.
%   framesPerChunk                  288000 for behavior. The recording is 10
%                                   such chunks.
%   readVoltageRows                 Lazy reader; see the performance note below.
%
%   sparsePulseHad extra field (panel D)
%   ------------------------------------
%   gfpPositive                     The Figure 3 loader, which supplies the
%                                   connectivity layer, does not read
%                                   /units/gfp_positive. Panel D needs it, so it
%                                   is added here.
%
%   Reading the voltage datasets efficiently
%   ----------------------------------------
%   normalized_voltage and subthreshold_voltage are [nCells x nFrames] in MATLAB
%   order and HDF5-chunked so that every chunk spans *all* cells for a window of
%   frames. Reading a block of frames for all cells is contiguous and fast;
%   reading one cell's full row costs the same as reading all cells, because it
%   must still touch every chunk. Iterate over frame blocks, not cells. Every
%   Figure 6 quantity that needs the full recording is a histogram, a mean or a
%   triggered-average accumulation over time, so block-wise accumulation is
%   exact. session.voltageChunkFrames reports the chunk extent along the frame
%   axis.
%
%   Example
%   -------
%     beh = loadFigure6SessionFromNWB('...\M_YQ0201_29_behavior.nwb');
%     npy = beh.gfpPositive(beh.figureCellIndices);       % 117 true of 258
%     block = beh.readVoltageRows('subthreshold_voltage', ...
%         1:beh.numberOfCells, [1 288000]);
%
%     opto = loadFigure6SessionFromNWB('...\M_YQ0201_29_sparseOpto.nwb');
%     adp = opto.intrinsicProperties(opto.figureCellIndices, 1);
%
%     conn = loadFigure6SessionFromNWB('...\M_YQ0201_27_sparsePulseHad.nwb');
%     ipsp = conn.readPairwiseSeries('blue_pulse_ipsp', 1:conn.numberOfCells, ...
%         1:conn.numberOfCells, []);
%
%   See also REPRODUCE_FIGURE6_FROM_NWB, LOADFIGURE5SESSIONFROMNWB,
%   LOADFIGURE3SESSIONFROMNWB.

if nargin < 1 || isempty(nwbFile)
    error('Figure6NWB:MissingFile', ...
        'Provide a behavior, sparseOpto or sparsePulseHad NWB file.');
end
validateattributes(nwbFile, {'char'}, {'row', 'nonempty'}, mfilename, 'nwbFile');
assert(exist(nwbFile, 'file') == 2, 'Figure6NWB:FileNotFound', ...
    'NWB file not found: %s', nwbFile);

sessionType = detectSessionType(nwbFile);
switch sessionType
    case 'behavior'
        session = loadBehaviorSession(nwbFile);
    case {'sparseOpto', 'sparsePulseHad'}
        session = loadStimulationSession(nwbFile, sessionType);
    otherwise
        error('Figure6NWB:UnknownSessionType', ...
            'Unrecognized session type: %s', sessionType);
end
session.sessionType = sessionType;
session.nwbFile = nwbFile;
end

% =========================================================================
function sessionType = detectSessionType(nwbFile)
% Probe for the marker dataset of each protocol rather than guessing from the
% file name. The intrinsic-property columns exist only in a sparseOpto file, the
% pairwise connectivity table only in a sparsePulseHad file, and the
% state-resolved correlation module only in a behavior file.
hasProperties = datasetExists(nwbFile, '/units/optical_rheobase_norm');
hasPairwise = datasetExists(nwbFile, ...
    '/analysis/pairwise_connectivity_metrics/synaptic_connection');
hasStateCorrelations = datasetExists(nwbFile, ...
    '/processing/state_dependent_correlations/spike_triggered_voltage_quiet/data');

nMatched = hasProperties + hasPairwise + hasStateCorrelations;
assert(nMatched > 0, 'Figure6NWB:UnrecognizedFile', ...
    ['%s is none of the three Figure 6 session types: no ', ...
     '/units/optical_rheobase_norm (sparseOpto), no pairwise connectivity ', ...
     'table (sparsePulseHad), and no state-dependent correlations (behavior).'], ...
    nwbFile);
assert(nMatched == 1, 'Figure6NWB:AmbiguousFile', ...
    '%s matches more than one session type at once.', nwbFile);

if hasProperties
    sessionType = 'sparseOpto';
elseif hasPairwise
    sessionType = 'sparsePulseHad';
else
    sessionType = 'behavior';
end
end

% =========================================================================
function session = loadStimulationSession(nwbFile, sessionType)
% Panels B, C and D reuse the Figure 5 loader wholesale. Only one field is
% missing for Figure 6's purposes.
assert(~isempty(which('loadFigure5SessionFromNWB')), 'Figure6NWB:MissingFigure5Loader', ...
    ['loadFigure5SessionFromNWB.m must be on the path; Figure 6 panels B, C ', ...
     'and D reuse it for the stimulation-session layer.']);
session = loadFigure5SessionFromNWB(nwbFile);

if ~isfield(session, 'gfpPositive')
    % The sparsePulseHad branch of the Figure 5 loader delegates to the
    % Figure 3 loader, which has no use for the NPY indicator and does not read
    % it. Panel D needs it.
    session.gfpPositive = logical(readColumn(nwbFile, '/units/gfp_positive', ...
        session.numberOfCells, 'gfp_positive'));
end
session.gfpPositive = session.gfpPositive(:);
session.npyLabels = {'NPY-', 'NPY+'};
session.gfpProvenanceNote = [ ...
    '/units/gfp_positive is the manuscript NPY indicator. For M-YQ0201-29 ', ...
    'sessions it is element-for-element identical to gfpMaskSD in ', ...
    'SD_to_FF_registration.mat, which the legacy scripts substitute for ', ...
    'gfp_idx, so no override is needed here.'];

if ~isfield(session, 'cellDistancePixels')
    session.cellDistancePixels = squareform(pdist(session.cellCoordinatePixels));
end
if strcmp(sessionType, 'sparsePulseHad') && ~isfield(session, 'pixelSizeMicrometers')
    session.pixelSizeMicrometers = 6.5;
end
end

% =========================================================================
function session = loadBehaviorSession(nwbFile)
paths = behaviorPaths();
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
assert(session.numberOfCells > 0, 'Figure6NWB:NoCells', 'The units table is empty.');

session.cellCoordinatePixels = [ ...
    readColumn(nwbFile, paths.xPixel, session.numberOfCells, 'x_pixel'), ...
    readColumn(nwbFile, paths.yPixel, session.numberOfCells, 'y_pixel')];
session.passedQC = logical(readColumn(nwbFile, paths.passedQC, ...
    session.numberOfCells, 'passed_qc'));
session.gfpPositive = logical(readColumn(nwbFile, paths.gfpPositive, ...
    session.numberOfCells, 'gfp_positive'));
session.figureCellIndices = find(session.passedQC);
session.npyLabels = {'NPY-', 'NPY+'};
session.gfpProvenanceNote = [ ...
    '/units/gfp_positive is the manuscript NPY indicator. For M-YQ0201-29 it ', ...
    'is element-for-element identical to gfpMaskSD in ', ...
    'SD_to_FF_registration.mat, which the legacy scripts substitute for ', ...
    'gfp_idx, so no override is needed here.'];
session.pixelSizeMicrometers = pixelSizeFromUnits(nwbFile, paths, ...
    session.cellCoordinatePixels);
session.cellDistancePixels = squareform(pdist(session.cellCoordinatePixels));

% --------------------------------------------------------------- timebase
session.voltageTimestampsSeconds = double(h5read(nwbFile, paths.voltageTimestamps));
session.voltageTimestampsSeconds = session.voltageTimestampsSeconds(:);
session.numberOfVoltageFrames = numel(session.voltageTimestampsSeconds);
session.framePeriodSeconds = 1.27e-3;

voltageInfo = h5info(nwbFile, paths.subthresholdVoltage);
assert(isequal(voltageInfo.Dataspace.Size, ...
    [session.numberOfCells session.numberOfVoltageFrames]), ...
    'Figure6NWB:VoltageLayout', ...
    'Expected MATLAB HDF5 layout cell x frame at %s.', paths.subthresholdVoltage);

session.recordingChunkStartSeconds = double(h5read(nwbFile, paths.chunkStart));
session.recordingChunkStopSeconds = double(h5read(nwbFile, paths.chunkStop));
session.numberOfChunks = numel(session.recordingChunkStartSeconds);
assert(session.numberOfChunks > 0 && ...
    mod(session.numberOfVoltageFrames, session.numberOfChunks) == 0, ...
    'Figure6NWB:ChunkStructure', ...
    'The voltage frame count is not a whole multiple of the chunk count.');
session.framesPerChunk = session.numberOfVoltageFrames / session.numberOfChunks;

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
assert(mod(session.numberOfFaceFrames, session.numberOfChunks) == 0, ...
    'Figure6NWB:FaceChunkStructure', ...
    'The face-video frame count is not a whole multiple of the chunk count.');
session.faceFramesPerChunk = session.numberOfFaceFrames / session.numberOfChunks;

session.whiskingMotionReconstruction = readWhiskingReconstruction(nwbFile, ...
    paths.whiskingReconstruction, session.numberOfFaceFrames);
session.numberOfWhiskingLambdas = size(session.whiskingMotionReconstruction, 2);
session.whiskingReconstructionNote = [ ...
    'Source mWhisk_rec, a sparse deconvolution of mWhisk at 15 regularisation ', ...
    'strengths, on the face timebase. Added to both behavior NWB files on ', ...
    '2026-09-19. h5read returns it as [nLambda x nFace]; it is transposed here ', ...
    'to the source time-by-repeat order. The regularisation constants are not ', ...
    'stored, so a lambda is identified only by its column index, which is all ', ...
    'the analysis uses.'];

session.maskProvenanceNote = [ ...
    'quiet_mask and whisking_mask follow the Gaussian-mixture rule of ', ...
    'dataPrepsForMAT.m lines 366-420, the same rule panels I, J and K use. ', ...
    'They came from an unseeded fitgmdist draw, so they are the same ', ...
    'definition as, but not necessarily bit-identical to, the masks in the ', ...
    'legacy run. Panel H uses a different median-threshold rule; see ', ...
    'reproduce_figure6_from_nwb.m.'];

% ----------------------------------------------------------- spike times
% Behavior sessions delivered no blue light, so there is no spike_times_no_blue.
session.spikeFrameIndex = readSpikeFrameIndices(nwbFile, paths.spikeTimes, ...
    paths.spikeTimesIndex, session.numberOfCells, ...
    session.voltageTimestampsSeconds, session.framePeriodSeconds, 'spike_times');
session.spikeCount = cellfun(@numel, session.spikeFrameIndex);

% -------------------------------------------------------------- readers
session.readVoltageRows = @(seriesName, cellIndex, frameRange) ...
    readCellFrameDataset(nwbFile, behaviorSeriesPath(paths, seriesName), ...
    session.numberOfCells, session.numberOfVoltageFrames, cellIndex, frameRange);
session.voltageChunkFrames = chunkFrames(nwbFile, paths.subthresholdVoltage);
end

% =========================================================================
function paths = behaviorPaths()
behavior = '/processing/behavior/behavioral_time_series/';
paths = struct( ...
    'cellId', '/units/id', ...
    'xPixel', '/units/x_pixel', ...
    'yPixel', '/units/y_pixel', ...
    'xMicrometers', '/units/x_um', ...
    'passedQC', '/units/passed_qc', ...
    'gfpPositive', '/units/gfp_positive', ...
    'spikeTimes', '/units/spike_times', ...
    'spikeTimesIndex', '/units/spike_times_index', ...
    'subthresholdVoltage', '/processing/optical_voltage/subthreshold_voltage/data', ...
    'normalizedVoltage', '/processing/optical_voltage/normalized_voltage/data', ...
    'voltageTimestamps', '/processing/optical_voltage/normalized_voltage/timestamps', ...
    'runningSpeed', [behavior 'running_speed/data'], ...
    'quietMask', [behavior 'quiet_mask/data'], ...
    'whiskingMask', [behavior 'whisking_mask/data'], ...
    'whiskingMotion', [behavior 'whisking_motion/data'], ...
    'whiskingMotionTimestamps', [behavior 'whisking_motion/timestamps'], ...
    'whiskingReconstruction', [behavior 'whisking_motion_reconstruction/data'], ...
    'chunkStart', '/intervals/recording_chunks/start_time', ...
    'chunkStop', '/intervals/recording_chunks/stop_time');
end

% =========================================================================
function reconstruction = readWhiskingReconstruction(nwbFile, datasetPath, nFaceFrames)
% Stored (nFace, nLambda) on disk. MATLAB's HDF5 interface reverses the axes, so
% h5read hands back [nLambda x nFace]; the source mWhisk_rec is time-by-repeat,
% so transpose. Getting this wrong does not error, it silently selects a
% different regularisation strength, so the orientation is asserted.
raw = double(h5read(nwbFile, datasetPath));
assert(ismatrix(raw), 'Figure6NWB:WhiskReconstructionRank', ...
    'Expected a 2-D whisking reconstruction at %s.', datasetPath);
if size(raw, 2) == nFaceFrames && size(raw, 1) ~= nFaceFrames
    reconstruction = raw';
elseif size(raw, 1) == nFaceFrames
    reconstruction = raw;
else
    error('Figure6NWB:WhiskReconstructionShape', ...
        ['%s is %d x %d; one axis must equal the %d face-video frames so the ', ...
         'source time-by-repeat order can be restored.'], datasetPath, ...
        size(raw, 1), size(raw, 2), nFaceFrames);
end
assert(size(reconstruction, 2) >= 1 && size(reconstruction, 2) < nFaceFrames, ...
    'Figure6NWB:WhiskReconstructionLambdas', ...
    'The whisking reconstruction should have far fewer repeats than samples.');
end

% =========================================================================
function values = readCellFrameDataset(nwbFile, datasetPath, nCells, nFrames, ...
    cellIndex, frameRange)
validateattributes(cellIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
cellIndex = cellIndex(:)';
assert(isequal(cellIndex, cellIndex(1):cellIndex(end)), 'Figure6NWB:CellRange', ...
    'HDF5 cell-time reads require contiguous cell indices.');
if isempty(frameRange)
    frameRange = [1 nFrames];
end
validateattributes(frameRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', nFrames});
assert(frameRange(1) <= frameRange(2), 'Figure6NWB:FrameRange', ...
    'frameRange must be ascending.');
start = [cellIndex(1) frameRange(1)];
count = [numel(cellIndex) diff(frameRange) + 1];
values = h5read(nwbFile, datasetPath, start, count);
end

function datasetPath = behaviorSeriesPath(paths, seriesName)
switch lower(strrep(strtrim(seriesName), ' ', '_'))
    case {'subthreshold_voltage', 'subthresholdvoltage', 'traces_all_sub_n', 'sub'}
        datasetPath = paths.subthresholdVoltage;
    case {'normalized_voltage', 'normalizedvoltage', 'traces_all_n', 'v'}
        datasetPath = paths.normalizedVoltage;
    otherwise
        error('Figure6NWB:UnknownVoltageSeries', ...
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
assert(numel(endOffset) == nCells, 'Figure6NWB:SpikeIndexSize', ...
    '%s_index must have one entry per unit.', label);
assert(endOffset(end) == numel(spikeSeconds), 'Figure6NWB:SpikeIndexTotal', ...
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
    'Figure6NWB:FrameLookup', ...
    ['%s could not be matched to voltage frames: largest residual %.3g s ', ...
     'exceeds half a frame (%.3g s).'], label, max(residual), framePeriodSeconds / 2);
end

% =========================================================================
function values = readColumn(nwbFile, datasetPath, expectedLength, label)
values = double(h5read(nwbFile, datasetPath));
values = values(:);
assert(numel(values) == expectedLength, 'Figure6NWB:ColumnLength', ...
    '%s must have %d entries, found %d.', label, expectedLength, numel(values));
end

function values = readRowVector(nwbFile, datasetPath, expectedLength, label)
values = double(h5read(nwbFile, datasetPath));
values = reshape(values, 1, []);
assert(numel(values) == expectedLength, 'Figure6NWB:RowLength', ...
    '%s must have %d samples, found %d.', label, expectedLength, numel(values));
end

function pixelSize = pixelSizeFromUnits(nwbFile, paths, coordinatePixels)
% The converters wrote both pixel and micrometer coordinates, so the recorded
% calibration is recoverable rather than assumed.
xMicrometers = double(h5read(nwbFile, paths.xMicrometers));
xMicrometers = xMicrometers(:);
scale = xMicrometers ./ coordinatePixels(:, 1);
scale = scale(isfinite(scale) & coordinatePixels(:, 1) ~= 0);
assert(~isempty(scale), 'Figure6NWB:PixelSize', ...
    'Could not recover the pixel size from /units/x_um and /units/x_pixel.');
pixelSize = median(scale);
assert(max(abs(scale - pixelSize)) < 1e-6 * max(1, pixelSize), ...
    'Figure6NWB:PixelSizeInconsistent', ...
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
assert(datasetExists(nwbFile, datasetPath), 'Figure6NWB:MissingDataset', ...
    'Required dataset is missing from %s: %s', nwbFile, datasetPath);
end
