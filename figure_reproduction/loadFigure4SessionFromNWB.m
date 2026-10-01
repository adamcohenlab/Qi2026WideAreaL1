function session = loadFigure4SessionFromNWB(nwbFile)
% LOADFIGURE4SESSIONFROMNWB Load Figure 4 spontaneous-activity inputs from NWB.
%
%   SESSION = LOADFIGURE4SESSIONFROMNWB(NWBFILE) reads a behavior
%   (spontaneous-activity) NWB file, eagerly loads the small per-cell and
%   behavioral arrays, and returns lazy readers for the large voltage and
%   correlation datasets. Source MATLAB array orientation is preserved on
%   output.
%
%   The loader uses built-in HDF5 functions only, so it is compatible with
%   MATLAB R2019b and does not require MatNWB at reproduction time. Dataset
%   paths follow the embedded MATLAB-to-NWB lookup table written by
%   convert_behavior_to_nwb.m and recorded in
%   M_YQ0201_27_behavior_NWB_progress.md and
%   M_YQ0201_29_behavior_NWB_progress.md.
%
%   Two cell axes
%   -------------
%   This is the one thing most likely to be misused. The file has two
%   different cell axes and they are not interchangeable:
%
%     * /units has one row per *detected* cell (320 in both sessions).
%       session.numberOfCells, session.cellCoordinatePixels,
%       session.passedQC, and session.spikeFrameIndex are on this axis.
%
%     * The state-resolved correlation arrays span only the cells that
%       passed quality control (247 for M-YQ0201-27, 258 for M-YQ0201-29).
%       Their indices are NOT unit IDs. session.correlationCellCount,
%       session.correlationCoordinatePixels, and
%       session.correlationDistancePixels are on this axis.
%
%   session.correlationUnitIndex maps correlation index -> 1-based /units
%   row, and session.unitToCorrelationIndex maps the other way (0 where a
%   unit is not in the correlation set). The mapping is read from
%   /analysis/correlation_cell_index and cross-checked against passed_qc.
%
%   The legacy figure scripts apply their cell-selection rule before any
%   analysis, so the correlation axis is the axis the Figure 4 panels use.
%   Panels A, E, F, and G must therefore be restricted to
%   session.correlationUnitIndex to match the published analysis, even
%   though the voltage traces are stored for all 320 cells.
%
%   Units and conventions
%   ---------------------
%   Voltage      : dimensionless, already divided by each cell's mean
%                  no-blue spike height ("SH" in the source scripts).
%   Coordinates  : session.cellCoordinatePixels is in camera pixels, as in
%                  the source cell_coord. Multiply by
%                  session.pixelSizeMicrometers (6.5 um/pixel) for microns.
%                  Distances returned by this loader are in pixels, which is
%                  the unit the source distance bins are expressed in.
%   Time         : NWB stores seconds. session.framePeriodSeconds is
%                  1.27e-3 s, the voltage camera exposure, which is the
%                  source scripts' dt.
%   Spike times  : /units/spike_times is in seconds. This loader converts
%                  them back to the 1-based voltage-frame indices used by
%                  the source scripts; see session.spikeFrameIndex.
%   Lag axis     : session.lagFrames is -nLag:nLag and
%                  session.zeroLagIndex is nLag+1 (141 for nLag = 140).
%   Masks        : quiet_mask and whisking_mask are stored as uint8 and are
%                  returned as logical.
%
%   Correlation array orientation
%   -----------------------------
%   Each array is returned as [trigger cell x target cell x lag], the source
%   MATLAB orientation: the row is the triggering/spiking cell and the
%   column is the cell whose voltage or spike train was averaged
%   (get_sta_mat_2 documents "column triggered by row"). MATLAB's HDF5
%   interface restores this order even though a generic HDF5 reader reports
%   the physical layout with the 281-sample lag axis first. The axis order
%   was also confirmed numerically: for spike_triggered_spike_rate at zero
%   lag, M(i,j)/M(j,i) tracks n_j/n_i, which only holds if the first index
%   is the trigger whose spike count normalizes the average.
%
%   Already-applied normalizations
%   ------------------------------
%   The stored arrays are the *normalized* products of fig4_xcorr.m, not raw
%   sums. Do not re-apply these steps:
%     * voltage_cross_correlation is divided by the per-state lag-1
%       standard deviation of each cell, and its zero-lag autocorrelation
%       diagonal has already been patched from the adjacent lag. Diagonals
%       at lag indices 140, 141, and 142 are all 1 to within 3e-6.
%     * spike_triggered_voltage is divided by the target cell's spike
%       height.
%
%   Reading the voltage datasets efficiently
%   ----------------------------------------
%   normalized_voltage and subthreshold_voltage are [nCells x nFrames] in
%   MATLAB order and HDF5-chunked as [nCells x 10000], so every chunk spans
%   *all* cells for a 10000-frame window. Consequences:
%
%     * Reading a block of frames for all cells is contiguous and fast,
%       measured at about 900 MB/s, so the whole 3.7 GB dataset streams in a
%       few seconds.
%     * Reading a single cell's full row costs the same as reading all cells,
%       because it must still touch every chunk. Looping over cells therefore
%       re-reads the entire dataset once per cell. For 247 cells that is about
%       13 minutes per session instead of 4 seconds.
%
%   So iterate over frame blocks with session.readVoltageRows(name, 1:nCells,
%   [firstFrame lastFrame]) and accumulate across blocks. Every Figure 4
%   quantity that needs the full recording (panels E, F, and G) is a histogram
%   accumulation over time, so block-wise accumulation is exact.
%   session.voltageChunkFrames reports the chunk extent along the frame axis.
%
%   Example
%   -------
%     s = loadFigure4SessionFromNWB('...\M_YQ0201_27_behavior.nwb');
%     xcorrQuiet = s.readCorrelationSeries('voltage_cross_correlation', 'quiet');
%     staRow = s.readCorrelationTrigger('spike_triggered_voltage', 'quiet', 209);
%
%   See also REPRODUCE_FIGURE4_FROM_NWB, LOADFIGURE3SESSIONFROMNWB.

if nargin < 1 || isempty(nwbFile)
    error('Figure4NWB:MissingFile', 'Provide a behavior (spontaneous-activity) NWB file.');
end
validateattributes(nwbFile, {'char'}, {'row', 'nonempty'}, mfilename, 'nwbFile');
assert(exist(nwbFile, 'file') == 2, 'Figure4NWB:FileNotFound', ...
    'NWB file not found: %s', nwbFile);

paths = figure4DatasetPaths();
requiredPaths = {paths.cellId, paths.xPixel, paths.yPixel, paths.passedQC, ...
    paths.gfpPositive, paths.spikeTimes, paths.spikeTimesIndex, ...
    paths.correlationIndex, paths.correlationUnitId, ...
    paths.subthresholdVoltage, paths.normalizedVoltage, paths.voltageTimestamps, ...
    paths.runningSpeed, paths.whiskingMotion, paths.whiskingMotionTimestamps, ...
    paths.quietMask, paths.whiskingMask, paths.chunkStart, paths.chunkStop};
for iPath = 1:numel(requiredPaths)
    assertDatasetExists(nwbFile, requiredPaths{iPath});
end

session = struct();
session.nwbFile = nwbFile;
session.paths = paths;

% ---------------------------------------------------------------- units axis
session.cellId = double(h5read(nwbFile, paths.cellId));
session.cellId = session.cellId(:);
session.numberOfCells = numel(session.cellId);
assert(session.numberOfCells > 0, 'Figure4NWB:NoCells', 'The units table is empty.');
session.cellCoordinatePixels = [double(h5read(nwbFile, paths.xPixel)), ...
    double(h5read(nwbFile, paths.yPixel))];
assert(size(session.cellCoordinatePixels, 1) == session.numberOfCells, ...
    'Figure4NWB:CoordinateSize', 'Coordinate columns must have one value per unit.');
session.passedQC = logical(h5read(nwbFile, paths.passedQC));
session.passedQC = session.passedQC(:);
session.gfpPositive = logical(h5read(nwbFile, paths.gfpPositive));
session.gfpPositive = session.gfpPositive(:);
assert(numel(session.passedQC) == session.numberOfCells && ...
    numel(session.gfpPositive) == session.numberOfCells, ...
    'Figure4NWB:UnitColumnSize', 'Unit columns must have one value per unit.');

% dx is recorded in the NWB lookup-table notes and implied by x_um/y_um. The
% source scripts use this calibration to convert cell_coord to micrometers.
session.pixelSizeMicrometers = 6.5;

% ---------------------------------------------------------- correlation axis
[session.correlationUnitIndex, session.unitToCorrelationIndex] = ...
    readCorrelationCellIndex(nwbFile, paths, session.passedQC);
session.correlationCellCount = numel(session.correlationUnitIndex);
session.correlationCoordinatePixels = ...
    session.cellCoordinatePixels(session.correlationUnitIndex, :);
session.correlationDistancePixels = ...
    squareform(pdist(session.correlationCoordinatePixels));

% ------------------------------------------------------------------ timebase
[session.framePeriodSeconds, session.numberOfLags, lagStartSeconds] = ...
    readLagTimebase(nwbFile, paths);
assert(mod(session.numberOfLags, 2) == 1, 'Figure4NWB:LagParity', ...
    'The lag axis must have an odd number of samples so that zero lag is centered.');
session.numberOfLagFrames = (session.numberOfLags - 1) / 2;
session.lagFrames = (-session.numberOfLagFrames:session.numberOfLagFrames)';
session.zeroLagIndex = session.numberOfLagFrames + 1;

% Build the lag axis as lagFrames * dt, matching the source's tau = (-nLag:nLag)*dt,
% rather than accumulating from starting_time. MatNWB stores a TimeSeries rate as
% single precision, so starting_time + (0:n-1)/rate drifts a few nanoseconds off
% center; the stored starting_time is checked against the exact axis instead.
session.lagSeconds = session.lagFrames * session.framePeriodSeconds;
assert(abs(lagStartSeconds - session.lagSeconds(1)) < session.framePeriodSeconds / 2, ...
    'Figure4NWB:ZeroLag', ...
    ['The stored starting_time (%.9g s) is not within half a frame of ', ...
     '-nLag*dt (%.9g s), so zero lag cannot be placed at the centered sample.'], ...
    lagStartSeconds, session.lagSeconds(1));

session.voltageTimestampsSeconds = double(h5read(nwbFile, paths.voltageTimestamps));
session.voltageTimestampsSeconds = session.voltageTimestampsSeconds(:);
session.numberOfVoltageFrames = numel(session.voltageTimestampsSeconds);
assert(all(diff(session.voltageTimestampsSeconds) > 0), 'Figure4NWB:VoltageClock', ...
    'Voltage timestamps must be strictly increasing.');

% ------------------------------------------------------- recording structure
chunkStart = double(h5read(nwbFile, paths.chunkStart));
chunkStop = double(h5read(nwbFile, paths.chunkStop));
session.recordingChunkStartSeconds = chunkStart(:);
session.recordingChunkStopSeconds = chunkStop(:);
session.numberOfChunks = numel(chunkStart);
assert(session.numberOfChunks > 0, 'Figure4NWB:NoChunks', ...
    'The recording_chunks table is empty.');
% The source scripts hard-code nFrame = 288000 frames per acquisition chunk and
% use it to place the chunk-boundary exclusion guard. Derive it instead, so the
% guard follows the file rather than a literal.
framesPerChunk = session.numberOfVoltageFrames / session.numberOfChunks;
assert(framesPerChunk == floor(framesPerChunk), 'Figure4NWB:ChunkFrames', ...
    'Voltage frames do not divide evenly into %d recording chunks.', session.numberOfChunks);
session.framesPerChunk = framesPerChunk;

% ------------------------------------------------------------------ behavior
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
session.whiskingMotionTimestampsSeconds = session.whiskingMotionTimestampsSeconds(:);
session.numberOfFaceFrames = numel(session.whiskingMotion);
assert(numel(session.whiskingMotionTimestampsSeconds) == session.numberOfFaceFrames, ...
    'Figure4NWB:FaceClock', 'whisking_motion and its timestamps differ in length.');

% The stored masks were written by dataPrepsForMAT.m, which fits its own
% unseeded two-component Gaussian mixture to log(mWhisk) to find the whisking
% threshold. They are therefore the same *definition* as the masks used to
% build the stored correlation arrays, but not necessarily the same draw.
session.maskProvenanceNote = ['quiet_mask and whisking_mask come from an ', ...
    'independent unseeded fitgmdist fit in dataPrepsForMAT.m; the definition ', ...
    'matches fig4_xcorr.m but the fit is a separate draw.'];

% -------------------------------------------------------------- spike timing
session.spikeFrameIndex = readSpikeFrameIndices(nwbFile, paths, ...
    session.numberOfCells, session.voltageTimestampsSeconds, ...
    session.framePeriodSeconds);
session.spikeCount = cellfun(@numel, session.spikeFrameIndex);

% ------------------------------------------------------------- lazy readers
session.readCorrelationSeries = @(seriesName, state) ...
    readCorrelationDataset(nwbFile, paths, session.correlationCellCount, ...
    session.numberOfLags, seriesName, state, [], [], []);
session.readCorrelationTrigger = @(seriesName, state, triggerIndex) ...
    squeeze(readCorrelationDataset(nwbFile, paths, session.correlationCellCount, ...
    session.numberOfLags, seriesName, state, ...
    [triggerIndex triggerIndex], [1 session.correlationCellCount], []));
session.readCorrelationLag = @(seriesName, state, lagIndex) ...
    squeeze(readCorrelationDataset(nwbFile, paths, session.correlationCellCount, ...
    session.numberOfLags, seriesName, state, [], [], [lagIndex lagIndex]));
session.readVoltageRows = @(seriesName, cellIndex, frameRange) ...
    readVoltageDataset(nwbFile, paths, session.numberOfCells, ...
    session.numberOfVoltageFrames, seriesName, cellIndex, frameRange);
session.voltageChunkFrames = voltageChunkFrames(nwbFile, paths);
end

% =========================================================================

function paths = figure4DatasetPaths()
correlations = '/processing/state_dependent_correlations';
behavior = '/processing/behavior/behavioral_time_series';
paths = struct( ...
    'cellId', '/units/id', ...
    'xPixel', '/units/x_pixel', ...
    'yPixel', '/units/y_pixel', ...
    'passedQC', '/units/passed_qc', ...
    'gfpPositive', '/units/gfp_positive', ...
    'spikeTimes', '/units/spike_times', ...
    'spikeTimesIndex', '/units/spike_times_index', ...
    'correlationIndex', '/analysis/correlation_cell_index/correlation_index', ...
    'correlationUnitId', '/analysis/correlation_cell_index/unit_id', ...
    'subthresholdVoltage', '/processing/optical_voltage/subthreshold_voltage/data', ...
    'normalizedVoltage', '/processing/optical_voltage/normalized_voltage/data', ...
    'voltageTimestamps', '/processing/optical_voltage/normalized_voltage/timestamps', ...
    'runningSpeed', [behavior '/running_speed/data'], ...
    'whiskingMotion', [behavior '/whisking_motion/data'], ...
    'whiskingMotionTimestamps', [behavior '/whisking_motion/timestamps'], ...
    'quietMask', [behavior '/quiet_mask/data'], ...
    'whiskingMask', [behavior '/whisking_mask/data'], ...
    'chunkStart', '/intervals/recording_chunks/start_time', ...
    'chunkStop', '/intervals/recording_chunks/stop_time', ...
    'voltageCrossCorrelation', [correlations '/voltage_cross_correlation'], ...
    'spikeTriggeredVoltage', [correlations '/spike_triggered_voltage'], ...
    'spikeTriggeredSpikeRate', [correlations '/spike_triggered_spike_rate']);
end

function [correlationUnitIndex, unitToCorrelationIndex] = ...
    readCorrelationCellIndex(nwbFile, paths, passedQC)
% The correlation arrays span only the QC cells. This table is the documented
% mapping from correlation index to unit ID, written by convert_behavior_to_nwb.m
% as unit_id = find(idxUse) - 1 in ascending unit order.
correlationIndex = double(h5read(nwbFile, paths.correlationIndex));
unitId = double(h5read(nwbFile, paths.correlationUnitId));
correlationIndex = correlationIndex(:);
unitId = unitId(:);
nCorrelation = numel(correlationIndex);
assert(numel(unitId) == nCorrelation, 'Figure4NWB:CorrelationTable', ...
    'correlation_index and unit_id differ in length.');
assert(isequal(correlationIndex, (0:nCorrelation-1)'), 'Figure4NWB:CorrelationOrder', ...
    'correlation_index must be 0..n-1 in ascending order.');
correlationUnitIndex = unitId + 1;
assert(all(correlationUnitIndex >= 1 & correlationUnitIndex <= numel(passedQC)), ...
    'Figure4NWB:CorrelationUnitId', 'Correlation unit IDs fall outside the units table.');
assert(isequal(correlationUnitIndex, sort(correlationUnitIndex)), ...
    'Figure4NWB:CorrelationSort', 'Correlation unit IDs must ascend.');
% Cross-check against passed_qc rather than trusting either source alone.
assert(isequal(correlationUnitIndex, find(passedQC)), 'Figure4NWB:CorrelationQC', ...
    ['correlation_cell_index does not match the units passed_qc column. ', ...
     'The correlation cell axis cannot be resolved.']);
unitToCorrelationIndex = zeros(numel(passedQC), 1);
unitToCorrelationIndex(correlationUnitIndex) = 1:nCorrelation;
end

function [framePeriodSeconds, nLagSamples, lagStartSeconds] = readLagTimebase(nwbFile, paths)
dataPath = [paths.voltageCrossCorrelation '_quiet/data'];
startingTimePath = [paths.voltageCrossCorrelation '_quiet/starting_time'];
assertDatasetExists(nwbFile, dataPath);
assertDatasetExists(nwbFile, startingTimePath);
lagStartSeconds = double(h5read(nwbFile, startingTimePath));
try
    rateHz = double(h5readatt(nwbFile, startingTimePath, 'rate'));
catch
    error('Figure4NWB:MissingRate', ...
        'Correlation TimeSeries starting_time lacks its required rate attribute: %s', ...
        startingTimePath);
end
assert(isscalar(rateHz) && isfinite(rateHz) && rateHz > 0, ...
    'Figure4NWB:Rate', 'Correlation TimeSeries rate must be positive.');

% Prefer the voltage-camera exposure recorded on /units/spike_times: it is
% written as an exact double (1.27e-3 s), whereas the TimeSeries rate is stored
% as single and only recovers dt to about 4e-11 s. Both are cross-checked so a
% future file with a genuinely different sampling interval cannot slip through.
try
    framePeriodSeconds = double(h5readatt(nwbFile, paths.spikeTimes, 'resolution'));
catch
    framePeriodSeconds = 1 / rateHz;
end
assert(isscalar(framePeriodSeconds) && isfinite(framePeriodSeconds) && ...
    framePeriodSeconds > 0, 'Figure4NWB:FramePeriod', ...
    'The voltage frame period must be a positive scalar.');
assert(abs(framePeriodSeconds - 1 / rateHz) < 1e-7, 'Figure4NWB:FramePeriodMismatch', ...
    ['The spike_times resolution (%.9g s) and the correlation TimeSeries rate ', ...
     '(%.9g s) disagree by more than 100 ns.'], framePeriodSeconds, 1 / rateHz);

info = h5info(nwbFile, dataPath);
datasetSize = info.Dataspace.Size;
assert(numel(datasetSize) == 3, 'Figure4NWB:CorrelationDimensions', ...
    'Expected a 3-D correlation dataset at %s.', dataPath);
nLagSamples = datasetSize(3);
end

function values = readRowVector(nwbFile, datasetPath, expectedLength, label)
values = double(h5read(nwbFile, datasetPath));
values = values(:)';
assert(numel(values) == expectedLength, 'Figure4NWB:BehaviorLength', ...
    '%s has %d samples but the voltage timebase has %d.', ...
    label, numel(values), expectedLength);
end

function spikeFrameIndex = readSpikeFrameIndices(nwbFile, paths, nCells, ...
    voltageSeconds, framePeriodSeconds)
% /units/spike_times holds seconds, produced by convert_behavior_to_nwb.m as
% voltageSeconds(sourceFrameIndex). Invert that lookup so the returned indices
% are the 1-based voltage-frame indices the source scripts use to index traces.
spikeSeconds = double(h5read(nwbFile, paths.spikeTimes));
spikeSeconds = spikeSeconds(:);
spikeIndex = double(h5read(nwbFile, paths.spikeTimesIndex));
spikeIndex = spikeIndex(:);
assert(numel(spikeIndex) == nCells, 'Figure4NWB:SpikeIndexSize', ...
    'spike_times_index must have one entry per unit.');
assert(spikeIndex(end) == numel(spikeSeconds), 'Figure4NWB:SpikeIndexTotal', ...
    'The last spike_times_index entry must equal the total spike count.');

nFrames = numel(voltageSeconds);
frames = interp1(voltageSeconds, (1:nFrames)', spikeSeconds, 'nearest', 'extrap');
residual = abs(voltageSeconds(frames) - spikeSeconds);
assert(max(residual) < framePeriodSeconds / 2, 'Figure4NWB:SpikeFrameLookup', ...
    ['Spike times could not be matched to voltage frames: largest residual ', ...
     '%.3g s exceeds half a frame (%.3g s).'], max(residual), framePeriodSeconds / 2);

spikeFrameIndex = cell(nCells, 1);
startRow = 1;
for iCell = 1:nCells
    stopRow = spikeIndex(iCell);
    spikeFrameIndex{iCell} = frames(startRow:stopRow)';
    startRow = stopRow + 1;
end
end

function values = readCorrelationDataset(nwbFile, paths, nCells, nLagSamples, ...
    seriesName, state, triggerRange, targetRange, lagRange)
% Returns [trigger x target x lag]; see the orientation note in the help text.
datasetPath = [correlationSeriesPath(paths, seriesName) '_' ...
    correlationStateSuffix(state) '/data'];
assertDatasetExists(nwbFile, datasetPath);
info = h5info(nwbFile, datasetPath);
datasetSize = info.Dataspace.Size;
assert(isequal(datasetSize(:)', [nCells nCells nLagSamples]), ...
    'Figure4NWB:CorrelationLayout', ...
    'Expected MATLAB HDF5 layout %d x %d x %d at %s, found %s.', ...
    nCells, nCells, nLagSamples, datasetPath, mat2str(datasetSize));

triggerRange = defaultRange(triggerRange, nCells, 'triggerRange');
targetRange = defaultRange(targetRange, nCells, 'targetRange');
lagRange = defaultRange(lagRange, nLagSamples, 'lagRange');

start = [triggerRange(1) targetRange(1) lagRange(1)];
count = [diff(triggerRange) + 1, diff(targetRange) + 1, diff(lagRange) + 1];
values = double(h5read(nwbFile, datasetPath, start, count));
end

function values = readVoltageDataset(nwbFile, paths, nCells, nFrames, ...
    seriesName, cellIndex, frameRange)
datasetPath = voltageSeriesPath(paths, seriesName);
validateattributes(cellIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
cellIndex = cellIndex(:)';
assert(isequal(cellIndex, cellIndex(1):cellIndex(end)), 'Figure4NWB:CellRange', ...
    'HDF5 cell-time reads require contiguous cell indices.');
if isempty(frameRange)
    frameRange = [1 nFrames];
end
validateattributes(frameRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', nFrames});
assert(frameRange(1) <= frameRange(2), 'Figure4NWB:FrameRange', ...
    'frameRange must be ascending.');
start = [cellIndex(1) frameRange(1)];
count = [numel(cellIndex) diff(frameRange) + 1];
values = h5read(nwbFile, datasetPath, start, count);
end

function nFramesPerChunk = voltageChunkFrames(nwbFile, paths)
% Reports the HDF5 chunk extent along the frame axis, which is the natural
% read-block size for the voltage datasets. See the performance note in the
% help text: chunks span all cells, so frame blocks are the only cheap slice.
info = h5info(nwbFile, paths.subthresholdVoltage);
if isempty(info.ChunkSize)
    nFramesPerChunk = 10000;
    return
end
nFramesPerChunk = info.ChunkSize(2);
end

function datasetPath = correlationSeriesPath(paths, seriesName)
switch lower(strrep(seriesName, ' ', '_'))
    case {'voltage_cross_correlation', 'voltagecrosscorrelation', 'xcorrv', 'v'}
        datasetPath = paths.voltageCrossCorrelation;
    case {'spike_triggered_voltage', 'spiketriggeredvoltage', 'stav'}
        datasetPath = paths.spikeTriggeredVoltage;
    case {'spike_triggered_spike_rate', 'spiketriggeredspikerate', 'starate'}
        datasetPath = paths.spikeTriggeredSpikeRate;
    otherwise
        error('Figure4NWB:UnknownCorrelationSeries', ...
            'Unknown correlation series: %s', seriesName);
end
end

function suffix = correlationStateSuffix(state)
if isnumeric(state)
    assert(isscalar(state) && any(state == [1 2]), 'Figure4NWB:State', ...
        'Numeric state must be 1 (quiet) or 2 (whisking).');
    suffix = {'quiet', 'whisking'};
    suffix = suffix{state};
    return
end
switch lower(char(state))
    case {'quiet', 'q'}
        suffix = 'quiet';
    case {'whisking', 'whisk', 'w'}
        suffix = 'whisking';
    otherwise
        error('Figure4NWB:State', ...
            'State must be quiet/whisking or 1/2, not %s.', char(state));
end
end

function datasetPath = voltageSeriesPath(paths, seriesName)
switch lower(strrep(seriesName, ' ', '_'))
    case {'subthreshold', 'subthreshold_voltage', 'subthresholdvoltage', 'vsub'}
        datasetPath = paths.subthresholdVoltage;
    case {'voltage', 'normalized_voltage', 'normalizedvoltage'}
        datasetPath = paths.normalizedVoltage;
    otherwise
        error('Figure4NWB:UnknownVoltageSeries', ...
            'Unknown voltage series: %s', seriesName);
end
end

function range = defaultRange(range, extent, label)
if isempty(range)
    range = [1 extent];
end
validateattributes(range, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', extent}, mfilename, label);
assert(range(1) <= range(2), 'Figure4NWB:Range', '%s must be ascending.', label);
range = double(range(:)');
end

function assertDatasetExists(nwbFile, path)
try
    h5info(nwbFile, path);
catch ME
    error('Figure4NWB:MissingDataset', 'Required dataset is missing: %s\n%s', ...
        path, ME.message);
end
end
