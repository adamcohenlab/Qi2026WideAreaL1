function session = loadFigure3SessionFromNWB(nwbFile)
% LOADFIGURE3SESSIONFROMNWB Load Figure 3 connectivity inputs from NWB.
%
%   SESSION = LOADFIGURE3SESSIONFROMNWB(NWBFILE) reads the static pairwise
%   connectivity metrics eagerly and provides lazy readers for the 281-sample
%   pairwise triggered-average arrays. Source MATLAB array orientation is
%   preserved on output: presynaptic/prejunctional cell x postsynaptic cell x
%   time offset.
%
%   The loader uses built-in HDF5 functions, so it is compatible with MATLAB
%   R2019b and does not require MatNWB at reproduction time. Dataset paths and
%   orientation follow the embedded MATLAB-to-NWB lookup table written by
%   convert_sparsePulseHad_to_nwb.m.

if nargin < 1 || isempty(nwbFile)
    error('Figure3NWB:MissingFile', 'Provide a sparsePulseHad NWB file.');
end
validateattributes(nwbFile, {'char'}, {'row', 'nonempty'}, mfilename, 'nwbFile');
assert(exist(nwbFile, 'file') == 2, 'Figure3NWB:FileNotFound', ...
    'NWB file not found: %s', nwbFile);

paths = figure3DatasetPaths();
requiredPaths = {paths.cellId, paths.xPixel, paths.yPixel, paths.passedQC, ...
    paths.estimatedSpikesPerPulse, paths.metricsId, paths.preCell, paths.postCell, ...
    paths.synapticConnection, paths.ipspAmplitude, paths.ipspDecay, paths.ipspDecaySE, ...
    paths.distancePixels, paths.gapJunctionConnection, paths.spikeletAmplitude, ...
    paths.bluePulseIPSP, paths.blueLightCrosstalk, paths.hadamardCrossSpikeSTA};
for iPath = 1:numel(requiredPaths)
    assertDatasetExists(nwbFile, requiredPaths{iPath});
end

session = struct();
session.nwbFile = nwbFile;
session.paths = paths;
session.cellId = h5read(nwbFile, paths.cellId);
session.cellId = double(session.cellId(:));
session.numberOfCells = numel(session.cellId);
assert(session.numberOfCells > 0, 'Figure3NWB:NoCells', 'The units table is empty.');

session.cellCoordinatePixels = [double(h5read(nwbFile, paths.xPixel)), ...
    double(h5read(nwbFile, paths.yPixel))];
session.passedQC = logical(h5read(nwbFile, paths.passedQC));
session.passedQC = session.passedQC(:);
session.estimatedSpikesPerBluePulse = double(h5read(nwbFile, paths.estimatedSpikesPerPulse));
session.estimatedSpikesPerBluePulse = session.estimatedSpikesPerBluePulse(:);
assert(numel(session.passedQC) == session.numberOfCells, ...
    'Figure3NWB:QCSize', 'passed_qc must have one value per unit.');

% dx is retained in the NWB lookup-table notes and in x_um/y_um. The Figure 3
% source uses this recorded calibration to turn distance_pixels into micrometers.
session.pixelSizeMicrometers = 6.5;
session.figureCellIndices = find(session.passedQC);
if isempty(session.figureCellIndices)
    error('Figure3NWB:NoQCCells', 'No units are marked passed_qc in %s.', nwbFile);
end

session.metrics = readPairwiseMetrics(nwbFile, paths, session.numberOfCells);
[session.pairwiseTimeSeconds, session.pairwiseFramePeriodSeconds] = ...
    readPairwiseTimebase(nwbFile, paths.bluePulseIPSP);

session.readPairwiseSeries = @(seriesName, preIndex, postIndex, timeRange) ...
    readNamedPairwiseSeries(nwbFile, paths, session.numberOfCells, ...
    seriesName, preIndex, postIndex, timeRange);
session.readCellTimeSeries = @(seriesName, cellIndex, frameRange) ...
    readNamedCellTimeSeries(nwbFile, paths, session.numberOfCells, ...
    seriesName, cellIndex, frameRange);
session.sortByX = @(cellIndex) sortCellsByX(session.cellCoordinatePixels, cellIndex);
end

function paths = figure3DatasetPaths()
paths = struct( ...
    'cellId', '/units/id', ...
    'xPixel', '/units/x_pixel', ...
    'yPixel', '/units/y_pixel', ...
    'passedQC', '/units/passed_qc', ...
    'estimatedSpikesPerPulse', '/units/estimated_spikes_per_blue_pulse', ...
    'metricsId', '/analysis/pairwise_connectivity_metrics/id', ...
    'preCell', '/analysis/pairwise_connectivity_metrics/presynaptic_unit_id', ...
    'postCell', '/analysis/pairwise_connectivity_metrics/postsynaptic_unit_id', ...
    'synapticConnection', '/analysis/pairwise_connectivity_metrics/synaptic_connection', ...
    'ipspAmplitude', '/analysis/pairwise_connectivity_metrics/ipsp_amplitude_norm_per_presynaptic_spike', ...
    'ipspDecay', '/analysis/pairwise_connectivity_metrics/ipsp_decay_source_units', ...
    'ipspDecaySE', '/analysis/pairwise_connectivity_metrics/ipsp_decay_standard_error_source_units', ...
    'distancePixels', '/analysis/pairwise_connectivity_metrics/distance_pixels', ...
    'gapJunctionConnection', '/analysis/pairwise_connectivity_metrics/gap_junction_connection', ...
    'spikeletAmplitude', '/analysis/pairwise_connectivity_metrics/spikelet_amplitude_normalized', ...
    'bluePulseIPSP', '/processing/connectivity_mapping/blue_pulse_ipsp_waveform/data', ...
    'blueLightCrosstalk', '/processing/connectivity_mapping/blue_light_crosstalk_estimate/data', ...
    'hadamardCrossSpikeSTA', '/processing/connectivity_mapping/hadamard_cross_spike_sta/data', ...
    'normalizedVoltage', '/processing/optical_voltage/normalized_voltage/data', ...
    'blueStimulation', '/stimulus/presentation/blue_stimulation_per_cell/data');
end

function metrics = readPairwiseMetrics(nwbFile, paths, nCells)
pre = double(h5read(nwbFile, paths.preCell));
post = double(h5read(nwbFile, paths.postCell));
nPairs = nCells * nCells;
assert(numel(h5read(nwbFile, paths.metricsId)) == nPairs && ...
    numel(pre) == nPairs && numel(post) == nPairs, ...
    'Figure3NWB:MetricRows', 'The pairwise metrics table must have nCells^2 rows.');
assert(all(pre >= 0 & pre < nCells & post >= 0 & post < nCells), ...
    'Figure3NWB:MetricIndex', 'Pairwise unit identifiers are outside the units table.');
tableIndex = sub2ind([nCells nCells], pre + 1, post + 1);
assert(numel(unique(tableIndex)) == nPairs, 'Figure3NWB:DuplicatePair', ...
    'The pairwise metrics table does not provide exactly one row per directed pair.');

% The static metric columns retain the source MATLAB matrix serialization
% (metric(:): presynaptic row varies fastest). The current NWB table's explicit
% pre/post ID columns vary in the opposite order, so they cannot safely be used
% to remap those columns. Restore the documented source order directly. This is
% also robust to a future NWB rewrite that corrects the table ID ordering.
sourceLinearIndex = (1:nPairs)';

metrics = struct();
metrics.synapticConnection = pairwiseColumnToMatrix(nwbFile, paths.synapticConnection, sourceLinearIndex, nCells, true);
metrics.ipspAmplitude = pairwiseColumnToMatrix(nwbFile, paths.ipspAmplitude, sourceLinearIndex, nCells, false);
metrics.ipspDecay = pairwiseColumnToMatrix(nwbFile, paths.ipspDecay, sourceLinearIndex, nCells, false);
metrics.ipspDecaySE = pairwiseColumnToMatrix(nwbFile, paths.ipspDecaySE, sourceLinearIndex, nCells, false);
metrics.distancePixels = pairwiseColumnToMatrix(nwbFile, paths.distancePixels, sourceLinearIndex, nCells, false);
metrics.gapJunctionConnection = pairwiseColumnToMatrix(nwbFile, paths.gapJunctionConnection, sourceLinearIndex, nCells, true);
metrics.spikeletAmplitude = pairwiseColumnToMatrix(nwbFile, paths.spikeletAmplitude, sourceLinearIndex, nCells, false);
end

function matrix = pairwiseColumnToMatrix(nwbFile, path, linearIndex, nCells, asLogical)
values = h5read(nwbFile, path);
values = values(:);
assert(numel(values) == nCells * nCells, 'Figure3NWB:MetricColumn', ...
    'Metric column at %s has an unexpected length.', path);
matrix = nan(nCells);
matrix(linearIndex) = double(values);
if asLogical
    matrix = logical(matrix);
end
end

function [timeSeconds, framePeriodSeconds] = readPairwiseTimebase(nwbFile, dataPath)
seriesPath = fileparts(dataPath);
startingTimePath = [seriesPath '/starting_time'];
assertDatasetExists(nwbFile, startingTimePath);
firstOffset = double(h5read(nwbFile, startingTimePath));
try
    rateHz = double(h5readatt(nwbFile, startingTimePath, 'rate'));
catch
    error('Figure3NWB:MissingPairRate', ...
        'Pairwise TimeSeries starting_time lacks its required rate attribute: %s', startingTimePath);
end
assert(isscalar(rateHz) && isfinite(rateHz) && rateHz > 0, ...
    'Figure3NWB:PairRate', 'Pairwise TimeSeries rate must be positive.');
framePeriodSeconds = 1 / rateHz;
info = h5info(nwbFile, dataPath);
datasetSize = info.Dataspace.Size;
assert(numel(datasetSize) == 3, 'Figure3NWB:PairDimensions', ...
    'Expected a 3-D pairwise dataset at %s.', dataPath);
nTime = datasetSize(3);
timeSeconds = firstOffset + (0:nTime-1)' * framePeriodSeconds;
end

function values = readNamedPairwiseSeries(nwbFile, paths, nCells, seriesName, preIndex, postIndex, timeRange)
switch lower(seriesName)
    case {'bluepulseipsp', 'blue_pulse_ipsp'}
        datasetPath = paths.bluePulseIPSP;
    case {'bluelightcrosstalk', 'blue_light_crosstalk'}
        datasetPath = paths.blueLightCrosstalk;
    case {'hadamardcrossspikesta', 'hadamard_cross_spike_sta'}
        datasetPath = paths.hadamardCrossSpikeSTA;
    otherwise
        error('Figure3NWB:UnknownPairwiseSeries', ...
            'Unknown pairwise series: %s', seriesName);
end
values = readPairwiseDataset(nwbFile, datasetPath, nCells, preIndex, postIndex, timeRange);
end

function values = readPairwiseDataset(nwbFile, datasetPath, nCells, preIndex, postIndex, timeRange)
validateattributes(preIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
validateattributes(postIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
preIndex = preIndex(:)';
postIndex = postIndex(:)';
assert(isequal(preIndex, preIndex(1):preIndex(end)) && ...
    isequal(postIndex, postIndex(1):postIndex(end)), 'Figure3NWB:PairRange', ...
    'HDF5 pairwise reads require contiguous pre- and postsynaptic index ranges.');
info = h5info(nwbFile, datasetPath);
datasetSize = info.Dataspace.Size;
assert(numel(datasetSize) == 3 && datasetSize(1) == nCells && datasetSize(2) == nCells, ...
    'Figure3NWB:PairDatasetLayout', ...
    'Expected MATLAB HDF5 layout pre x post x time at %s.', datasetPath);
nTime = datasetSize(3);
if isempty(timeRange)
    timeRange = [1 nTime];
end
validateattributes(timeRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', nTime});
assert(timeRange(1) <= timeRange(2), 'Figure3NWB:PairTimeRange', ...
    'timeRange must be ascending.');

% MATLAB's HDF5 interface restores the source MATLAB dimension order even
% though generic HDF5 readers report the physical on-disk order as time x post
% x pre. The returned values therefore already have pre x post x time order.
start = [preIndex(1) postIndex(1) timeRange(1)];
count = [numel(preIndex) numel(postIndex) timeRange(2) - timeRange(1) + 1];
values = h5read(nwbFile, datasetPath, start, count);
end

function values = readNamedCellTimeSeries(nwbFile, paths, nCells, seriesName, cellIndex, frameRange)
switch lower(seriesName)
    case {'voltage', 'normalizedvoltage', 'normalized_voltage'}
        datasetPath = paths.normalizedVoltage;
    case {'blue', 'bluestimulation', 'blue_stimulation'}
        datasetPath = paths.blueStimulation;
    otherwise
        error('Figure3NWB:UnknownCellSeries', 'Unknown cell time series: %s', seriesName);
end
validateattributes(cellIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nCells});
cellIndex = cellIndex(:)';
assert(isequal(cellIndex, cellIndex(1):cellIndex(end)), 'Figure3NWB:CellRange', ...
    'HDF5 cell-time reads require contiguous cell indices.');
info = h5info(nwbFile, datasetPath);
datasetSize = info.Dataspace.Size;
assert(numel(datasetSize) == 2 && datasetSize(1) == nCells, ...
    'Figure3NWB:CellDatasetLayout', 'Expected MATLAB HDF5 layout cell x time at %s.', datasetPath);
nFrames = datasetSize(2);
if isempty(frameRange)
    frameRange = [1 nFrames];
end
validateattributes(frameRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', nFrames});
start = [cellIndex(1) frameRange(1)];
count = [numel(cellIndex) frameRange(2) - frameRange(1) + 1];
values = h5read(nwbFile, datasetPath, start, count);
end

function sortedIndex = sortCellsByX(coordinates, cellIndex)
if nargin < 2 || isempty(cellIndex)
    cellIndex = 1:size(coordinates, 1);
end
[~, localOrder] = sort(coordinates(cellIndex, 1));
sortedIndex = cellIndex(localOrder);
end

function assertDatasetExists(nwbFile, path)
try
    h5info(nwbFile, path);
catch ME
    error('Figure3NWB:MissingDataset', 'Required dataset is missing: %s\n%s', path, ME.message);
end
end
