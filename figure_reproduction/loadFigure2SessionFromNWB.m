function session = loadFigure2SessionFromNWB(nwbFile)
% LOADFIGURE2SESSIONFROMNWB Load Figure 2 inputs from a sparseOpto NWB file.
%
%   SESSION = LOADFIGURE2SESSIONFROMNWB(NWBFILE) reads the small, shared
%   Figure 2 inputs eagerly and exposes function handles for cell/time reads.
%   This avoids loading the full multi-gigabyte voltage and stimulation arrays.
%
%   The loader uses MATLAB's built-in HDF5 reader so it remains compatible with
%   MATLAB R2019b and does not require a separately installed MatNWB reader.
%   Dataset paths follow the authoritative MATLAB-to-NWB correspondence in
%   M_YQ0201_27_sparseOpto_NWB_progress.md.
%
%   Returned eager fields
%   ---------------------
%   cellId, cellCoordinatePixels, passedQC, properties, spikeTimes,
%   blueEvokedSpikeTriggeredAverage, voltageTimestamps,
%   framesPerTrial, numberOfTrials, framePeriodSeconds.
%
%   Returned lazy readers
%   ---------------------
%   readCellTimeSeries(seriesName, cellIndex, frameRange)
%       seriesName is 'voltage', 'subthreshold', or 'blue'. frameRange is
%       [firstFrame lastFrame], using MATLAB's one-based inclusive indexing.
%
%   stimulatedTrialIndices(cellIndex)
%       Returns trials in which this cell received blue stimulation.
%
%   readTrial(cellIndex, trialIndex)
%       Returns normalized voltage, subthreshold voltage, blue waveform, and
%       local spike-frame indices for the requested protocol trial.
%
%   trialAverageFiringRate(cellIndex)
%       Returns the per-frame probability of an inferred spike over that cell's
%       blue-stimulated trials, plus the corresponding trial indices.

if nargin < 1 || isempty(nwbFile)
    error('Figure2NWB:MissingFile', 'Provide the path to a sparseOpto NWB file.');
end
validateattributes(nwbFile, {'char'}, {'row', 'nonempty'}, mfilename, 'nwbFile');
assert(exist(nwbFile, 'file') == 2, 'Figure2NWB:FileNotFound', ...
    'NWB file not found: %s', nwbFile);

paths = figure2DatasetPaths();
required = struct2cell(paths);
for iPath = 1:numel(required)
    try
        h5info(nwbFile, required{iPath});
    catch ME
        error('Figure2NWB:MissingDataset', ...
            'Required Figure 2 dataset is missing: %s\n%s', required{iPath}, ME.message);
    end
end

session = struct();
session.nwbFile = nwbFile;
session.paths = paths;
session.cellId = h5read(nwbFile, paths.cellId);
session.cellId = session.cellId(:);
session.numberOfCells = numel(session.cellId);
assert(session.numberOfCells > 0, 'Figure2NWB:NoCells', 'The NWB units table is empty.');

session.cellCoordinatePixels = [h5read(nwbFile, paths.xPixel), h5read(nwbFile, paths.yPixel)];
session.passedQC = logical(h5read(nwbFile, paths.passedQC));
session.passedQC = session.passedQC(:);
session.properties = struct( ...
    'afterDepolarization', h5read(nwbFile, paths.afterDepolarization), ...
    'rheobase', h5read(nwbFile, paths.rheobase), ...
    'maximumFiringRateHz', h5read(nwbFile, paths.maximumFiringRateHz), ...
    'adaptation', h5read(nwbFile, paths.adaptation), ...
    'membraneTimeConstantMs', h5read(nwbFile, paths.membraneTimeConstantMs));
session.propertyMatrix = [session.properties.afterDepolarization(:), ...
    session.properties.rheobase(:), ...
    session.properties.maximumFiringRateHz(:), ...
    session.properties.adaptation(:), ...
    session.properties.membraneTimeConstantMs(:)];
session.propertyLabels = {'ADP (norm)', 'V_R_h (norm)', 'F_m_a_x (Hz)', ...
    'Adaptation', 'tau_m (ms)'};

session.voltageTimestamps = h5read(nwbFile, paths.voltageTimestamps);
session.voltageTimestamps = session.voltageTimestamps(:);
session.numberOfFrames = numel(session.voltageTimestamps);
assert(session.numberOfFrames > 1, 'Figure2NWB:InsufficientTime', ...
    'The normalized-voltage time series has fewer than two samples.');
sampleCount = min(session.numberOfFrames - 1, 100000);
session.framePeriodSeconds = median(diff(session.voltageTimestamps(1:sampleCount + 1)));
session.samplingRateHz = 1 / session.framePeriodSeconds;

trialId = h5read(nwbFile, paths.trialId);
session.numberOfTrials = numel(trialId);
assert(session.numberOfTrials > 0 && mod(session.numberOfFrames, session.numberOfTrials) == 0, ...
    'Figure2NWB:TrialFrames', ...
    'Voltage frames (%d) are not evenly divisible by protocol trials (%d).', ...
    session.numberOfFrames, session.numberOfTrials);
session.framesPerTrial = session.numberOfFrames / session.numberOfTrials;

session.blueEvokedSpikeTriggeredAverage = readCellTimeDataset( ...
    nwbFile, paths.blueEvokedSpikeTriggeredAverage, session.numberOfCells, ...
    1:session.numberOfCells, []);

session.spikeTimes = readSpikeTimes(nwbFile, paths.spikeTimes, ...
    paths.spikeTimesIndex, session.numberOfCells);

session.readCellTimeSeries = @(seriesName, cellIndex, frameRange) ...
    readNamedCellTimeSeries(nwbFile, paths, session.numberOfCells, ...
    seriesName, cellIndex, frameRange);
session.stimulatedTrialIndices = @(cellIndex) stimulatedTrialIndices( ...
    nwbFile, paths.blue, session.numberOfCells, session.framesPerTrial, cellIndex);
session.readTrial = @(cellIndex, trialIndex) readTrial( ...
    nwbFile, paths, session, cellIndex, trialIndex);
session.trialAverageFiringRate = @(cellIndex) trialAverageFiringRate(session, cellIndex);
end

function paths = figure2DatasetPaths()
paths = struct( ...
    'cellId', '/units/id', ...
    'xPixel', '/units/x_pixel', ...
    'yPixel', '/units/y_pixel', ...
    'passedQC', '/units/passed_qc', ...
    'afterDepolarization', '/units/adp_blue_norm', ...
    'rheobase', '/units/optical_rheobase_norm', ...
    'maximumFiringRateHz', '/units/max_firing_rate_hz', ...
    'adaptation', '/units/pulse_adaptation', ...
    'membraneTimeConstantMs', '/units/membrane_time_constant_ms', ...
    'voltage', '/processing/optical_voltage/normalized_voltage/data', ...
    'subthreshold', '/processing/optical_voltage/subthreshold_voltage/data', ...
    'blueEvokedSpikeTriggeredAverage', ...
        '/processing/optical_voltage/blue_evoked_spike_triggered_average/data', ...
    'blue', '/stimulus/presentation/blue_stimulation_per_cell/data', ...
    'voltageTimestamps', ...
        '/processing/optical_voltage/normalized_voltage/timestamps', ...
    'spikeTimes', '/units/spike_times', ...
    'spikeTimesIndex', '/units/spike_times_index', ...
    'trialId', '/intervals/trials/id');
end

function values = readNamedCellTimeSeries(nwbFile, paths, numberOfCells, seriesName, cellIndex, frameRange)
switch lower(seriesName)
    case 'voltage'
        datasetPath = paths.voltage;
    case 'subthreshold'
        datasetPath = paths.subthreshold;
    case 'blue'
        datasetPath = paths.blue;
    otherwise
        error('Figure2NWB:UnknownSeries', ...
            'seriesName must be ''voltage'', ''subthreshold'', or ''blue''.');
end
values = readCellTimeDataset(nwbFile, datasetPath, numberOfCells, cellIndex, frameRange);
end

function values = readCellTimeDataset(nwbFile, datasetPath, numberOfCells, cellIndex, frameRange)
validateattributes(cellIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', numberOfCells});
cellIndex = cellIndex(:)';
info = h5info(nwbFile, datasetPath);
datasetSize = info.Dataspace.Size;
assert(numel(datasetSize) == 2, 'Figure2NWB:UnexpectedDimensions', ...
    'Expected a two-dimensional cell-by-time dataset at %s.', datasetPath);
cellDimension = find(datasetSize == numberOfCells, 1, 'first');
assert(~isempty(cellDimension), 'Figure2NWB:CellDimension', ...
    'Could not identify the cell dimension at %s.', datasetPath);
timeDimension = 3 - cellDimension;
numberOfFrames = datasetSize(timeDimension);

if isempty(frameRange)
    frameRange = [1 numberOfFrames];
end
validateattributes(frameRange, {'numeric'}, {'vector', 'numel', 2, 'integer', ...
    '>=', 1, '<=', numberOfFrames});
assert(frameRange(1) <= frameRange(2), 'Figure2NWB:FrameRange', ...
    'frameRange must be ascending.');

% HDF5 cannot read arbitrary scattered indices in one call. All current Figure 2
% reads are contiguous; enforce that restriction rather than silently reading a
% different range.
assert(isequal(cellIndex, cellIndex(1):cellIndex(end)), ...
    'Figure2NWB:CellRange', 'Cell indices must be contiguous for HDF5 reads.');

start = ones(1, 2);
count = datasetSize;
start(cellDimension) = cellIndex(1);
count(cellDimension) = numel(cellIndex);
start(timeDimension) = frameRange(1);
count(timeDimension) = frameRange(2) - frameRange(1) + 1;
values = h5read(nwbFile, datasetPath, start, count);

if cellDimension == 2
    values = values';
end
end

function spikeTimes = readSpikeTimes(nwbFile, spikeTimesPath, spikeTimesIndexPath, numberOfCells)
allSpikeTimes = h5read(nwbFile, spikeTimesPath);
spikeStops = double(h5read(nwbFile, spikeTimesIndexPath));
assert(numel(spikeStops) == numberOfCells, 'Figure2NWB:SpikeIndexLength', ...
    'The spike-times index does not contain one stop index per cell.');
spikeTimes = cell(numberOfCells, 1);
first = 1;
for cellIndex = 1:numberOfCells
    last = spikeStops(cellIndex);
    spikeTimes{cellIndex} = allSpikeTimes(first:last);
    first = last + 1;
end
assert(first == numel(allSpikeTimes) + 1, 'Figure2NWB:SpikeIndexContent', ...
    'The spike-times index does not span the spike-times data.');
end

function trialIndices = stimulatedTrialIndices(nwbFile, bluePath, numberOfCells, framesPerTrial, cellIndex)
blue = readCellTimeDataset(nwbFile, bluePath, numberOfCells, cellIndex, []);
trialMatrix = reshape(blue, framesPerTrial, []);
trialIndices = find(any(trialMatrix > 0, 1));
end

function [voltage, subthreshold, blue, localSpikeFrames] = readTrial(nwbFile, paths, session, cellIndex, trialIndex)
validateattributes(cellIndex, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', session.numberOfCells});
validateattributes(trialIndex, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', session.numberOfTrials});
firstFrame = (trialIndex - 1) * session.framesPerTrial + 1;
lastFrame = trialIndex * session.framesPerTrial;
frameRange = [firstFrame lastFrame];
voltage = readCellTimeDataset(nwbFile, paths.voltage, session.numberOfCells, cellIndex, frameRange);
subthreshold = readCellTimeDataset(nwbFile, paths.subthreshold, session.numberOfCells, cellIndex, frameRange);
blue = readCellTimeDataset(nwbFile, paths.blue, session.numberOfCells, cellIndex, frameRange);
globalSpikeFrames = spikeTimesToFrames(session.spikeTimes{cellIndex}, session.voltageTimestamps);
localSpikeFrames = globalSpikeFrames(globalSpikeFrames >= firstFrame & globalSpikeFrames <= lastFrame) - firstFrame + 1;
end

function [firingRatePerFrame, trialIndices] = trialAverageFiringRate(session, cellIndex)
trialIndices = session.stimulatedTrialIndices(cellIndex);
assert(~isempty(trialIndices), 'Figure2NWB:NoStimulatedTrials', ...
    'Cell %d has no blue-stimulated trials.', cellIndex);
globalSpikeFrames = spikeTimesToFrames(session.spikeTimes{cellIndex}, session.voltageTimestamps);
spikeTrial = floor((globalSpikeFrames - 1) / session.framesPerTrial) + 1;
spikeFrameWithinTrial = mod(globalSpikeFrames - 1, session.framesPerTrial) + 1;
keep = ismember(spikeTrial, trialIndices);
spikeTrial = spikeTrial(keep);
spikeFrameWithinTrial = spikeFrameWithinTrial(keep);

% The legacy implementation stores spike occurrences in a logical matrix.
% Matching that representation prevents repeated frame indices from changing the
% trial average if they are present in a source spike-time vector.
uniqueTrialFrames = unique([spikeTrial(:), spikeFrameWithinTrial(:)], 'rows');
firingRatePerFrame = accumarray(uniqueTrialFrames(:, 2), 1, ...
    [session.framesPerTrial 1], @sum, 0) ./ numel(trialIndices);
end

function frameIndex = spikeTimesToFrames(spikeTimes, voltageTimestamps)
if isempty(spikeTimes)
    frameIndex = zeros(0, 1);
    return
end
frameCoordinates = (1:numel(voltageTimestamps))';
frameIndex = interp1(voltageTimestamps, frameCoordinates, spikeTimes(:), 'nearest', 'extrap');
frameIndex = round(frameIndex);
frameIndex = frameIndex(frameIndex >= 1 & frameIndex <= numel(voltageTimestamps));
end
