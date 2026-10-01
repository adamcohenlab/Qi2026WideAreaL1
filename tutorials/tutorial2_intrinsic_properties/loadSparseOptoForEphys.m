function data = loadSparseOptoForEphys(nwbFile, varargin)
%LOADSPARSEOPTOFOREPHYS Load a sparseOpto NWB file for Tutorial 2.
%
%   data = LOADSPARSEOPTOFOREPHYS(nwbFile) returns a struct whose fields carry
%   the *source pipeline's* variable names, so the tutorial bodies read like
%   sparseOpto_ephysProp_meas_cleanup.m.
%
%   data = LOADSPARSEOPTOFOREPHYS(nwbFile, 'Name', value, ...) accepts:
%
%     'ScanBlue'   true (default) makes one streaming pass over blue_all to
%                  build the per-cell stimulated-chunk map and to extract the
%                  single shared stimulus waveform. About 30 to 75 s. Set false
%                  for metadata only; nothing downstream will work.
%     'Verbose'    true (default) prints a load report.
%
%   FIELDS
%     nCells, T, nFrame, nChunk, dt, dx
%     cell_coord   nCells x 2, [x_pixel y_pixel], source column order.
%     gfpIdx       nCells x 1 logical, GFP/NPY label.
%     idxUse       nCells x 1 logical, the source QC mask.
%     spk_t        nCells x 1 cell array of 1-based frame indices. This is the
%                  deposited spk_t, which is the source spk_t_all: NOT motion
%                  corrected. See SPIKE TIMES below.
%     stimChunk    nCells x nChunk logical, true where that cell was stimulated.
%     blue         nFrame x 1 single, the stimulus waveform. One waveform is
%                  shared by every cell and every stimulated chunk; the load
%                  pass asserts this rather than assuming it.
%     readChunk    readChunk(series, chunkIndex) -> nCells x nFrame single.
%                  series is 'voltage', 'subthreshold' or 'blue'.
%     deposited    struct of the stored per-cell properties, used as the
%                  comparison target in each stage:
%                  adp_blue, vRheobase, frMax, frPulseAdp, membraneC, and
%                  selfBlueStaFN_ADP (nCells x 281).
%
%   NORMALISATION
%   The deposit stores spike-height-normalised voltage. From dataPrepsForMAT.m:
%
%       traces_all_n     = traces_all_bc        ./ spkHgtNoBlue
%       traces_all_sub_n = traces_all_sub_sgfilt./ spkHgtNoBlue
%
%   Every place the source script divides by spkHgtNoBlue is a linear operation
%   on the trace (a spike-triggered average, or a trial average), so it commutes
%   with the division and the tutorial simply computes on the normalised trace.
%   spkHgtNoBlue is never needed and is not in the deposit. The source guard
%   spkHgtNoBlue(spkHgtNoBlue==0) = nan never fires on either session: no cell
%   has zero spike height, and there is no non-finite sample in either asset.
%
%   SPIKE TIMES
%   The source script motion-corrects the spike set before deriving blue-evoked
%   spikes, via spike_times_motion_correction_aggressive(spk_t, spk_t_mat,
%   mcTraceAvgCor', 4). That needs mcTrace_all, which the sparseOpto deposit
%   does not carry, so the tutorial works from the uncorrected spike set. This
%   is the one irreducible departure from the published numbers; it is measured
%   per property in TUTORIAL2_VALIDATION.md.
%
%   The source also strips a chunk-boundary window t_noise = (0:nFrame:T) +
%   (-140:140) from spk_t_mat_c. That window is frames 1:140 and 7861:8000 of
%   each chunk, and the stimulus occupies 393:7480, so t_noise cannot intersect
%   any window this tutorial measures. It is therefore not implemented, and
%   tutorial2a asserts the disjointness rather than leaving it to the reader.
%
%   SUBTHRESHOLD
%   data.readChunk('subthreshold', k) returns the DEPOSITED traces_all_sub_n,
%   which was built with a median-filter substitution and a Savitzky-Golay pass.
%   That is not the subthreshold this analysis uses. Rheobase wants the plain
%   linear interpolation of sparseOpto_ephysProp_meas_cleanup.m, because the
%   trial average that follows already removes the noise the filter targets and
%   the filter would distort the slow ramp the cubic fit reads. Build it with
%   subthresholdByInterpolation; the deposited series is a cross-check only.
%
%   See also TUTORIAL2A_STIMULUSPROTOCOL, SUBTHRESHOLDBYINTERPOLATION.

%% ---------------------------------------------------------------- options
opt = struct('ScanBlue', true, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
assert(exist(nwbFile, 'file') == 2, 'Tutorial2:MissingFile', ...
    'NWB file not found: %s', nwbFile);

timer = tic;

paths = struct( ...
    'voltage',      '/processing/optical_voltage/normalized_voltage/data', ...
    'subthreshold', '/processing/optical_voltage/subthreshold_voltage/data', ...
    'blue',         '/stimulus/presentation/blue_stimulation_per_cell/data');

%% ------------------------------------------------------------ dimensions
info = h5info(nwbFile, paths.voltage);
sz = info.Dataspace.Size;          % MATLAB order: [nCells nFrames]

data = struct();
data.nwbFile = nwbFile;
data.paths = paths;
data.nCells = sz(1);
data.T = sz(2);
data.dt = 1.27e-3;                 % confirmed voltage-camera exposure
data.dx = 6.5;                     % um per pixel

% Chunk structure comes from the deposit, not from a literal 8000, so the
% tutorial is not tied to one session's acquisition length.
chunkStart = h5read(nwbFile, '/intervals/recording_chunks/start_time');
data.nChunk = numel(chunkStart);
assert(mod(data.T, data.nChunk) == 0, 'Tutorial2:ChunkStructure', ...
    'Frame count %d is not a whole multiple of the %d recording chunks.', ...
    data.T, data.nChunk);
data.nFrame = data.T / data.nChunk;

%% ------------------------------------------------------------ unit table
xPixel = h5read(nwbFile, '/units/x_pixel');
yPixel = h5read(nwbFile, '/units/y_pixel');
data.cell_coord = [xPixel(:), yPixel(:)];
data.gfpIdx = logical(h5read(nwbFile, '/units/gfp_positive'));
data.idxUse = logical(h5read(nwbFile, '/units/passed_qc'));

%% ---------------------------------------------------------- spike frames
% Spike times are stored in session-relative seconds against a timestamp vector
% that preserves the gaps between acquisition chunks, so the seconds-to-frame
% map is a nearest-neighbour lookup, not a division by dt.
voltageTimestamps = h5read(nwbFile, ...
    '/processing/optical_voltage/normalized_voltage/timestamps');
spikeTimes = h5read(nwbFile, '/units/spike_times');
spikeIndex = h5read(nwbFile, '/units/spike_times_index');

frameCoordinates = (1:numel(voltageTimestamps))';
allFrames = interp1(voltageTimestamps, frameCoordinates, spikeTimes(:), ...
    'nearest', 'extrap');
assert(all(isfinite(allFrames)), 'Tutorial2:SpikeFrameMap', ...
    'Spike-time to frame mapping produced non-finite indices.');

data.spk_t = cell(data.nCells, 1);
first = 1;
for ii = 1:data.nCells
    last = double(spikeIndex(ii));
    data.spk_t{ii} = reshape(allFrames(first:last), 1, []);
    first = last + 1;
end
assert(first == numel(spikeTimes) + 1, 'Tutorial2:SpikeIndex', ...
    'units/spike_times_index does not span units/spike_times.');

%% --------------------------------------------------- deposited products
dep = struct();
dep.adp_blue   = double(h5read(nwbFile, '/units/adp_blue_norm'));
dep.vRheobase  = double(h5read(nwbFile, '/units/optical_rheobase_norm'));
dep.frMax      = double(h5read(nwbFile, '/units/max_firing_rate_hz'));
dep.frPulseAdp = double(h5read(nwbFile, '/units/pulse_adaptation'));
dep.membraneC  = double(h5read(nwbFile, '/units/membrane_time_constant_ms'));
dep.selfBlueStaFN_ADP = double(h5read(nwbFile, ...
    '/processing/optical_voltage/blue_evoked_spike_triggered_average/data'));
% The series is written time-by-cell in HDF5 and comes back cell-by-time in
% MATLAB order, which is the source orientation. nBack = nFront = 140.
assert(size(dep.selfBlueStaFN_ADP, 1) == data.nCells, 'Tutorial2:StaOrientation', ...
    'Deposited STA is %s; expected %d rows.', ...
    mat2str(size(dep.selfBlueStaFN_ADP)), data.nCells);
dep.nBack = (size(dep.selfBlueStaFN_ADP, 2) - 1) / 2;
dep.nFront = dep.nBack;
data.deposited = dep;

%% -------------------------------------------------------------- readers
data.readChunk = @(series, k) readChunkImpl(nwbFile, paths, data.nCells, ...
    data.nFrame, data.nChunk, series, k);

%% ------------------------------------------------- blue scan, stim map
if opt.ScanBlue
    [data.stimChunk, data.blue, blueDeviation] = ...
        scanBlue(nwbFile, paths.blue, data.nCells, data.nFrame, data.nChunk);
    data.blueDeviation = blueDeviation;
else
    data.stimChunk = [];
    data.blue = [];
    data.blueDeviation = NaN;
end

%% ---------------------------------------------------------------- report
if opt.Verbose
    nSpk = cellfun(@numel, data.spk_t);
    fprintf('Loaded %s\n', nwbFile);
    fprintf('  %d cells, %d frames, %d chunks of %d frames, dt = %.2f ms\n', ...
        data.nCells, data.T, data.nChunk, data.nFrame, data.dt*1e3);
    fprintf('  spikes: %d total, median %g/cell [%d, %d]\n', ...
        sum(nSpk), median(nSpk), min(nSpk), max(nSpk));
    if opt.ScanBlue
        nStim = sum(data.stimChunk, 2);
        fprintf('  stimulated chunks/cell: median %d [%d, %d]\n', ...
            median(nStim), min(nStim), max(nStim));
        fprintf('  shared stimulus waveform, max deviation %.3g\n', blueDeviation);
    else
        fprintf('  blue NOT scanned (ScanBlue false)\n');
    end
    fprintf('  %.1f s\n', toc(timer));
end
end

% -------------------------------------------------------------------------
function block = readChunkImpl(nwbFile, paths, nCells, nFrame, nChunk, series, k)
%READCHUNKIMPL Read one acquisition chunk, all cells.
%
% Every series is chunked [nCells 10000] in MATLAB order, so a whole-population
% time block is the cheap read and a single cell's full trace is the expensive
% one. See the I/O note in README.md before changing this to a per-cell loop.
name = validatestring(series, {'voltage', 'subthreshold', 'blue'});
validateattributes(k, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', nChunk});
block = h5read(nwbFile, paths.(name), [1 (k-1)*nFrame+1], [nCells nFrame]);
end

% -------------------------------------------------------------------------
function [stimChunk, blue, deviation] = scanBlue(nwbFile, bluePath, nCells, nFrame, nChunk)
%SCANBLUE One pass over blue_all: which cells are stimulated in which chunk,
% and the stimulus waveform itself.
%
% The protocol turns out to be a single waveform delivered to a rotating subset
% of cells, so this pass replaces every later read of blue_all. The equality is
% asserted here, at a cost of one max() per stimulated cell, rather than assumed
% at each use.
stimChunk = false(nCells, nChunk);
blue = [];
deviation = 0;
for k = 1:nChunk
    block = h5read(nwbFile, bluePath, [1 (k-1)*nFrame+1], [nCells nFrame]);
    stimulated = max(block, [], 2) > 0;
    stimChunk(:, k) = stimulated;
    if ~any(stimulated); continue; end
    rows = block(stimulated, :);
    if isempty(blue)
        blue = reshape(rows(1, :), [], 1);
    end
    deviation = max(deviation, max(max(abs(rows - blue'))));
end
assert(~isempty(blue), 'Tutorial2:NoStimulation', ...
    'No chunk in this asset carries blue stimulation.');
assert(deviation < 1e-5, 'Tutorial2:BlueNotShared', ...
    ['Blue stimulation is not a single shared waveform (max deviation %.3g). ' ...
     'Tutorial 2 assumes one protocol; rewrite the stages to read blue per ' ...
     'cell before trusting any result.'], deviation);
end
