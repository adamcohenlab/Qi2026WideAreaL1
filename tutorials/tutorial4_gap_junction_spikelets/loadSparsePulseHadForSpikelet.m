function data = loadSparsePulseHadForSpikelet(nwbFile, varargin)
%LOADSPARSEPULSEHADFORSPIKELET Load a sparsePulseHad NWB file for Tutorial 4.
%
%   data = LOADSPARSEPULSEHADFORSPIKELET(nwbFile) returns a struct whose
%   fields carry the *source pipeline's* variable names, so the tutorial
%   bodies read like `sparsePulseHad_spikelet_control.m` and can call the
%   original helper functions (get_sta_mat_single, get_sta_mat_self)
%   unchanged.
%
%   Options:
%     'LoadTraces'   true (default) loads the nCells-by-nFrames normalized
%                    voltage, about 4 GB as single. false gives a
%                    metadata-and-spikes load.
%     'LoadSpikes'   true (default) reads /units/spike_times and converts it
%                    back to 1-based voltage-frame indices.
%     'LoadDepositedTriggers'  false (default). true also reads the ragged
%                    hadamard_orthogonal_spike_times column, a 39 M element
%                    read that only stage 4a needs.
%     'Verbose'      true (default).
%
%   FIELDS
%     traces_all     nCells x nFrames single, the deposited normalized
%                    voltage. See NORMALISATION.
%     blue_all       nCells x nFrames single, source scale retained.
%     spk_t          nCells x 1 cell of 1-based voltage-frame indices.
%     cell_coord     nCells x 2, [x_pixel y_pixel].
%     distMat        nCells x nCells, pixels, read from the deposit.
%     dt, dx, nCells, nFrame, nChunk, T
%     deposited      struct of stored products used as ground truth:
%                    spkletHadAmpCor, pMat_spklet, readHadSta(pre,post),
%                    and optionally spk_t_had_ortho.
%
%   NORMALISATION
%   The deposit stores normalized voltage. From dataPrepsForMAT.m:
%
%       traces_all_n    = traces_all_bc ./ spkHgtNoBlueAvg
%       spkletHadAmpCor = spkletH ./ spkHgtNoBlueAvg'
%
%   and `sparsePulseHad_spikelet_cm_correction.m` divides by the same vector
%   before correcting and multiplies it back afterwards, so that round trip
%   cancels. The whole chain therefore carries exactly one division by the
%   post cell's spike height, and computing from the deposited normalized
%   trace applies it once. Every spikelet quantity below is in normalized
%   units and `spkHgtNoBlueAvg` is never needed, which is just as well
%   because it is not in the deposit.
%
%   There is no nSpikesPulse factor anywhere on this path, unlike the IPSP
%   path in Tutorial 3. The spikelet null is trigger-count matched by
%   construction rather than rate normalised, so observed and null share the
%   post cell's scale exactly and it cancels out of every p-value.
%
%   SPIKE TIMES ARE THE UNCORRECTED SET
%   `sparsePulseHad_spikelet_control.m` runs
%   `spike_times_motion_correction_aggressive` before building the Hadamard
%   spike matrix, and its input `mcTrace_all` is not in the deposit.
%   /units/spike_times is the uncorrected `spk_t_all`. Measured on
%   M-YQ0201-27, this costs essentially nothing inside the Hadamard epochs:
%   the correction only ever removes events, never shifts them, and
%   correction plus masking drops a median 0.01 per cent of them, because the
%   +/-140 frame guard around chunk boundaries already excludes the frames
%   where motion-flagged spikes live. TUTORIAL4_VALIDATION.md has the
%   per-pair trigger-set comparison against the deposit.
%
%   See also TUTORIAL4A_ORTHOGONALTRIGGERS, BUILDORTHOGONALTRIGGERS.

%% ---------------------------------------------------------------- options
opt = struct('LoadTraces', true, 'LoadSpikes', true, ...
    'LoadDepositedTriggers', false, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
assert(exist(nwbFile, 'file') == 2, 'Tutorial4:MissingFile', ...
    'NWB file not found: %s', nwbFile);

timer = tic;

%% ------------------------------------------------------------ dimensions
info = h5info(nwbFile, '/processing/optical_voltage/normalized_voltage/data');
sz = info.Dataspace.Size;              % MATLAB order: [nCells nFrames]
nCells = sz(1);
T = sz(2);

data = struct();
data.nwbFile = nwbFile;
data.nCells = nCells;
data.T = T;
data.dt = 1.27e-3;                     % confirmed voltage-camera exposure
data.dx = 6.5;                         % um per pixel

chunkStart = h5read(nwbFile, '/intervals/recording_chunks/start_time');
data.nChunk = numel(chunkStart);
assert(mod(T, data.nChunk) == 0, 'Tutorial4:ChunkStructure', ...
    'Frame count %d is not a whole multiple of the %d recording chunks.', T, data.nChunk);
data.nFrame = T / data.nChunk;

%% ------------------------------------------------------------ unit table
xPixel = h5read(nwbFile, '/units/x_pixel');
yPixel = h5read(nwbFile, '/units/y_pixel');
data.cell_coord = [xPixel(:), yPixel(:)];
data.gfpIdx = logical(h5read(nwbFile, '/units/gfp_positive'));
data.idxUse = logical(h5read(nwbFile, '/units/passed_qc'));

data.distMat = reshape(double(h5read(nwbFile, ...
    '/analysis/pairwise_connectivity_metrics/distance_pixels')), nCells, nCells);
distCheck = squareform(pdist(data.cell_coord));
relErr = max(abs(data.distMat(:) - distCheck(:))) / max(distCheck(:));
assert(relErr < 1e-5, 'Tutorial4:DistanceMismatch', ...
    'Deposited distMat disagrees with pdist(cell_coord); relative error %.2e.', relErr);

%% ----------------------------------------------------- blue stimulation
data.blue_all = h5read(nwbFile, '/stimulus/presentation/blue_stimulation_per_cell/data');

%% ----------------------------------------------------------- spike times
data.tVolt = h5read(nwbFile, '/processing/optical_voltage/normalized_voltage/timestamps');
if opt.LoadSpikes
    data.spk_t = readRaggedFrames(nwbFile, '/units/spike_times', data.tVolt);
    assert(numel(data.spk_t) == nCells, 'Tutorial4:SpikeCount', ...
        'Expected %d spike-time rows, found %d.', nCells, numel(data.spk_t));
else
    data.spk_t = {};
end

%% --------------------------------------------------- deposited products
base = '/analysis/pairwise_connectivity_metrics/';
readPair = @(name) reshape(double(h5read(nwbFile, [base name])), nCells, nCells);
dep = struct();
dep.spkletHadAmpCor = readPair('spikelet_amplitude_normalized');
dep.pMat_spklet = logical(readPair('gap_junction_connection'));

staPath = '/processing/connectivity_mapping/hadamard_cross_spike_sta/data';
staInfo = h5info(nwbFile, staPath);
nTau = staInfo.Dataspace.Size(3);
dep.nTau = nTau;
dep.tBack = floor((nTau-1)/2);
dep.tFront = nTau - dep.tBack - 1;
dep.staPath = staPath;
% One pair, and one whole presynaptic row. The tensor is chunked
% [nCells nCells 25], so a single-pair read still touches whole chunks; the
% row read costs the same and is what stage 4a actually wants.
dep.readHadSta = @(pre,post) squeeze(h5read(nwbFile, staPath, [pre post 1], [1 1 nTau]));
dep.readHadStaRow = @(pre) squeeze(h5read(nwbFile, staPath, [pre 1 1], [1 nCells nTau]));

% Coverage. Patched into the dataset description on 2026-09-21: only pairs
% closer than rLim = 400 um were ever computed, everything else is exact
% zeros, and the diagonal is a self average rather than a cross average.
dep.rLimDeposit = 400;
dep.computedMask = data.distMat < dep.rLimDeposit / data.dx;

if opt.LoadDepositedTriggers
    dep.spk_t_had_ortho = reshape( ...
        readRaggedFrames(nwbFile, [base 'hadamard_orthogonal_spike_times'], data.tVolt), ...
        nCells, nCells);           % source order: presynaptic index varies fastest
end
data.deposited = dep;

%% ---------------------------------------------------------------- traces
if opt.LoadTraces
    data.traces_all = h5read(nwbFile, '/processing/optical_voltage/normalized_voltage/data');
else
    data.traces_all = [];
end

%% ---------------------------------------------------------------- report
if opt.Verbose
    fprintf('Loaded %s\n', nwbFile);
    fprintf('  %d cells, %d frames, %d chunks of %d frames, dt = %.2f ms\n', ...
        nCells, T, data.nChunk, data.nFrame, data.dt*1e3);
    if opt.LoadSpikes
        n = cellfun(@numel, data.spk_t);
        fprintf('  spikes: %d total, median %d/cell [%d, %d]\n', ...
            sum(n), median(n), min(n), max(n));
    end
    fprintf('  deposited pMat_spklet: %d positive directed pairs\n', sum(dep.pMat_spklet(:)));
    fprintf('  deposited STA covers %d of %d pairs (< %d um)\n', ...
        sum(dep.computedMask(:)), nCells^2, dep.rLimDeposit);
    if opt.LoadTraces
        w = whos('data');
        fprintf('  traces loaded, struct is %.2f GB\n', w.bytes/2^30);
    else
        fprintf('  traces NOT loaded (LoadTraces false)\n');
    end
    fprintf('  %.1f s\n', toc(timer));
end
end

%% ------------------------------------------------------------- helpers
function frames = readRaggedFrames(nwbFile, dataPath, tVolt)
%READRAGGEDFRAMES Read a ragged seconds column back to 1-based frame indices.
%
%   The deposit stores event times in session-relative seconds, converted
%   from source frame indices against these same timestamps, so the inverse
%   is a nearest-timestamp lookup and is exact. Rows are returned as row
%   vectors of doubles, matching what the source cell arrays held.
seconds = h5read(nwbFile, dataPath);
index = double(h5read(nwbFile, [dataPath '_index']));
edges = [0; index(:)];

% Nearest-timestamp lookup on a strictly increasing axis. discretize with
% midpoint edges is O(n log m) and avoids interp1's copy of the whole axis.
mid = (tVolt(1:end-1) + tVolt(2:end)) / 2;
allFrames = discretize(seconds, [-inf; mid(:); inf]);

frames = cell(numel(index), 1);
for k = 1:numel(index)
    frames{k} = reshape(allFrames(edges(k)+1 : edges(k+1)), 1, []);
end
end
