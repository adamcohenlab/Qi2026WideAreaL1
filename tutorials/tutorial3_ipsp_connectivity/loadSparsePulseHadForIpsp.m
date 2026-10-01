function data = loadSparsePulseHadForIpsp(nwbFile, varargin)
%LOADSPARSEPULSEHADFORIPSP Load a sparsePulseHad NWB file for Tutorial 3.
%
%   data = LOADSPARSEPULSEHADFORIPSP(nwbFile) returns a struct whose fields
%   carry the *source pipeline's* variable names, so the tutorial bodies read
%   like the original scripts and can call the original helper functions
%   (get_sta_mat_raw, get_sta_mat_2, get_sta_mat_single) unchanged.
%
%   data = LOADSPARSEPULSEHADFORIPSP(nwbFile, 'Name', value, ...) accepts:
%
%     'LoadTraces'   true (default) loads the full nCells-by-nFrames voltage
%                    matrix into memory, about 4 GB as single. Set false to get
%                    metadata, blue stimulation and t_blue only.
%     'Verbose'      true (default) prints a short load report.
%
%   FIELDS
%     traces_all     nCells x nFrames single. This is the deposited
%                    normalized_voltage, i.e. the source traces_all_bc divided
%                    by the per-cell spkHgtNoBlue. See NORMALISATION below.
%     blue_all       nCells x nFrames single, source scale retained.
%     cell_coord     nCells x 2, [x_pixel y_pixel], source column order.
%     distMat        nCells x nCells, pixels, read from the deposit.
%     t_blue         nCells x 1 cell array of sparse-pulse trigger frames,
%                    CONVENTION A (see below).
%     nSpikesPulse   deposited estimate, for validation in stage 1.
%     dt, dx, nCells, nFrame, nChunk, T
%     deposited      struct of stored analysis products used as ground truth:
%                    pMat, ipspAmpTestN2, ipspDecay, seTau, and
%                    readWaveform(pre,post) / readCrosstalk(pre,post).
%
%   NORMALISATION
%   The deposit stores normalized voltage, not raw. From dataPrepsForMAT.m:
%
%       traces_all_n               = traces_all_bc ./ spkHgtNoBlue
%       crossBluePulseStaFNCorFull = crossBluePulseStaFCorFull ...
%                                      ./ spkHgtNoBlue'(post) ./ nSpikesPulse(pre)
%
%   so an STA computed from data.traces_all is already divided by the post
%   cell's spkHgtNoBlue but NOT by the pre cell's nSpikesPulse. Any waveform
%   that must be compared with the deposited tensor, INCLUDING EVERY NULL
%   WAVEFORM, has to be divided by nSpikesPulse(pre) as well. Getting this
%   wrong silently invalidates every p-value. spkHgtNoBlue itself is not in the
%   deposit, so the tutorial works in normalised units throughout; detection and
%   kinetics are invariant to that per-pair scale and amplitudes come out in the
%   normalised units the manuscript reports.
%
%   TRIGGER CONVENTION
%   t_blue uses convention A, find(diff(blue_i>0 & stimTypeMask)>0), so tau = 0
%   is the last frame BEFORE the pulse. This is the convention behind the
%   deposited waveform tensor and behind crsTalkTemp's 2-frame onset. The source
%   pipeline used a different, one-frame-later convention for the null set only;
%   that inconsistency is resolved here in favour of A. See
%   ISSUE_tblue_trigger_convention.md.
%
%   See also TUTORIAL3A_BLUEPULSESPIKECOUNT.

%% ---------------------------------------------------------------- options
opt = struct('LoadTraces', true, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
assert(exist(nwbFile, 'file') == 2, 'Tutorial3:MissingFile', ...
    'NWB file not found: %s', nwbFile);

timer = tic;

%% ------------------------------------------------------------ dimensions
info = h5info(nwbFile, '/processing/optical_voltage/normalized_voltage/data');
sz = info.Dataspace.Size;          % MATLAB order: [nCells nFrames]
nCells = sz(1);
T = sz(2);

data = struct();
data.nCells = nCells;
data.T = T;
data.dt = 1.27e-3;                 % confirmed voltage-camera exposure
data.dx = 6.5;                     % um per pixel
data.nwbFile = nwbFile;

% Chunk structure comes from the deposit rather than the source literals
% nFrame = 4960 / 640 chunks, so the tutorial is not tied to one session.
chunkStart = h5read(nwbFile, '/intervals/recording_chunks/start_time');
data.nChunk = numel(chunkStart);
assert(mod(T, data.nChunk) == 0, 'Tutorial3:ChunkStructure', ...
    'Frame count %d is not a whole multiple of the %d recording chunks.', T, data.nChunk);
data.nFrame = T / data.nChunk;

%% ------------------------------------------------------------ unit table
xPixel = h5read(nwbFile, '/units/x_pixel');
yPixel = h5read(nwbFile, '/units/y_pixel');
data.cell_coord = [xPixel(:), yPixel(:)];
data.nSpikesPulse = double(h5read(nwbFile, '/units/estimated_spikes_per_blue_pulse'));
data.nSpikesPulse = data.nSpikesPulse(:);
data.gfpIdx = logical(h5read(nwbFile, '/units/gfp_positive'));
data.idxUse = logical(h5read(nwbFile, '/units/passed_qc'));

% distMat is deposited; cross-check it against the coordinates so a coordinate
% mix-up cannot pass silently.
data.distMat = reshape(double(h5read(nwbFile, ...
    '/analysis/pairwise_connectivity_metrics/distance_pixels')), nCells, nCells);
distCheck = squareform(pdist(data.cell_coord));
relErr = max(abs(data.distMat(:) - distCheck(:))) / max(distCheck(:));
assert(relErr < 1e-5, 'Tutorial3:DistanceMismatch', ...
    'Deposited distMat disagrees with pdist(cell_coord); relative error %.2e.', relErr);

%% ----------------------------------------------------- blue stimulation
data.blue_all = h5read(nwbFile, '/stimulus/presentation/blue_stimulation_per_cell/data');

%% ------------------------------------------- t_blue, convention A
% Reproduced verbatim from synapticConn_baselineGen.m / bluePulseSpikeCount.m.
% stimTypeIdx{1} is the sparse-pulse epoch: everything except the last 11 frames
% of each chunk and the whole of every fourth (Hadamard) chunk.
nFrame = data.nFrame;
stimTypeIdx = { ...
    setdiff((1:T), [reshape((-10:0) + (nFrame:nFrame:T)',1,[]) ...
                    reshape((1:nFrame) + (3*nFrame:nFrame*4:T-nFrame)',1,[])]), ...
    (1:nFrame) + (3*nFrame:nFrame*4:T)' };
stimTypeIdx{1} = setdiff(stimTypeIdx{1}, stimTypeIdx{2}(:));

stimTypeMask = false(1,T);
stimTypeMask(stimTypeIdx{1}(:)) = true;
data.stimTypeMask = stimTypeMask;

data.t_blue = cell(nCells,1);
for ii = 1:nCells
    blue_i = data.blue_all(ii,:);
    data.t_blue{ii} = find(diff(blue_i>0 & stimTypeMask) > 0);
end
data.t_blueConvention = ['A: find(diff(blue_i>0 & stimTypeMask)>0); tau = 0 is ' ...
    'the last frame before the pulse. Source: sparsePulseHad_genStaMat.m, ' ...
    'synapticConn_baselineGen.m, bluePulseSpikeCount.m.'];

% This is the ONLY trigger definition the loader supplies, and it is the one
% behind every deposited product: nSpikesPulse, the pair waveform tensor, and
% the far-pulse baseline. Stage 3d deliberately builds a DIFFERENT one locally
% for its null set, one frame later, because the published detection script
% does. That quirk is confined to the detection script in the source pipeline,
% so it is confined to tutorial3d here too rather than being hoisted up to a
% shared field where other stages might pick it up by accident. See the comment
% at its point of use in tutorial3d_detectAndFitIpsp.m.

% tRm / removedTimes: chunk-start frames excluded from far-pulse null sampling.
data.removedTimes = (1:nFrame:T) + (0:160)';

%% --------------------------------------------------- deposited products
base = '/analysis/pairwise_connectivity_metrics/';
readPair = @(name) reshape(double(h5read(nwbFile, [base name])), nCells, nCells);
dep = struct();
dep.pMat = logical(readPair('synaptic_connection'));
dep.ipspAmpTestN2 = readPair('ipsp_amplitude_norm_per_presynaptic_spike');
dep.ipspDecay = readPair('ipsp_decay_source_units');
dep.seTau = readPair('ipsp_decay_standard_error_source_units');

% Lazy single-pair readers for the 3-D tensors. MATLAB order is
% [pre post time]; the tensors are chunked [nCells nCells 25], so a single-pair
% read still touches whole chunks. Read once and cache if you need many pairs.
wavePath = '/processing/connectivity_mapping/blue_pulse_ipsp_waveform/data';
xtalkPath = '/processing/connectivity_mapping/blue_light_crosstalk_estimate/data';
waveInfo = h5info(nwbFile, wavePath);
nTau = waveInfo.Dataspace.Size(3);
dep.nTau = nTau;
dep.tBack = floor((nTau-1)/2);
dep.tFront = nTau - dep.tBack - 1;
dep.readWaveform = @(pre,post) squeeze(h5read(nwbFile, wavePath, [pre post 1], [1 1 nTau]));
dep.readCrosstalk = @(pre,post) squeeze(h5read(nwbFile, xtalkPath, [pre post 1], [1 1 nTau]));
dep.wavePath = wavePath;
data.deposited = dep;

data.tBack = dep.tBack;
data.tFront = dep.tFront;

%% ---------------------------------------------------------------- traces
if opt.LoadTraces
    data.traces_all = h5read(nwbFile, '/processing/optical_voltage/normalized_voltage/data');
else
    data.traces_all = [];
end

%% ---------------------------------------------------------------- report
if opt.Verbose
    nEvent = cellfun(@numel, data.t_blue);
    fprintf('Loaded %s\n', nwbFile);
    fprintf('  %d cells, %d frames, %d chunks of %d frames, dt = %.2f ms\n', ...
        nCells, T, data.nChunk, data.nFrame, data.dt*1e3);
    fprintf('  sparse-pulse epoch: %.1f%% of frames\n', 100*mean(stimTypeMask));
    fprintf('  t_blue (convention A): median %d pulses/cell [%d, %d]\n', ...
        median(nEvent), min(nEvent), max(nEvent));
    fprintf('  deposited pMat: %d positive directed pairs\n', sum(dep.pMat(:)));
    if opt.LoadTraces
        w = whos('data');
        fprintf('  traces loaded, struct is %.2f GB\n', w.bytes/2^30);
    else
        fprintf('  traces NOT loaded (LoadTraces false)\n');
    end
    fprintf('  %.1f s\n', toc(timer));
end
end
