%% Tutorial 4a — orthogonal Hadamard stimulation, and what counts as a trigger
%
% Two cells are electrically coupled if one's spike puts a small, fast,
% passively filtered copy of itself into the other: a spikelet. Finding
% spikelets optically across 320 cells needs a way to be sure that a
% deflection locked to cell i's spike is not simply cell j responding to the
% same light that made cell i fire.
%
% That is what the Hadamard block is for. In every fourth recording chunk
% each cell's blue stimulation follows one row of a Hadamard pattern, so at
% any instant roughly half the field is illuminated and which half keeps
% changing. For any pair there are therefore plenty of frames where the pre
% cell is lit and the post cell is dark. Restrict to those, and a shared
% response to shared light is not available as an explanation.
%
% This stage builds the trigger sets and checks them against the deposit,
% which stores them directly.
%
% Runtime: about 60 s, plus 10 s for the load.
%
% Requires: get_sta_mat_single on the path (see the folder README).

clear; close all

%% ----------------------------------------------------------- configuration
% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

cfg = struct();
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
cfg.nBack = 20;              % tutorial STA window; the deposit stores +/-140
cfg.nFront = 20;
cfg.nCheckPairs = 8;         % pairs for the waveform cross-check
cfg.exampleChunk = 4;        % a Hadamard chunk to draw

%% -------------------------------------------------------------- 1. load
% LoadDepositedTriggers pulls in the ragged per-pair spike column, 39 M
% events. Only this stage needs it.
data = loadSparsePulseHadForSpikelet(cfg.nwbFile, 'LoadDepositedTriggers', true);

%% ------------------------------------------- 2. build the trigger sets
trig = buildOrthogonalTriggers(data);

% Look at one pair over one Hadamard chunk. The grey band is where the pre
% cell is lit and the post cell is not; grey dots are all of the pre cell's
% spikes inside a shrunk Hadamard step, red dots are the ones that survive
% the post cell's edge guard and the chunk-boundary guard.
examplePre = 152;
examplePost = 285;
plotSpikeletDiagnostics('orthogonality', trig, examplePre, examplePost, ...
    cfg.exampleChunk, data.nFrame);

%% ----------------------------------- 3. check the triggers against the deposit
% /analysis/pairwise_connectivity_metrics/hadamard_orthogonal_spike_times
% holds the published trigger set for every pair closer than 400 um. This is
% the strongest check available on the construction above, because it is the
% same quantity rather than a downstream summary.
%
% One difference is expected and cannot be removed. The source ran
% `spike_times_motion_correction_aggressive` before building the Hadamard
% spike matrix, and its input `mcTrace_all` is not in the deposit, so the
% reconstruction starts from the uncorrected /units/spike_times. The
% reconstruction should therefore be a strict SUPERSET of the deposit: the
% correction only removes events.
fprintf('\n--- trigger sets versus the deposit ---\n');
pairMask = data.deposited.computedMask & ~logical(eye(data.nCells));
[it, jt] = find(pairMask);
nExact = 0; nDep = 0; nGot = 0; nMissing = 0; nExtra = 0;
worstJaccard = 1; worstPair = [0 0];
tic
for k = 1:numel(it)
    a = double(data.deposited.spk_t_had_ortho{it(k), jt(k)});
    b = trig.pairTriggers(it(k), jt(k));
    nDep = nDep + numel(a);
    nGot = nGot + numel(b);
    if isequal(sort(a(:))', sort(b(:))')
        nExact = nExact + 1;
    else
        nMissing = nMissing + numel(setdiff(a, b));
        nExtra   = nExtra   + numel(setdiff(b, a));
        j = numel(intersect(a,b)) / numel(union(a,b));
        if j < worstJaccard, worstJaccard = j; worstPair = [it(k) jt(k)]; end
    end
end
fprintf('  %d pairs compared in %.0f s\n', numel(it), toc);
fprintf('  exact set match     : %d (%.2f%%)\n', nExact, 100*nExact/numel(it));
fprintf('  deposited events    : %d\n', nDep);
fprintf('  reconstructed events: %d (%+.4f%%)\n', nGot, 100*(nGot-nDep)/nDep);
fprintf('  missing / extra     : %d / %d\n', nMissing, nExtra);
if nMissing == 0
    fprintf('  -> strict superset, as predicted. The extras are the events the\n');
    fprintf('     motion correction removed, %.4f%% of the total.\n', 100*nExtra/nDep);
else
    warning('Tutorial4:MissingTriggers', ...
        ['%d deposited events are missing from the reconstruction. That is not ', ...
         'the motion correction; check the guard construction.'], nMissing);
end
if worstPair(1)
    fprintf('  worst pair          : %d -> %d, Jaccard %.5f\n', ...
        worstPair(1), worstPair(2), worstJaccard);
end

%% ----------------------------- 4. check the averaging against the deposit
% The waveform check has to be separated from the trigger check, or the
% 0.01 per cent trigger difference hides whether the averaging itself is
% right. So: feed the DEPOSITED trigger times into the tutorial's own
% spike-subtracted trace and compare with the deposited waveform. Anything
% above floating-point noise here is a real disagreement.
%
% The deposit stores +/-140 samples and this stage computes +/-20. Those are
% comparable because the chunk-boundary guard leaves no trigger within 140
% frames of either end of the recording, so neither window drops events that
% the other keeps, and each offset is averaged independently.
fprintf('\n--- waveform reconstruction versus the deposit ---\n');
nEvent = cellfun(@numel, data.deposited.spk_t_had_ortho);
candidate = find(pairMask & nEvent > 3000);
rng(1, 'twister');
candidate = candidate(randperm(numel(candidate), min(cfg.nCheckPairs, numel(candidate))));
[checkPre, checkPost] = ind2sub([data.nCells data.nCells], candidate);

postCells = unique(checkPost);
tracesSubCheck = subtractSpikes(data.traces_all(postCells,:), data.spk_t(postCells), ...
    'Verbose', false);
rowOf = zeros(data.nCells,1);
rowOf(postCells) = 1:numel(postCells);

tau = -cfg.nBack:cfg.nFront;
centre = data.deposited.tBack + 1 + tau;
maxDiff = zeros(numel(checkPre),1);
for k = 1:numel(checkPre)
    ii = checkPre(k); jj = checkPost(k);
    trace = double(tracesSubCheck(rowOf(jj),:))';
    ours = squeeze(get_sta_mat_single( ...
        {double(data.deposited.spk_t_had_ortho{ii,jj})}, trace, [cfg.nBack cfg.nFront], 1));
    ref = data.deposited.readHadSta(ii,jj);
    ref = ref(centre);
    maxDiff(k) = max(abs(ours(:) - ref(:)));
    fprintf('  %3d -> %3d  %6d events  peak |waveform| %.3e  max |diff| %.3e\n', ...
        ii, jj, nEvent(ii,jj), max(abs(ref)), maxDiff(k));
    if k == 1
        plotSpikeletDiagnostics('staCheck', tau, ours, ref, ii, jj, data.dt);
    end
end
fprintf('  worst over %d pairs: %.3e\n', numel(checkPre), max(maxDiff));
if max(maxDiff) < 1e-6
    fprintf('  -> exact to single precision. The averaging is the published one.\n');
end

%% ------------------------------------------------ 5. what this bought us
% The orthogonality requirement is not free: it throws away every spike that
% happened while the post cell was also lit, which is about half of them,
% and the guards remove more. The count that survives is a property of the
% pair, which is why the deposit stores it per pair rather than per cell.
nHad = cellfun(@numel, trig.tHad);
nPair = zeros(numel(it),1);
for k = 1:numel(it)
    nPair(k) = numel(trig.pairTriggers(it(k), jt(k)));
end
fprintf('\n--- trigger budget ---\n');
fprintf('  spikes per cell, all epochs        : median %d\n', median(cellfun(@numel, data.spk_t)));
fprintf('  inside a shrunk Hadamard step      : median %d (%.0f%%)\n', ...
    median(nHad), 100*median(nHad)/median(cellfun(@numel, data.spk_t)));
fprintf('  surviving a pair''s orthogonality   : median %d (%.0f%% of the Hadamard spikes)\n', ...
    median(nPair), 100*median(nPair)/median(nHad));

fprintf('\nStage 4a done. Next: tutorial4b_spikeletWaveform.\n');
