%% Tutorial 3b -- the pair waveform, and the two things subtracted from it
%
% The raw quantity behind every connectivity claim in this study is simple:
% average cell j's voltage around the blue pulses delivered to cell i. If i
% inhibits j there is a small negative deflection a few milliseconds in.
%
% The difficulty is that this average is dominated by things that have nothing
% to do with synapses. Scattered blue light and direct optical crosstalk reach
% every cell on every pulse, producing a POSITIVE deflection far larger than
% any IPSP, for every pair regardless of connectivity. Two corrections are
% applied before anything is fitted:
%
%   1. FAR-PULSE BASELINE. Average cell j's response to pulses delivered
%      >800 um away, where no synaptic effect is plausible, and subtract it.
%      This removes the field-wide artifact. Indexed by the RECORDED cell.
%
%   2. NEIGHBOUR-ANNULUS MEDIAN. For each stimulated cell i, take the median
%      response across cells 300-500 um from i and subtract that too. This
%      removes what remains of i's own scattered light at intermediate range.
%      Indexed by the STIMULATED cell.
%
% The result is crossBluePulseStaFCorFull, the input to detection in 3d.
%
% Sources: synapticConn_baselineGen.m (stage 1 and the Monte Carlo baseline)
% and crosstalkTemplate_nullSet.m (the annulus median).

clear; clc

%% 1. Configuration
cfg = struct;
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
cfg.preCells = [3 17 42 88 151 260];   % presynaptic cells to reproduce
cfg.plotPre  = 3;                      % which one to plot, or []

% Optional. results.mat is NOT part of the deposit. When available it supplies
% spkHgtNoBlueAvg, which lets section 7 additionally reproduce the PUBLISHED
% annulus median exactly and so validate the STA and far-pulse baseline
% independently of the normalisation-order question. Set to '' to skip.
cfg.resultsMat = '';

% Radii, in pixels, exactly as in the source.
cfg.rCutNear = 400 / 6.5;   % pair inclusion, and near-exclusion for the baseline
cfg.rCutFar  = 800 / 6.5;   % "far enough that a synapse is implausible"
cfg.rNear    = 300 / 6.5;   % annulus inner radius
cfg.rNeigh   = 500 / 6.5;   % annulus outer radius

% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

%% 2. Load
data = loadSparsePulseHadForIpsp(cfg.nwbFile);
nCells = data.nCells;
tBack = data.tBack; tFront = data.tFront;
nTau = tBack + tFront + 1;
tau = (-tBack:tFront) * data.dt;
pre = cfg.preCells(:)';

%% 3. Raw pair waveform for the selected presynaptic cells
% The source computes the whole tensor with
%     get_sta_mat_2(t_blue(:,1), traces_all', [tBack tFront], 1)
% which loops over presynaptic cells internally. get_sta_mat_single does the
% same arithmetic for a chosen subset of trigger sets, so
%     get_sta_mat_single(t_blue(S), traces_all', ...)
% equals rows S of the full tensor. Section 7 checks that equivalence.
fprintf('\nRaw pair waveform for %d presynaptic cells...\n', numel(pre));
timer = tic;
crossBluePulseStaF = get_sta_mat_single(data.t_blue(pre), data.traces_all', ...
    [tBack tFront], 1);
fprintf('  %.1f s\n', toc(timer));

%% 4. Correction 1: far-pulse baseline
% Needed for EVERY recorded cell, not just the selected ones, because the
% annulus median in section 5 averages over many recorded cells.
fprintf('\nFar-pulse baseline for all %d cells...\n', nCells);
timer = tic;
bluePulseBaseline = farPulseBaseline(1:nCells, data, tBack, tFront, ...
    cfg.rCutNear, cfg.rCutFar, true);
fprintf('  %.1f s total\n', toc(timer));

% Indexed by recorded cell, so it broadcasts along the presynaptic axis.
crossBluePulseStaFCor = crossBluePulseStaF - shiftdim(bluePulseBaseline, -1);

%% 5. Correction 2: neighbour-annulus median
% For stimulated cell ii, the reference is the median response across recorded
% cells 300-500 um FROM THAT STIMULATED CELL. It does not depend on the
% recorded cell, so one reference is subtracted from that whole row.
%
% (The source wraps this in a second loop over the recorded cell, but nothing
% in the loop body uses that index, so it recomputes the same value 320 times.
% Collapsing it is numerically identical. See README.)
%
% NORMALISATION ORDER. The median runs ACROSS recorded cells, and those cells
% have different spike heights. So it matters whether the median is taken
% before or after dividing by spkHgtNoBlue, and the two do not commute:
%
%   published :  median of UNNORMALISED, then divide     (what detection used)
%   intended  :  divide first, then median               (equal weight per cell)
%
% The published pipeline did the first; it should have done the second. The
% tutorial implements the INTENDED version -- the one deliberate departure from
% published behaviour in Tutorial 3. It is also what falls out naturally here,
% since the deposit is already normalised. See
% ISSUE_annulus_median_normalisation.md.
nearMask = data.distMat(pre,:) < cfg.rCutNear;
annulus = @(ii) data.distMat(ii,:) < cfg.rNeigh & data.distMat(ii,:) > cfg.rNear;

crossBluePulseStaFCorFull = zeros(numel(pre), nCells, nTau);
for k = 1:numel(pre)
    m = annulus(pre(k));
    if sum(m) == 0
        trRef = zeros(1, nTau);
    else
        trRef = squeeze(nanmedian(crossBluePulseStaFCor(k, m, :), 2))';
    end
    crossBluePulseStaFCorFull(k,:,:) = squeeze(crossBluePulseStaFCor(k,:,:)) - trRef;
end
crossBluePulseStaFCorFull = crossBluePulseStaFCorFull .* nearMask;

%% 6. Normalise onto the deposited scale
% data.traces_all is normalized_voltage, already divided by the recorded cell's
% spkHgtNoBlue. The remaining factor is the stimulated cell's nSpikesPulse.
% Getting this wrong is the single easiest way to invalidate everything
% downstream -- see the loader's NORMALISATION note.
crossBluePulseStaFNCorFull = crossBluePulseStaFCorFull ./ data.nSpikesPulse(pre);
crossBluePulseStaFNCorFull(isinf(crossBluePulseStaFNCorFull)) = nan;

%% 7. Check against the deposit
fprintf('\n--- reproduction of blue_pulse_ipsp_waveform ---\n');

% 7a. get_sta_mat_single vs get_sta_mat_2 on one presynaptic cell.
checkPre = pre(1);
tBlueFull = cell(nCells,1);
tBlueFull(:) = {zeros(1,0)};
tBlueFull{checkPre} = data.t_blue{checkPre};
viaMat2 = get_sta_mat_2(tBlueFull, data.traces_all', [tBack tFront], 1);
viaSingle = squeeze(crossBluePulseStaF(1,:,:));
d = max(abs(squeeze(viaMat2(checkPre,:,:)) - viaSingle), [], 'all');
fprintf('get_sta_mat_single vs get_sta_mat_2, cell %d: max abs diff %.3e\n', checkPre, d);

% 7b. CORRECTNESS CHECK -- the tutorial's own computation against the deposit.
% As of the 2026-09-26 reconversion the deposit carries the corrected annulus
% median (normalise, then take the median), which is what section 5 computes.
% So this is now a straight reproduction test and it must pass.
fprintf('\n(a) tutorial (normalise-then-median) vs the deposit\n');
fprintf('%6s %8s %14s %14s %10s\n', 'pre','nPairs','maxAbsDiff','relToRange','corr');
allRel = nan(numel(pre),1);
for k = 1:numel(pre)
    ii = pre(k);
    postIdx = find(data.distMat(ii,:) < cfg.rCutNear);
    mine = squeeze(crossBluePulseStaFNCorFull(k, postIdx, :));
    theirs = nan(numel(postIdx), nTau);
    for q = 1:numel(postIdx), theirs(q,:) = data.deposited.readWaveform(ii, postIdx(q)); end
    good = isfinite(mine) & isfinite(theirs);
    allRel(k) = max(abs(mine(good)-theirs(good))) / (max(theirs(good))-min(theirs(good)));
    fprintf('%6d %8d %14.3e %14.3e %10.6f\n', ii, numel(postIdx), ...
        max(abs(mine(good)-theirs(good))), allRel(k), corr(mine(good), theirs(good)));
end
fprintf('worst: %.3e\n', max(allRel));
if max(allRel) < 1e-4
    fprintf(['PASS: the STA, the far-pulse baseline and the annulus median all\n' ...
             'reproduce the deposit.\n']);
else
    fprintf('FAIL -- see TUTORIAL3_VALIDATION.md.\n');
end

% 7c. HISTORICAL COMPARISON, not a test. Rebuild the median the way the
% pipeline did before 2026-09-23 -- on unnormalised data -- to show how much
% that correction moved the waveform. Needs the per-cell spike heights, which
% are not in the deposit, so this branch is skipped without results.mat.
haveH = ~isempty(cfg.resultsMat) && exist(cfg.resultsMat, 'file') == 2;
if haveH
    R = load(cfg.resultsMat, 'spkHgtNoBlueAvg');
    h = double(R.spkHgtNoBlueAvg(:)).';    % 1 x nCells, indexed by RECORDED cell
    h(h == 0) = nan;

    publishedStyle = zeros(numel(pre), nCells, nTau);
    for k = 1:numel(pre)
        m = annulus(pre(k));
        C = squeeze(crossBluePulseStaFCor(k,:,:)) .* h.';   % back to unnormalised
        if sum(m) == 0
            M = zeros(1, nTau);
        else
            M = nanmedian(C(m,:), 1);
        end
        publishedStyle(k,:,:) = (C - M) ./ h.';
    end
    publishedStyle = publishedStyle .* nearMask ./ data.nSpikesPulse(pre);
    publishedStyle(isinf(publishedStyle)) = nan;

    fprintf('\n(b) the pre-2026-09-23 median (unnormalised) vs the deposit\n');
    fprintf('%6s %8s %14s %14s %10s\n', 'pre','nPairs','maxAbsDiff','relToRange','corr');
    relPub = nan(numel(pre),1);
    for k = 1:numel(pre)
        ii = pre(k);
        postIdx = find(data.distMat(ii,:) < cfg.rCutNear);
        a = squeeze(publishedStyle(k, postIdx, :));
        b = nan(numel(postIdx), nTau);
        for q = 1:numel(postIdx), b(q,:) = data.deposited.readWaveform(ii, postIdx(q)); end
        g = isfinite(a) & isfinite(b);
        relPub(k) = max(abs(a(g)-b(g))) / (max(b(g))-min(b(g)));
        fprintf('%6d %8d %14.3e %14.3e %10.6f\n', ii, numel(postIdx), ...
            max(abs(a(g)-b(g))), relPub(k), corr(a(g), b(g)));
    end
    fprintf('worst: %.3e of range\n', max(relPub));
    fprintf(['\nThat is the size of the annulus-median correction, NOT an error.\n' ...
             'Taking the median across neighbouring cells before dividing by each\n' ...
             'cell''s spike height weights the reference towards cells with large\n' ...
             'raw signals. Over the full 9124-pair family the correction moved 82\n' ...
             'connection calls and shifted amplitudes by a median 6.2%%. See\n' ...
             'ISSUE_annulus_median_normalisation.md.\n']);
else
    fprintf('\n(b) skipped: results.mat not available, so the older unnormalised\n');
    fprintf('    median cannot be reconstructed (needs spkHgtNoBlueAvg).\n');
end

%% 8. Optional: what each correction removed
if ~isempty(cfg.plotPre)
    k = find(pre == cfg.plotPre, 1);
    assert(~isempty(k), 'cfg.plotPre must be one of cfg.preCells.');
    ii = pre(k);
    postIdx = find(data.distMat(ii,:) < cfg.rCutNear);
    plotPairWaveformDiagnostics(ii, postIdx, tau, ...
        squeeze(crossBluePulseStaF(k, postIdx, :)), ...
        squeeze(crossBluePulseStaFCor(k, postIdx, :)), ...
        squeeze(crossBluePulseStaFCorFull(k, postIdx, :)), ...
        data.deposited.pMat(ii, postIdx));
end

%% 9. What to take away
% What is left after both corrections still contains a large positive crosstalk
% transient -- the corrections remove the FIELD-WIDE and MID-RANGE artifact, not
% the local one. Near pairs are exactly the pairs where scattered light is
% strongest, so the residual artifact is largest precisely where the synaptic
% signal lives.
%
% That is why detection cannot be "look for a negative deflection". Tutorial 3c
% builds an explicit model of the remaining crosstalk, and 3d fits it jointly
% with a strictly negative synaptic term, judged against a null built from the
% same far pulses used here.
