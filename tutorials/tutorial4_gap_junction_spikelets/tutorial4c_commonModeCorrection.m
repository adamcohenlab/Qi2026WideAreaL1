%% Tutorial 4c — removing the distance trend and the common mode
%
% Stage 4b left two artifacts on top of the spikelet: a smooth trend that
% grows towards short distances, and a pedestal shared by the whole field.
% This stage removes them, in that order, and checks the result against the
% deposited `spikelet_amplitude_normalized`.
%
% The design constraint that shapes everything here: both corrections are
% estimated from data that contains the signal being looked for. If the
% distance trend is fitted through the peak, then a real spikelet in the
% nearest distance bin becomes part of the trend and is subtracted away, and
% the method quietly guarantees itself a null result. So both fits exclude
% the peak region, and the correction is deliberately not judged by how flat
% it makes the output.
%
% Runtime: about 2 min, dominated by the one-off spike subtraction.

clear; close all

%% ----------------------------------------------------------- configuration
% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

cfg = struct();
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
cfg.nBack = 20;
cfg.nFront = 20;
cfg.preCells = [25 100 152 316];    % presynaptic cells to correct and check
cfg.examplePre = 152;

%% ------------------------------------------- 1. load, triggers, subtraction
data = loadSparsePulseHadForSpikelet(cfg.nwbFile);
trig = buildOrthogonalTriggers(data);
fprintf('\nspike subtraction (one-off, about 70 s):\n');
tracesSub = subtractSpikes(data.traces_all, data.spk_t);
data.traces_all = [];
traceMat = tracesSub';   % T x nCells, kept single: 4 GB, not 8
clear tracesSub

%% ------------------------------------------------- 2. correct, cell by cell
% The correction is row-wise: everything it needs is one presynaptic cell's
% relationship to all 320 post cells. That is also why it cannot be run from
% the deposited waveform tensor, which stops at 400 um and stores the rest
% as zeros — the far field it needs is not there.
tau = -cfg.nBack:cfg.nFront;
results = struct('pre', {}, 'staRow', {}, 'out', {}, 'ampRaw', {}, 'ampCor', {});

for k = 1:numel(cfg.preCells)
    ii = cfg.preCells(k);
    staRow = spikeletStaRow(ii, trig, traceMat, ...
        'nBack', cfg.nBack, 'nFront', cfg.nFront, 'Verbose', true);

    out = correctCommonMode(staRow, data.distMat(ii,:), ...
        'nBack', cfg.nBack, 'dx', data.dx, 'selfIndex', ii);

    results(k).pre = ii;
    results(k).staRow = staRow;
    results(k).out = out;
    results(k).ampRaw = spikeletAmplitude(staRow - mean(staRow,2), cfg.nBack);
    results(k).ampCor = spikeletAmplitude(out.trRes, cfg.nBack);
end

%% ------------------------------------------------------ 3. look at one cell
% Four panels: the waveforms after the distance trend comes out, the fitted
% trend itself per distance bin, the same waveforms with the far-field
% common mode drawn over them, and what is left.
kEx = find(cfg.preCells == cfg.examplePre, 1);
plotSpikeletDiagnostics('correction', results(kEx).out, ...
    data.distMat(cfg.examplePre,:), data.dx, cfg.examplePre, data.dt);
plotSpikeletDiagnostics('amplitudeVsDistance', data.distMat(cfg.examplePre,:), ...
    results(kEx).ampRaw, results(kEx).ampCor, data.dx);

% The clamp is worth seeing. c is the per-post-cell regression coefficient
% onto the far-field waveform, fitted off-peak and then clamped at zero. A
% negative coefficient would mean the post cell carried an inverted copy of
% the common mode, which is not physical here, and subtracting a negative
% multiple would ADD a bump at offset zero and manufacture spikelets.
c = results(kEx).out.c;
fprintf('\n--- common-mode coefficients, pre cell %d ---\n', cfg.examplePre);
fprintf('  median %.3f, IQR [%.3f %.3f]\n', median(c), quantile(c,0.25), quantile(c,0.75));
fprintf('  clamped to zero: %d of %d post cells\n', sum(c == 0), numel(c));

%% -------------------------------------- 4. check against the deposit
% `spikelet_amplitude_normalized` is the published corrected amplitude,
% `spkletHadAmpCor`, and unlike the waveform tensor it covers all 320 x 320
% pairs. The comparison band is the one the detection family uses, under
% 400 um.
%
% Exact agreement is not expected and would be suspicious: stage 4a showed
% the reconstruction carries 0.0119 per cent more trigger events than the
% deposit, because the source's spike motion correction cannot be
% reproduced without `mcTrace_all`. That difference propagates here.
fprintf('\n--- corrected amplitude versus the deposit ---\n');
fprintf('%5s %7s | %9s | %10s %10s | %10s\n', ...
    'pre', 'nNear', 'corr', 'max|diff|', 'med|diff|', 'med|dep|');
for k = 1:numel(results)
    ii = results(k).pre;
    near = data.distMat(ii,:)' * data.dx < 400 & (1:data.nCells)' ~= ii;
    a = results(k).ampCor(near);
    b = data.deposited.spkletHadAmpCor(ii, near)';
    r = corrcoef(a, b);
    fprintf('%5d %7d | %9.6f | %10.3e %10.3e | %10.3e\n', ...
        ii, sum(near), r(1,2), max(abs(a-b)), median(abs(a-b)), median(abs(b)));
end

figure('Color','w','Name','corrected amplitude vs deposit');
allOurs = []; allDep = [];
for k = 1:numel(results)
    ii = results(k).pre;
    near = data.distMat(ii,:)' * data.dx < 400 & (1:data.nCells)' ~= ii;
    allOurs = [allOurs; results(k).ampCor(near)];            %#ok<AGROW>
    allDep  = [allDep;  data.deposited.spkletHadAmpCor(ii, near)']; %#ok<AGROW>
end
plot(allDep, allOurs, '.', 'MarkerSize', 10); hold on
lim = [min([allDep;allOurs]) max([allDep;allOurs])];
plot(lim, lim, 'r'); axis equal; xlim(lim); ylim(lim)
xlabel('Deposited spkletHadAmpCor'); ylabel('Recomputed')
title(sprintf('%d pairs under 400 \\mum, %d presynaptic cells', ...
    numel(allDep), numel(results)))

%% -------------------------------------------- 5. what the correction costs
% The correction is not free, and the honest way to see that is to ask what
% it does to pairs that are almost certainly uncoupled. It should pull the
% far field to zero, leave the near field with whatever is real, and it
% should NOT simply scale everything down.
fprintf('\n--- effect of the correction, pooled over %d presynaptic cells ---\n', ...
    numel(results));
edges = [0 100 200 400 800 1600 3000];
for k = 1:numel(edges)-1
    raw = []; cor = [];
    for m = 1:numel(results)
        ii = results(m).pre;
        d = data.distMat(ii,:)' * data.dx;
        sel = d >= edges(k) & d < edges(k+1) & (1:data.nCells)' ~= ii;
        raw = [raw; results(m).ampRaw(sel)];    %#ok<AGROW>
        cor = [cor; results(m).ampCor(sel)];    %#ok<AGROW>
    end
    if ~isempty(raw)
        fprintf('  %4d - %4d um : n = %4d, median raw %+.3e -> corrected %+.3e\n', ...
            edges(k), edges(k+1), numel(raw), median(raw), median(cor));
    end
end
fprintf(['\n  The near bins keep a clearly positive median and the far bins\n', ...
         '  collapse. That asymmetry is the point: a correction that flattened\n', ...
         '  everything would have removed the signal along with the artifact.\n']);
%
% The far bins do not land exactly on zero; they go slightly negative, by a
% few times 1e-4 against near-field amplitudes of 1e-2. That is the
% over-subtraction anticipated in stage 4b: the common mode is estimated
% beyond 800 um, where the distance trend has flattened but not quite
% vanished, so a little of the trend is absorbed into the pedestal and
% subtracted from everyone. It biases the far field low, which is the safe
% direction for a one-sided test, and stage 4d does not rely on the
% corrected far field anyway — its null is built from the uncorrected
% traces.

save(fullfile(tempdir, 'tutorial4c_state.mat'), 'cfg', 'results', '-v7.3');
fprintf('\nStage 4c done. Next: tutorial4d_spikeletDetection.\n');
