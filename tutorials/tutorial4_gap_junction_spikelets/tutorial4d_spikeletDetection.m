%% Tutorial 4d — is the spikelet real? The far-cell null, and BH
%
% Stage 4c produced a corrected amplitude per directed pair. It is a small
% number, and "small and positive" is not evidence on its own: the
% corrections are estimates, and a residue of either artifact would also be
% small and positive. What is needed is a distribution of amplitudes that
% this pair would have produced if there were no coupling, built from the
% same data so that everything not under test is matched rather than
% modelled.
%
% The surrogate: keep the post cell, keep the number of triggers, keep the
% optical condition, and change only WHICH CELL FIRED — to a cell too far
% from the post cell to be coupled to it. Anything shared between a near pre
% cell and a far one, including network synchrony and the post cell's own
% noise, appears in both and cancels.
%
% Runtime: about 2 min of setup, then roughly 3.4 s per pair at nRand =
% 1000 and 17 s per pair at the published nRand = 5000.

clear; close all

%% ----------------------------------------------------------- configuration
% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

cfg = struct();
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
cfg.nBack = 20;
cfg.nFront = 20;
cfg.nBack_c = 10;            % null window; source literal
cfg.preCells = [152 25];     % presynaptic cells to test
cfg.nPairsPerPre = 6;        % post cells per presynaptic cell
cfg.nRand = 1000;            % published run uses 5000; see runFullSpikeletDetection
cfg.rLim = 400;              % um; "far" for the null pool, and the BH family radius
cfg.seed = 0;

% The published BH cutoff over the full 9,124-pair family, measured from the
% source intermediates (crossHadOrtho_longDist.mat + crossHadOrthoSpklet_BH.mat)
% and recorded in TUTORIAL4_VALIDATION.md. It is used here ONLY to compare a
% subset's p-values against the deposited decisions; BH itself cannot be run
% on a subset, see section 5.
cfg.publishedBhCutoff = 0.024795;

%% ------------------------------------------- 1. load, triggers, subtraction
data = loadSparsePulseHadForSpikelet(cfg.nwbFile);
trig = buildOrthogonalTriggers(data);
fprintf('\nspike subtraction (one-off, about 70 s):\n');
tracesSub = subtractSpikes(data.traces_all, data.spk_t);
data.traces_all = [];
traceMat = tracesSub';   % T x nCells, kept single: 4 GB, not 8
clear tracesSub

%% ----------------------------------- 2. corrected amplitudes for those cells
corrected = containers.Map('KeyType','double','ValueType','any');
for ii = cfg.preCells
    staRow = spikeletStaRow(ii, trig, traceMat, ...
        'nBack', cfg.nBack, 'nFront', cfg.nFront, 'Verbose', true);
    out = correctCommonMode(staRow, data.distMat(ii,:), ...
        'nBack', cfg.nBack, 'dx', data.dx, 'selfIndex', ii);
    corrected(ii) = spikeletAmplitude(out.trRes, cfg.nBack);
end

%% --------------------------------------------------- 3. the null, pair by pair
rng(cfg.seed, 'twister');
rows = struct('pre',{},'post',{},'dist',{},'nSpk',{},'nCtrl',{}, ...
    'ampCor',{},'nullMean',{},'nullSd',{},'p',{},'ours',{},'deposited',{});

fprintf('\n--- per-pair test, nRand = %d ---\n', cfg.nRand);
fprintf('%5s %5s %7s %7s %8s | %10s %10s %10s | %7s %5s %5s\n', ...
    'pre','post','dist','nSpk','nCtrl','ampCor','nullMean','nullSD','p','ours','dep');
for ii = cfg.preCells
    amp = corrected(ii);
    near = find(data.distMat(ii,:)' * data.dx < cfg.rLim & (1:data.nCells)' ~= ii);
    pick = near(round(linspace(1, numel(near), min(cfg.nPairsPerPre, numel(near)))));

    for jj = pick(:)'
        ns = spikeletNullSet(ii, jj, trig, traceMat.', data.distMat, ...
            'nRand', cfg.nRand, 'rLim', cfg.rLim, 'nBack_c', cfg.nBack_c, ...
            'dx', data.dx, 'Seed', cfg.seed);

        % Right-tailed empirical p, written exactly as dataPrepsForMAT.m
        % does: the count of nulls at or above the observation, over
        % nRand + 1. Note the numerator has no +1, so p can be exactly 0.
        p = sum(amp(jj) < ns.nullAmp) / (cfg.nRand + 1);

        r = struct('pre',ii,'post',jj,'dist',data.distMat(ii,jj)*data.dx, ...
            'nSpk',ns.nSpk,'nCtrl',ns.nCtrl,'ampCor',amp(jj), ...
            'nullMean',mean(ns.nullAmp),'nullSd',std(ns.nullAmp),'p',p, ...
            'ours',p < cfg.publishedBhCutoff, ...
            'deposited',data.deposited.pMat_spklet(ii,jj));
        rows(end+1) = r; %#ok<SAGROW>

        fprintf('%5d %5d %7.0f %7d %8d | %10.2e %10.2e %10.2e | %7.4f %5d %5d\n', ...
            ii, jj, r.dist, r.nSpk, r.nCtrl, r.ampCor, r.nullMean, r.nullSd, ...
            p, r.ours, r.deposited);

        if ii == cfg.preCells(1) && jj == pick(1)
            plotSpikeletDiagnostics('nullPair', ns.nullAmp, amp(jj), p, ii, jj);
        end
    end
end

agree = sum([rows.ours] == [rows.deposited]);
fprintf('\n  decisions matching the deposit: %d of %d\n', agree, numel(rows));

%% ----------------------------------- 4. what the null does and does not hold
% Worth reading off the table above: nullMean is close to zero and much
% smaller than nullSD. The null is NOT centred on the common mode that stage
% 4c removed, even though it is built from uncorrected traces.
%
% That is not an accident and it is the reason the asymmetry is coherent.
% The common mode is specific to the PRE cell: it is that cell's spikes
% being correlated with network activity. The null pools triggers across
% every far cell, tens of thousands of events from hundreds of cells, which
% averages that cell-specific structure away. So the null measures the
% dispersion of an amplitude estimate at this trigger count, and the
% correction removes a bias the null could not have captured. They are
% doing different jobs, in series.
nm = [rows.nullMean]; nsd = [rows.nullSd];
fprintf('\n--- null diagnostics ---\n');
fprintf('  |nullMean| / nullSD : median %.3f, max %.3f\n', ...
    median(abs(nm./nsd)), max(abs(nm./nsd)));
%
% The second thing the table shows is why the post cell has to be held
% fixed. Trigger counts barely vary across these pairs, yet nullSD varies
% severalfold. It is dominated by the post cell's own noise, not by how
% many snippets were averaged. A parametric threshold set once for the
% whole dataset would therefore be far too strict for the quiet cells and
% far too loose for the noisy ones; a per-pair null gets this right without
% having to model it.
fprintf('  nSpk range   : %d to %d (%.2fx)\n', ...
    min([rows.nSpk]), max([rows.nSpk]), max([rows.nSpk])/min([rows.nSpk]));
fprintf('  nullSD range : %.2e to %.2e (%.2fx)\n', ...
    min(nsd), max(nsd), max(nsd)/min(nsd));

%% ------------------------------------------------ 5. BH, and why not here
% The p-values above are per pair. The published decision is not: it is
% Benjamini-Hochberg at q = 0.05 over the whole family of directed pairs
% closer than 400 um, which is 9,124 pairs for this session. A subset cannot
% reproduce it, because the BH cutoff depends on the ranks of every p-value
% in the family. That is why section 3 compares against a cutoff measured
% from the full published run rather than recomputing one.
%
% For a real run use runFullSpikeletDetection, which does the whole family.
%
% One implementation note. dataPrepsForMAT.m does not use the standard
% step-up BH; it walks up the sorted p-values and stops at the FIRST one
% that fails p(k) < k/m*q, taking that as the cutoff. Standard BH takes the
% LAST k that passes. The step-down version is the more conservative of the
% two whenever they differ. On this dataset they do not differ at all: both
% reject exactly 4,516 pairs, at cutoffs of 0.024795 and 0.024595. The
% tutorial uses the canonical `benjaminiHochberg` helper and records the
% check, rather than reproducing a non-standard loop that happens to agree.
distMask = data.distMat < cfg.rLim/data.dx & ~logical(eye(data.nCells));
fprintf('\n--- the family ---\n');
fprintf('  directed pairs under %d um : %d\n', cfg.rLim, sum(distMask(:)));
fprintf('  deposited positives        : %d (%.1f%%)\n', ...
    sum(data.deposited.pMat_spklet(:)), 100*sum(data.deposited.pMat_spklet(:))/sum(distMask(:)));
fprintf('  at %.1f s/pair and nRand = 5000, a full run is about %.0f h serial\n', ...
    17, 17*sum(distMask(:))/3600);

fprintf('\nStage 4d done. For the published family, see runFullSpikeletDetection.\n');
