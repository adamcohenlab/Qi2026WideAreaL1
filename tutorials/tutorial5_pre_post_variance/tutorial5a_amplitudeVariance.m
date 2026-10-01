%% Tutorial 5a -- is IPSP amplitude set by the presynaptic or the postsynaptic cell?
%
% Tutorial 3 detected inhibitory connections and measured each one's IPSP
% amplitude. Here those amplitudes are the data. The question is where their
% pair-to-pair variation comes from. If a strongly inhibiting presynaptic
% cell tends to be strong onto every target, amplitude follows the
% presynaptic cell. If some postsynaptic cells receive large IPSPs from
% every input, it follows the postsynaptic cell.
%
% THE MODEL
% For connected pair (i -> j) in animal k,
%
%   log(amplitude_ij) = baseline + animal_k + slope * distance_ij
%                       + a_i + b_j + e_ij
%
%   a_i  one number per presynaptic cell, shared by all of its outputs
%   b_j  one number per postsynaptic cell, shared by all of its inputs
%   e_ij whatever is left, specific to the pair
%
% a, b and e are "random effects": instead of estimating each cell's value
% as a separate parameter, the model estimates how much they vary across
% cells. Those three variances are the answer:
%
%   sigma2_pre   how much amplitude differs between presynaptic cells
%   sigma2_post  how much it differs between postsynaptic cells
%   sigma2_resid how much differs between pairs beyond both
%
% The pre and post effects are "crossed": every cell can be presynaptic in
% some pairs and postsynaptic in others, and the model separates the two
% roles because each cell contributes several pairs in each role.
%
% Why the other terms:
%   log    amplitudes are positive and right-skewed; on a log scale the
%          effects become additive and the scatter roughly normal.
%   animal a fixed offset per mouse, so a difference in overall imaging
%          sensitivity between animals is not counted as cell identity.
%   distance  amplitude may fall with distance, and nearby cells share
%          partners; leaving distance out would let it masquerade as cell
%          identity. Whether distance enters linearly or as log(distance)
%          is chosen from the data (section 5).
%
% Then a parametric bootstrap puts 95% intervals on the three variances and
% tests whether sigma2_pre and sigma2_post differ (section 8).
%
% The manuscript reports this analysis for IPSP amplitude, pooling
% M-YQ0201-27 and M-YQ0201-29. Tutorial 5b does the same for the decay
% time constant.
%
% Runtime: a few seconds to fit; the 2,000-replicate bootstrap takes about
% 6 minutes for both animals.
%
% Source: codex/psp matrix analysis/scripts/
%         runMilestone12VarianceComponentBootstrapTutorial.m

clear; close all

%% 1. Configuration
cfg = struct();
cfg.nwbFiles = { ...   % edit: your local copies
    'path/to/M_YQ0201_27_sparsePulseHad.nwb'
    'path/to/M_YQ0201_29_sparsePulseHad.nwb'};
cfg.labels = {'M27'; 'M29'};   % one short label per session, in the same order
cfg.minDistanceUm = 60;        % connections closer than this are excluded
cfg.nBootstrap = 2000;         % published value; 200 is enough to learn from
cfg.seed = 12012;              % published seed
cfg.saveFile = '';             % e.g. 'tutorial5a_results.mat'; '' saves nothing

%% 2. Load the pairwise metrics
% Only the pairwise table and two unit columns are read, so this takes about
% a second per session even though the files are 25-30 GB.
sessions = cell(numel(cfg.nwbFiles), 1);
for s = 1:numel(cfg.nwbFiles)
    sessions{s} = loadSparsePulseHadForLme(cfg.nwbFiles{s}, ...
        'MinDistanceUm', cfg.minDistanceUm);
end

%% 3. One row per connected pair
% Pairs with a non-positive amplitude cannot be logged and are left out
% (2 of 802 connections in M-YQ0201-27, none in M-YQ0201-29).
[tbl, dRef] = buildConnectedPairTable(sessions, cfg.labels, 'ipspAmpTestN2');

fprintf('\n%-6s %8s %12s %13s %22s\n', 'animal', 'pairs', 'pre cells', ...
    'post cells', 'cells in >=2 pairs pre/post');
for s = 1:numel(cfg.labels)
    rows = tbl.animal == cfg.labels{s};
    nPer = @(g) accumarray(double(removecats(g(rows))), 1);
    fprintf('%-6s %8d %12d %13d %14d / %d\n', cfg.labels{s}, nnz(rows), ...
        numel(nPer(tbl.pre)), numel(nPer(tbl.post)), ...
        nnz(nPer(tbl.pre) >= 2), nnz(nPer(tbl.post) >= 2));
end
% Cells seen in only one pair still count, but only cells seen in several
% pairs tell the model how consistent a cell is, which is what separates
% sigma2_pre and sigma2_post from sigma2_resid.

%% 4. Look at the data
figure('Name', 'IPSP amplitude vs distance', 'Color', 'w');
colors = lines(numel(cfg.labels));
hold on
for s = 1:numel(cfg.labels)
    rows = tbl.animal == cfg.labels{s};
    scatter(tbl.distanceUm(rows), tbl.response(rows), 12, colors(s, :), ...
        'filled', 'MarkerFaceAlpha', 0.5);
end
hold off
set(gca, 'XScale', 'log', 'YScale', 'log', 'Box', 'off', 'TickDir', 'out');
xlabel('Distance between cells (\mum)');
ylabel('IPSP amplitude (spike heights per presynaptic spike)');
legend(cfg.labels, 'Location', 'southwest');
title('Every detected connection');

%% 5. Linear or log distance?
% Both candidates are fitted by maximum likelihood on the same rows and the
% lower AIC wins. Published result: log distance, by 6.2 AIC units.
[formula, distanceComparison] = selectDistanceForm(tbl);
fprintf('\nDistance form, compared by AIC (ML fits):\n');
disp(distanceComparison(:, {'Candidate','AIC','LogLikelihood'}));
fprintf('Selected: %s\n', formula);

%% 6. Fit the chosen model and read off the variances
% REML (restricted maximum likelihood) is the standard way to estimate
% variance components: plain ML underestimates them slightly because it
% ignores that the fixed effects were estimated from the same data.
model = fitlme(tbl, formula, 'FitMethod', 'REML');
point = varianceComponents(model);

fprintf('\nFixed effects:\n');
disp(model.Coefficients(:, {'Name','Estimate','SE','pValue'}));
fprintf(['Variance components of log amplitude (share of the non-fixed ' ...
    'variation):\n']);
fprintf('  presynaptic   %.4f  (%4.1f%%)\n', point.sigmaPre2, 100 * point.fractionPre);
fprintf('  postsynaptic  %.4f  (%4.1f%%)\n', point.sigmaPost2, 100 * point.fractionPost);
fprintf('  residual      %.4f  (%4.1f%%)\n', point.sigmaResid2, 100 * point.fractionResid);

%% 7. What the per-cell effects look like
% The model's best guess of each cell's own offset (a_i or b_j). They are
% pulled toward zero for cells seen in few pairs, so their spread is smaller
% than sigma; the plot is for intuition, and the variances above are the
% estimates. A wider presynaptic histogram means amplitude follows the
% presynaptic cell more than the postsynaptic one.
[blup, blupNames] = randomEffects(model);
isPre = strcmp(blupNames.Group, 'pre');
figure('Name', 'Per-cell effects on log amplitude', 'Color', 'w');
edges = linspace(min(blup), max(blup), 30);
histogram(blup(isPre), edges, 'Normalization', 'probability', ...
    'FaceColor', [0.85 0.4 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on
histogram(blup(~isPre), edges, 'Normalization', 'probability', ...
    'FaceColor', [0.3 0.5 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('Cell''s offset in log amplitude');
ylabel('Fraction of cells');
legend({sprintf('presynaptic (n = %d)', nnz(isPre)), ...
    sprintf('postsynaptic (n = %d)', nnz(~isPre))});
title('Estimated per-cell effects');

%% 8. How certain are the variances? Parametric bootstrap
% Simulate the experiment again from the fitted model, on the same cells and
% pairs, refit, and repeat. The distance form chosen in section 5 is held
% fixed. Progress is printed every 10%.
fprintf('\nBootstrap, %d replicates:\n', cfg.nBootstrap);
boot = bootstrapVarianceComponents(model, tbl, formula, cfg.nBootstrap, cfg.seed);
summary = summarizeVarianceComponentBootstrap(point, boot);

%% 9. Result
fprintf('\n%d of %d refits valid, %.0f s.\n', boot.nValid, boot.nBootstrap, ...
    boot.elapsedSeconds);
disp(summary.components);
fprintf('pre - post = %.4f, 95%% basic interval [%.4f, %.4f]\n', ...
    summary.delta.Estimate, summary.delta.BasicLower95, summary.delta.BasicUpper95);
if summary.nBeyondZero == 0
    comparator = '<=';   % no replicate crossed zero: the p-value is a bound
else
    comparator = '=';
end
fprintf('%d of %d replicates on the other side of zero: two-sided p %s %.4f\n', ...
    summary.nBeyondZero, summary.nValid, comparator, summary.pTwoSided);

plotVarianceComponents(summary, boot, 'IPSP amplitude', ...
    'Variance of log amplitude');

% Reading the result: an interval that excludes zero for the difference
% means the two sources of variation are distinguishable given these cells.
% The intervals assume the model's form (normal random effects on the log
% scale) and treat the chosen distance form as given.

%% 10. Optional save
if ~isempty(cfg.saveFile)
    save(cfg.saveFile, 'cfg', 'tbl', 'dRef', 'formula', 'distanceComparison', ...
        'point', 'boot', 'summary');
    fprintf('Saved %s\n', cfg.saveFile);
end
