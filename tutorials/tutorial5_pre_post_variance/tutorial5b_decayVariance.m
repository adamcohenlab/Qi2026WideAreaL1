%% Tutorial 5b -- is IPSP decay set by the presynaptic or the postsynaptic cell?
%
% The same question as tutorial 5a, asked of the IPSP decay time constant
% (tau) instead of the amplitude, with the same crossed model:
%
%   log(tau_ij) = baseline + animal_k + slope * distance_ij + a_i + b_j + e_ij
%
% Read tutorial 5a first; this script only explains what is new.
%
% WHAT IS NEW: PRECISION WEIGHTS
% Every tau comes out of a curve fit and carries a standard error. Some
% decays are measured precisely, others are barely constrained. In an
% unweighted model the noisy ones add scatter that lands in the residual
% and can blur the pre/post comparison. Weighting lets precise measurements
% count more.
%
% Weighting by precision alone (1 / SE^2) would go too far: a handful of
% extremely precise pairs would dominate, although real pairs differ from
% one another by more than their measurement error. So each pair's weight is
%
%   w_ij = 1 / (SE_ij^2 + s0^2)
%
% where SE_ij is the standard error of log(tau_ij) and s0^2 is the genuine
% pair-to-pair scatter that remains however precisely a pair is measured.
% s0^2 is estimated once from the data ("method of moments"): the residual
% variance of the unweighted model minus the typical measurement variance.
% Once s0^2 is large relative to SE^2, the weights flatten out, and pairs
% differ in influence only where measurement noise is a real part of their
% error.
%
% The weights are rescaled to average 1, so the residual variance keeps the
% meaning "scatter of a typical pair". The model then assumes a pair's
% residual variance is sigma2_resid / w_ij.
%
% The manuscript reports this weighted model for decay. The unweighted
% estimates are printed alongside for comparison (section 6).
%
% Runtime: a few seconds to fit; the 2,000-replicate bootstrap takes about
% 7 minutes for both animals.
%
% Source: codex/psp matrix analysis/scripts/
%         runMilestone12MomentWeightedDecayBootstrapTutorial.m

clear; close all

%% 1. Configuration
cfg = struct();
cfg.nwbFiles = { ...   % edit: your local copies
    'path/to/M_YQ0201_27_sparsePulseHad.nwb'
    'path/to/M_YQ0201_29_sparsePulseHad.nwb'};
cfg.labels = {'M27'; 'M29'};
cfg.minDistanceUm = 60;
cfg.nBootstrap = 2000;         % published value
cfg.seed = 12014;              % published seed
cfg.saveFile = '';

%% 2. Load the pairwise metrics
sessions = cell(numel(cfg.nwbFiles), 1);
for s = 1:numel(cfg.nwbFiles)
    sessions{s} = loadSparsePulseHadForLme(cfg.nwbFiles{s}, ...
        'MinDistanceUm', cfg.minDistanceUm);
end

%% 3. One row per connected pair, with its measurement error
% The deposit does not record tau's time unit, but the fit that produced it
% (tutorial3d, tauDecayBounds = [0.006 0.4]) works in seconds. It would not
% matter anyway: the model uses log(tau), a unit change only shifts the
% baseline, and SE/tau below is unit-free.
[tbl, dRef] = buildConnectedPairTable(sessions, cfg.labels, 'ipspDecay');

% Standard error of log(tau) from the standard error of tau (delta method:
% the derivative of log(tau) is 1/tau).
tbl.seLogTau = tbl.seTau ./ tbl.response;
tbl.measurementVariance = tbl.seLogTau .^ 2;

usable = isfinite(tbl.logResponse) & tbl.response > 0 & ...
    isfinite(tbl.seTau) & tbl.seTau > 0 & ...
    isfinite(tbl.measurementVariance) & tbl.measurementVariance > 0;
fprintf('\n%d of %d connected pairs have a usable decay standard error.\n', ...
    nnz(usable), height(tbl));
tbl = tbl(usable, :);

%% 4. Linear or log distance?
% Chosen exactly as in 5a, from unweighted ML fits on the rows used below.
% Published result: log distance, by 9.2 AIC units.
[formula, distanceComparison] = selectDistanceForm(tbl);
disp(distanceComparison(:, {'Candidate','AIC','LogLikelihood'}));
fprintf('Selected: %s\n', formula);

%% 5. Build the weights
unweightedModel = fitlme(tbl, formula, 'FitMethod', 'REML');
unweighted = varianceComponents(unweightedModel);

typicalMeasurementVariance = median(tbl.measurementVariance);
s0Squared = max(0, unweighted.sigmaResid2 - typicalMeasurementVariance);
rawWeight = 1 ./ (tbl.measurementVariance + s0Squared);
tbl.weight = rawWeight ./ mean(rawWeight);

effectiveN = sum(tbl.weight)^2 / sum(tbl.weight.^2);
fprintf(['\nUnweighted residual variance %.3f, median measurement variance %.3f, ' ...
    'so s0^2 = %.3f.\n'], unweighted.sigmaResid2, typicalMeasurementVariance, s0Squared);
fprintf('Weights range %.3g to %.3g (mean 1). Effective number of pairs %.0f of %d.\n', ...
    min(tbl.weight), max(tbl.weight), effectiveN, height(tbl));
% In the published data s0^2 (0.66) is more than ten times the median
% measurement variance (0.057), so most pairs weigh about the same and the
% weighting mainly discounts the few very noisy fits (effective n 1,364 of
% 1,469).

figure('Name', 'Decay precision weights', 'Color', 'w', 'Position', [100 100 800 330]);
subplot(1, 2, 1);
scatter(tbl.seLogTau, tbl.weight, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
set(gca, 'XScale', 'log', 'Box', 'off', 'TickDir', 'out');
xlabel('Standard error of log \tau');
ylabel('Weight');
title('Noisy fits are discounted; precise ones plateau');
subplot(1, 2, 2);
histogram(tbl.weight, 40, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('Weight');
ylabel('Pairs');
title(sprintf('Effective n = %.0f of %d', effectiveN, height(tbl)));

%% 6. Fit the weighted model
model = fitlme(tbl, formula, 'FitMethod', 'REML', 'Weights', tbl.weight);
point = varianceComponents(model);

% Decay fits that stopped at a bound. Tutorial 3d fits tau between 6 ms and
% 400 ms. A pair whose tau sits exactly on either limit was not measured, only
% bounded, and its standard error is large (SE of log tau about 1 against
% about 0.2 for the rest). The weights discount these pairs to about 0.46
% instead of about 1.1, which is the spike near 0.5 in the weight histogram.
atBound = tbl.response >= 0.4 - 1e-9 | tbl.response <= 0.006 + 1e-9;
fprintf('\n%d of %d pairs have tau at a fit bound (%d at 400 ms, %d at 6 ms).\n', ...
    nnz(atBound), height(tbl), nnz(tbl.response >= 0.4 - 1e-9), ...
    nnz(tbl.response <= 0.006 + 1e-9));

fprintf('\nFixed effects (weighted):\n');
disp(model.Coefficients(:, {'Name','Estimate','SE','pValue'}));
fprintf('Variance components of log tau:   weighted          unweighted\n');
fprintf('  presynaptic               %.4f (%4.1f%%)    %.4f (%4.1f%%)\n', ...
    point.sigmaPre2, 100 * point.fractionPre, unweighted.sigmaPre2, 100 * unweighted.fractionPre);
fprintf('  postsynaptic              %.4f (%4.1f%%)    %.4f (%4.1f%%)\n', ...
    point.sigmaPost2, 100 * point.fractionPost, unweighted.sigmaPost2, 100 * unweighted.fractionPost);
fprintf('  residual                  %.4f (%4.1f%%)    %.4f (%4.1f%%)\n', ...
    point.sigmaResid2, 100 * point.fractionResid, unweighted.sigmaResid2, 100 * unweighted.fractionResid);

%% 6b. Why the weighting changes the presynaptic variance so much
% Unweighted, presynaptic and postsynaptic variance come out much closer
% together (pooled 0.107 vs 0.163); weighted, presynaptic variance nearly
% vanishes (0.021 vs 0.111). The check below shows that
% the bound-limited fits flagged above are the reason. Drop them, and even
% the unweighted model puts presynaptic variance well below postsynaptic,
% while the weighted estimates barely move. So the weights are not creating
% the result: they do, smoothly, roughly what removing the failed fits does
% outright. Point estimates only; the bootstrap below uses all pairs.
kept = tbl(~atBound, :);
keptUnweighted = varianceComponents(fitlme(kept, formula, 'FitMethod', 'REML'));
keptS0 = max(0, keptUnweighted.sigmaResid2 - median(kept.measurementVariance));
keptWeight = 1 ./ (kept.measurementVariance + keptS0);
keptWeighted = varianceComponents(fitlme(kept, formula, 'FitMethod', 'REML', ...
    'Weights', keptWeight ./ mean(keptWeight)));
fprintf('\nWithout the %d bound-limited pairs:  weighted    unweighted\n', nnz(atBound));
fprintf('  presynaptic                        %.4f      %.4f\n', ...
    keptWeighted.sigmaPre2, keptUnweighted.sigmaPre2);
fprintf('  postsynaptic                       %.4f      %.4f\n', ...
    keptWeighted.sigmaPost2, keptUnweighted.sigmaPost2);
fprintf('  residual                           %.4f      %.4f\n', ...
    keptWeighted.sigmaResid2, keptUnweighted.sigmaResid2);

%% 7. Parametric bootstrap, with the same weights in every refit
% The weights are computed once, above, and held fixed. They are not
% rebuilt from each simulated dataset.
fprintf('\nBootstrap, %d replicates:\n', cfg.nBootstrap);
boot = bootstrapVarianceComponents(model, tbl, formula, cfg.nBootstrap, cfg.seed, ...
    {'Weights', tbl.weight});
summary = summarizeVarianceComponentBootstrap(point, boot);

%% 8. Result
fprintf('\n%d of %d refits valid, %.0f s.\n', boot.nValid, boot.nBootstrap, ...
    boot.elapsedSeconds);
disp(summary.components);
fprintf('pre - post = %.4f, 95%% basic interval [%.4f, %.4f]\n', ...
    summary.delta.Estimate, summary.delta.BasicLower95, summary.delta.BasicUpper95);
if summary.nBeyondZero == 0
    comparator = '<=';
else
    comparator = '=';
end
fprintf('%d of %d replicates on the other side of zero: two-sided p %s %.4f\n', ...
    summary.nBeyondZero, summary.nValid, comparator, summary.pTwoSided);
fprintf('Presynaptic variance estimated as zero in %d of %d replicates.\n', ...
    summary.nBoundaryPre, summary.nValid);

plotVarianceComponents(summary, boot, 'IPSP decay time constant (weighted)', ...
    'Variance of log \tau');

% Reading the result: for decay the ordering is the reverse of amplitude.
% The postsynaptic cell accounts for most of the non-fixed variation, and
% the presynaptic variance is small enough that some bootstrap refits put it
% at exactly zero, which is why its interval starts at 0. The honest
% summary is that presynaptic identity explains little or none of the
% variation in decay, not that it explains a small but definite amount.

%% 9. Optional save
if ~isempty(cfg.saveFile)
    save(cfg.saveFile, 'cfg', 'tbl', 'dRef', 'formula', 'distanceComparison', ...
        'unweighted', 's0Squared', 'point', 'boot', 'summary');
    fprintf('Saved %s\n', cfg.saveFile);
end
