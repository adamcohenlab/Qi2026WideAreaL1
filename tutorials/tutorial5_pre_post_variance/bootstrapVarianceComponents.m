function boot = bootstrapVarianceComponents(model, tbl, formula, nBootstrap, seed, fitOptions)
%BOOTSTRAPVARIANCECOMPONENTS Parametric bootstrap of a fitted crossed LME.
%
%   boot = BOOTSTRAPVARIANCECOMPONENTS(model, tbl, formula, nBootstrap, seed)
%   repeats the following nBootstrap times:
%
%     1. simulate a new response for every row from the fitted model: new
%        random intercepts for every pre and post cell drawn with the fitted
%        variances, new residual noise, plus the fitted fixed effects;
%     2. refit the SAME formula by REML to the simulated data, keeping every
%        pair, cell and distance exactly as observed;
%     3. record the three variance components.
%
%   The spread of the refitted values shows how much each estimate would
%   wobble if the experiment were repeated on the same cells and pairs.
%   The distance form is not re-selected inside the loop, so the intervals
%   are conditional on that choice.
%
%   fitOptions (optional) is a cell array of extra fitlme name/value pairs,
%   used by tutorial5b to pass the same observation weights to every refit.
%   random() already simulates with those weights (R2019b
%   LinearMixedModel.random divides the draw by sqrt(w)), so simulation and
%   refit agree.
%
%   Failed refits are kept, flagged invalid, and never silently redrawn.
%
%   OUTPUT: struct with a replicate table (sigmaPre2, sigmaPost2,
%   sigmaResid2, deltaVariance = sigmaPre2 - sigmaPost2, isValid, message)
%   plus nValid, seed and elapsed seconds.

if nargin < 6
    fitOptions = {};
end

sigmaPre2 = nan(nBootstrap, 1);
sigmaPost2 = nan(nBootstrap, 1);
sigmaResid2 = nan(nBootstrap, 1);
isValid = false(nBootstrap, 1);
message = repmat({''}, nBootstrap, 1);

rng(seed, 'twister');
timer = tic;
reportEvery = max(1, round(nBootstrap / 10));
for b = 1:nBootstrap
    try
        simulated = tbl;
        simulated.logResponse = random(model);
        refit = fitlme(simulated, formula, 'FitMethod', 'REML', fitOptions{:});
        vc = varianceComponents(refit);
        values = [vc.sigmaPre2, vc.sigmaPost2, vc.sigmaResid2];
        if all(isfinite(values)) && all(values >= 0) && isfinite(refit.LogLikelihood)
            sigmaPre2(b) = vc.sigmaPre2;
            sigmaPost2(b) = vc.sigmaPost2;
            sigmaResid2(b) = vc.sigmaResid2;
            isValid(b) = true;
        else
            message{b} = 'nonfinite or negative variance';
        end
    catch err
        message{b} = err.message;
    end
    if mod(b, reportEvery) == 0
        fprintf('  bootstrap %4d/%d  (%.0f s)\n', b, nBootstrap, toc(timer));
    end
end

boot = struct();
boot.replicates = table(sigmaPre2, sigmaPost2, sigmaResid2, ...
    sigmaPre2 - sigmaPost2, isValid, message, 'VariableNames', ...
    {'sigmaPre2','sigmaPost2','sigmaResid2','deltaVariance','isValid','message'});
boot.nBootstrap = nBootstrap;
boot.nValid = nnz(isValid);
boot.seed = seed;
boot.formula = formula;
boot.elapsedSeconds = toc(timer);
end
