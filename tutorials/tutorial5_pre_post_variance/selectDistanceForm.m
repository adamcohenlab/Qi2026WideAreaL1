function [formula, comparison] = selectDistanceForm(tbl)
%SELECTDISTANCEFORM Choose linear or log distance by AIC, fitted by ML.
%
%   [formula, comparison] = SELECTDISTANCEFORM(tbl) fits the crossed
%   random-intercept model twice on the same rows, once with distance entering
%   linearly and once as log(distance), and returns the formula with the lower
%   AIC together with a table of both fits.
%
%   The two candidates differ in their fixed effects, so they are compared
%   with maximum likelihood (ML). REML likelihoods are not comparable between
%   models whose fixed effects differ. The chosen formula is then refitted
%   by REML by the caller, because REML gives less biased variance
%   components, and those are what this tutorial reports.
%
%   With a single session the animal term is dropped, since it would have
%   only one level.

candidateNames = {'linear distance'; 'log distance'};
distanceTerms = {'distanceCentered'; 'logDistanceRatio'};
if numel(categories(removecats(tbl.animal))) > 1
    fixedPart = 'logResponse ~ 1 + animal + ';
else
    fixedPart = 'logResponse ~ 1 + ';
end

% Cell-array inputs, because strcat trims trailing spaces from char inputs.
formulas = strcat({fixedPart}, distanceTerms, {' + (1|pre) + (1|post)'});
aic = nan(2, 1);
logLikelihood = nan(2, 1);
for k = 1:2
    model = fitlme(tbl, formulas{k}, 'FitMethod', 'ML');
    aic(k) = model.ModelCriterion.AIC;
    logLikelihood(k) = model.LogLikelihood;
end
comparison = table(candidateNames, formulas, aic, logLikelihood, ...
    'VariableNames', {'Candidate','Formula','AIC','LogLikelihood'});
[~, best] = min(aic);
formula = formulas{best};
end
