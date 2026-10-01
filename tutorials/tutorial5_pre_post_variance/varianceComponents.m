function vc = varianceComponents(model)
%VARIANCECOMPONENTS Presynaptic, postsynaptic and residual variance of a fit.
%
%   vc = VARIANCECOMPONENTS(model) returns a struct with fields sigmaPre2,
%   sigmaPost2 and sigmaResid2 from a LinearMixedModel containing the random
%   intercepts (1|pre) and (1|post), plus each as a fraction of their sum
%   (fractionPre, fractionPost, fractionResid).
%
%   The groups are found by name, not by position, so the order of the
%   random-effect terms in the formula does not matter.
%
%   The fractions are shares of the variation left AFTER the fixed effects
%   (animal and distance) are accounted for. They answer "of the pair-to-pair
%   scatter not explained by distance or animal, how much follows the
%   presynaptic cell, how much the postsynaptic cell, and how much neither".

[psi, mse, stats] = covarianceParameters(model);
vc = struct('sigmaPre2', NaN, 'sigmaPost2', NaN, 'sigmaResid2', mse);
for k = 1:numel(psi)
    group = strtrim(char(stats{k}.Group(1, :)));
    if strcmpi(group, 'pre')
        vc.sigmaPre2 = psi{k}(1, 1);
    elseif strcmpi(group, 'post')
        vc.sigmaPost2 = psi{k}(1, 1);
    end
end
assert(isfinite(vc.sigmaPre2) && isfinite(vc.sigmaPost2), ...
    'Tutorial5:MissingGroups', 'The model needs random intercepts named pre and post.');

total = vc.sigmaPre2 + vc.sigmaPost2 + vc.sigmaResid2;
vc.fractionPre = vc.sigmaPre2 / total;
vc.fractionPost = vc.sigmaPost2 / total;
vc.fractionResid = vc.sigmaResid2 / total;
end
