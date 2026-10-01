function summary = summarizeVarianceComponentBootstrap(point, boot)
%SUMMARIZEVARIANCECOMPONENTBOOTSTRAP 95% intervals and a p-value for pre vs post.
%
%   summary = SUMMARIZEVARIANCECOMPONENTBOOTSTRAP(point, boot) takes the
%   point estimates from varianceComponents and the output of
%   bootstrapVarianceComponents, and returns:
%
%     components  table: point estimate and 95% percentile interval for the
%                 presynaptic, postsynaptic and residual variances.
%     delta       point estimate of sigmaPre2 - sigmaPost2 and its 95% basic
%                 interval (see below).
%     pTwoSided   bootstrap p-value for "sigmaPre2 and sigmaPost2 differ".
%     pOneSided   the one-sided version in the direction of the estimate.
%     nBeyondZero replicates whose difference fell on the other side of zero.
%
%   PERCENTILE vs BASIC INTERVAL
%   A variance cannot be negative, so its bootstrap distribution is skewed
%   near zero; the plain 2.5th-97.5th percentile interval respects that
%   boundary. The difference of two variances has no boundary, and there the
%   "basic" interval, which reflects the bootstrap errors around the point
%   estimate, corrects for bias in the bootstrap distribution. The two
%   differ little here.
%
%   THE P-VALUE
%   A parametric bootstrap simulates from the fitted model, not from a model
%   in which pre and post are equal, so the p-value is obtained by inverting
%   the interval: it is the smallest two-sided level at which the percentile
%   interval for the difference would exclude zero. With k replicates on the
%   far side of zero out of n, p = 2 (k + 1) / (n + 1), where the +1 terms
%   count the observed data as one draw. With k = 0 and n = 2000 this gives
%   p = 0.0010, the smallest value 2000 replicates can support; report it as
%   p <= 0.001, never as p = 0.

valid = boot.replicates(boot.replicates.isValid, :);
n = height(valid);
if n < 0.95 * boot.nBootstrap
    warning('Tutorial5:FewValidReplicates', ...
        'Only %d of %d bootstrap refits are valid; treat the intervals with care.', ...
        n, boot.nBootstrap);
end

pointValues = [point.sigmaPre2; point.sigmaPost2; point.sigmaResid2];
bounds = prctile([valid.sigmaPre2, valid.sigmaPost2, valid.sigmaResid2], [2.5 97.5], 1);
summary.components = table({'Presynaptic'; 'Postsynaptic'; 'Residual'}, ...
    pointValues, bounds(1, :)', bounds(2, :)', ...
    'VariableNames', {'Component','Estimate','Lower95','Upper95'});

pointDelta = point.sigmaPre2 - point.sigmaPost2;
errorQuantiles = prctile(valid.deltaVariance - pointDelta, [2.5 97.5]);
summary.delta = table(pointDelta, pointDelta - errorQuantiles(2), ...
    pointDelta - errorQuantiles(1), 'VariableNames', ...
    {'Estimate','BasicLower95','BasicUpper95'});

if pointDelta >= 0
    k = nnz(valid.deltaVariance <= 0);
else
    k = nnz(valid.deltaVariance >= 0);
end
summary.nBeyondZero = k;
summary.nValid = n;
summary.pOneSided = (k + 1) / (n + 1);
summary.pTwoSided = min(1, 2 * (k + 1) / (n + 1));
summary.nBoundaryPre = nnz(valid.sigmaPre2 <= 1e-10);
summary.nBoundaryPost = nnz(valid.sigmaPost2 <= 1e-10);
end
