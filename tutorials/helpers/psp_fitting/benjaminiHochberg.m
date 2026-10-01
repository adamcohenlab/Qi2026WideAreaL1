% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function out = benjaminiHochberg(pValues, targetFdr)
%BENJAMINIHOCHBERG Benjamini-Hochberg adjusted q-values and decisions.
if nargin < 2 || isempty(targetFdr)
    targetFdr = 0.05;
end
validateattributes(pValues, {'numeric'}, {'vector','real','>=',0,'<=',1}, mfilename, 'pValues');
validateattributes(targetFdr, {'numeric'}, {'scalar','real','>',0,'<',1}, mfilename, 'targetFdr');

p = pValues(:);
m = numel(p);
[pSorted, order] = sort(p);
rank = (1:m)';
thresholds = targetFdr * rank / m;
last = find(pSorted <= thresholds, 1, 'last');
rejectSorted = false(m, 1);
if ~isempty(last)
    rejectSorted(1:last) = true;
    cutoff = pSorted(last);
else
    cutoff = NaN;
end
qSorted = flipud(cummin(flipud((m ./ rank) .* pSorted)));
qSorted = min(qSorted, 1);
qValues = nan(m, 1); qValues(order) = qSorted;
isRejected = false(m, 1); isRejected(order) = rejectSorted;

out.qValues = qValues;
out.isRejected = isRejected;
out.threshold = cutoff;
out.sortedIndices = order;
out.sortedPValues = pSorted;
end
