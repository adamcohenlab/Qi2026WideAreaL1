% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function pValue = empiricalPValue(Tobs, Tnull)
%EMPIRICALPVALUE Add-one empirical one-sided p-value for a test statistic.
validateattributes(Tobs, {'numeric'}, {'scalar','real','finite'}, mfilename, 'Tobs');
validateattributes(Tnull, {'numeric'}, {'vector','real'}, mfilename, 'Tnull');
Tnull = Tnull(isfinite(Tnull));
if isempty(Tnull)
    error('empiricalPValue:NoNullStatistics', 'At least one finite null statistic is required.');
end
pValue = (1 + sum(Tnull >= Tobs)) / (numel(Tnull) + 1);
end
