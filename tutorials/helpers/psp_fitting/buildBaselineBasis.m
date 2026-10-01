% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function G = buildBaselineBasis(t, stimTime, order)
%BUILDBASELINEBASIS Low-order baseline basis for IPSP detection.
%   G = BUILDBASELINEBASIS(T, STIMTIME, ORDER) returns an intercept and a
%   time term centred at STIMTIME and scaled to the supplied time span.
%   ORDER is 1 (default) or 2.  A deliberately low-order basis prevents
%   the baseline model from absorbing slow inhibitory PSPs.

if nargin < 3 || isempty(order)
    order = 1;
end
validateattributes(t, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 't');
validateattributes(stimTime, {'numeric'}, {'scalar','real','finite'}, mfilename, 'stimTime');
validateattributes(order, {'numeric'}, {'scalar','integer','>=',1,'<=',2}, mfilename, 'order');

t = t(:);
span = max(t) - min(t);
if span <= 0
    error('buildBaselineBasis:DegenerateTime', 't must contain at least two distinct time values.');
end
z = (t - stimTime) ./ span;
G = [ones(size(t)), z];
if order == 2
    G = [G, z.^2];
end
end
