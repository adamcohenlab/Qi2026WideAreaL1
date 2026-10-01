% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function phi = makePspTemplate(t, stimTime, latency, tauRise, tauDecay, measurementKernel)
%MAKEPSPTEMPLATE Unit-peak causal double-exponential IPSP template.
%   The returned template is positive; detection uses -A*phi, A >= 0.
%   measurementKernel is optional.  When supplied, it is a sampled causal
%   response/exposure kernel and is convolved using the sampling interval.

if nargin < 6
    measurementKernel = [];
end
validateattributes(t, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 't');
validateattributes(stimTime, {'numeric'}, {'scalar','real','finite'}, mfilename, 'stimTime');
validateattributes(latency, {'numeric'}, {'scalar','real','finite','>=',0}, mfilename, 'latency');
validateattributes(tauRise, {'numeric'}, {'scalar','real','finite','>',0}, mfilename, 'tauRise');
validateattributes(tauDecay, {'numeric'}, {'scalar','real','finite','>',0}, mfilename, 'tauDecay');
if tauRise >= tauDecay
    error('makePspTemplate:InvalidTimeConstants', 'tauRise must be strictly smaller than tauDecay.');
end

t = t(:);
u = t - (stimTime + latency);
phi = zeros(size(t));
isCausal = u >= 0;
phi(isCausal) = exp(-u(isCausal)./tauDecay) - exp(-u(isCausal)./tauRise);

if ~isempty(measurementKernel)
    validateattributes(measurementKernel, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'measurementKernel');
    if numel(t) < 2
        error('makePspTemplate:Sampling', 'At least two time samples are required for convolution.');
    end
    dt = median(diff(t));
    if dt <= 0 || any(abs(diff(t) - dt) > 1e-8 * max(1, abs(dt)))
        error('makePspTemplate:Sampling', 't must be increasing and uniformly sampled when using a measurement kernel.');
    end
    phi = conv(phi, measurementKernel(:), 'same') * dt;
end

peak = max(phi);
if peak <= 0 || ~isfinite(peak)
    error('makePspTemplate:NoPeak', 'The template has no positive peak on the supplied time grid.');
end
phi = phi ./ peak;
end
