% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function theta = pspThetaQToPhysical(q, cfg)
%PSPTHETAQTOPHYSICAL Map unconstrained coordinates to ordered PSP kinetics.
% q = [logit(onset fraction), log(tauRise), log(tauDecay-tauRise-minTauGap)].

q = double(q(:));
if numel(q) ~= 3
    error('pspThetaQToPhysical:Q', 'q must contain exactly three coordinates.');
end
onsetFraction = 1 ./ (1 + exp(-max(-50, min(50, q(1)))));
theta = struct('latency', cfg.onsetBounds(1) + diff(cfg.onsetBounds) * onsetFraction, ...
    'tauRise', exp(q(2)), 'tauDecay', NaN);
theta.tauDecay = theta.tauRise + cfg.minTauGap + exp(q(3));
end
