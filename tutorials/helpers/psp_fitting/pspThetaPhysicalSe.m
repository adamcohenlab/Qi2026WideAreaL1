% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function se = pspThetaPhysicalSe(q, Vq, cfg)
%PSPTHETAPHYSICALSE Delta-method standard errors on physical theta scales.

q = q(:);
h = cfg.hessianStep .* max(1, abs(q));
J = zeros(3);
for j = 1:3
    e = zeros(3,1); e(j) = h(j);
    plus = thetaVector(pspThetaQToPhysical(q+e, cfg));
    minus = thetaVector(pspThetaQToPhysical(q-e, cfg));
    J(:,j) = (plus - minus) / (2*h(j));
end
Vtheta = J * Vq * J.';
se = sqrt(max(0, diag(Vtheta)));
end

function v = thetaVector(theta)
v = [theta.latency; theta.tauRise; theta.tauDecay];
end
