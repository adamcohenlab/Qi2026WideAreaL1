% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function out = profilePspThetaUnwhitened(q, pairData, cfg)
%PROFILEPSPTHETAUNWHITENED Profile constrained linear terms at continuous PSP kinetics.
%   The model is y = G*beta + alpha*c - A*phi(theta) + error, with
%   0 <= alpha <= alphaMax and A >= 0.  This is intentionally the
%   unwhitened continuation of the detector.  pairData must already use
%   the detector's final mask, selected raw crosstalk delay, and cap.

q = double(q(:));
if numel(q) ~= 3
    error('profilePspThetaUnwhitened:Q', 'q must contain [onset, log(tauRise), log(tauGap)].');
end
required = {'y','G','c','alphaMax','t','stimTime'};
if ~isstruct(pairData) || ~all(isfield(pairData, required))
    error('profilePspThetaUnwhitened:PairData', 'pairData lacks a required fitted-waveform field.');
end
cfg = normaliseCfg(cfg);

theta = pspThetaQToPhysical(q, cfg);
out = emptyOut(theta, q);
if ~isfinite(theta.latency) || theta.latency < cfg.onsetBounds(1) || theta.latency > cfg.onsetBounds(2) || ...
        theta.tauRise < cfg.tauRiseBounds(1) || theta.tauRise > cfg.tauRiseBounds(2) || ...
        theta.tauDecay < cfg.tauDecayBounds(1) || theta.tauDecay > cfg.tauDecayBounds(2) || ...
        theta.tauDecay - theta.tauRise < cfg.minTauGap
    out.rss = cfg.invalidRss;
    return
end

try
    phi = makePspTemplate(pairData.t, pairData.stimTime, theta.latency, theta.tauRise, ...
        theta.tauDecay, cfg.measurementKernel);
    X = [double(pairData.G), double(pairData.c(:)), -double(phi(:))];
    y = double(pairData.y(:));
    nBaseline = size(pairData.G, 2);
    lower = [-Inf(nBaseline,1); 0; 0];
    upper = [ Inf(nBaseline,1); double(pairData.alphaMax); Inf];
    linear = lsqlin(X, y, [], [], [], [], lower, upper, [], cfg.lsqlinOptions);
    residual = y - X * linear;
catch ME
    out.message = ME.message;
    out.rss = cfg.invalidRss;
    return
end

out.rss = residual' * residual;
out.residual = residual;
out.phi = phi;
out.betaHat = linear(1:nBaseline);
out.alphaHat = linear(nBaseline + 1);
out.amplitudeHat = linear(nBaseline + 2);
out.alphaHitCap = isfinite(pairData.alphaMax) && ...
    abs(out.alphaHat - pairData.alphaMax) <= cfg.activeTolerance * max(1, pairData.alphaMax);
out.amplitudeAtZero = out.amplitudeHat <= cfg.activeTolerance;
out.fittedBaseline = pairData.G * out.betaHat;
out.fittedCrosstalk = out.alphaHat * pairData.c(:);
out.fittedPsp = -out.amplitudeHat * phi(:);
out.fittedTotal = out.fittedBaseline + out.fittedCrosstalk + out.fittedPsp;
out.isValid = true;
end

function out = emptyOut(theta, q)
out = struct('q',q,'theta',theta,'rss',Inf,'residual',[],'phi',[], ...
    'betaHat',[],'alphaHat',NaN,'amplitudeHat',NaN,'alphaHitCap',false, ...
    'amplitudeAtZero',false,'fittedBaseline',[],'fittedCrosstalk',[], ...
    'fittedPsp',[],'fittedTotal',[],'isValid',false,'message','');
end

function cfg = normaliseCfg(cfg)
needed = {'onsetBounds','tauRiseBounds','tauDecayBounds','minTauGap','measurementKernel'};
for k = 1:numel(needed)
    if ~isfield(cfg, needed{k})
        error('profilePspThetaUnwhitened:Cfg', 'cfg.%s is required.', needed{k});
    end
end
if ~isfield(cfg, 'invalidRss'), cfg.invalidRss = 1e30; end
if ~isfield(cfg, 'activeTolerance'), cfg.activeTolerance = 1e-8; end
if ~isfield(cfg, 'lsqlinOptions') || isempty(cfg.lsqlinOptions)
    cfg.lsqlinOptions = optimoptions('lsqlin', 'Display', 'off');
end
end
