% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function fit = fitThetaIndependent(pairData, cfg)
%FITTHETAINDEPENDENT Continuously refine one detected pair's PSP kinetics.
%   Uses the winning detection-bank template as the first start and profiles
%   the detector-constrained linear coefficients at every candidate q.

cfg = normaliseCfg(cfg);
q0 = thetaToQ(pairData.bankTheta, cfg);
starts = makeStarts(q0, cfg);
nStart = size(starts, 2);
startRss = nan(nStart,1); startExitFlag = nan(nStart,1);
startOutput = cell(nStart,1); qBest = q0; bestRss = Inf;

for s = 1:nStart
    objective = @(q) profilePspThetaUnwhitened(q, pairData, cfg).rss;
    try
        [qCandidate, rssCandidate, exitFlag, output] = fmincon(objective, starts(:,s), ...
            [], [], [], [], cfg.qLower, cfg.qUpper, [], cfg.fminconOptions);
        startRss(s) = rssCandidate;
        startExitFlag(s) = exitFlag;
        startOutput{s} = output;
        if isfinite(rssCandidate) && rssCandidate < bestRss
            qBest = qCandidate;
            bestRss = rssCandidate;
        end
    catch ME
        startOutput{s} = struct('message', ME.message);
    end
end

profile = profilePspThetaUnwhitened(qBest, pairData, cfg);
[Vq, hessian, hessianInfo] = localCovariance(qBest, profile.rss, pairData, cfg);
fit = profile;
fit.qHat = qBest(:);
fit.Vq = Vq;
fit.hessian = hessian;
fit.hessianInfo = hessianInfo;
fit.seQ = sqrt(max(0, diag(Vq)));
fit.seTheta = pspThetaPhysicalSe(qBest, Vq, cfg);
fit.startQ = starts;
fit.startRss = startRss;
fit.startExitFlag = startExitFlag;
fit.startOutput = startOutput;
fit.converged = any(startExitFlag > 0) && profile.isValid;
fit.boundaryFlags = boundaryFlags(qBest, profile, cfg);
end

function [Vq, H, info] = localCovariance(q, rss, pairData, cfg)
% RSS has local Hessian 2*J'*J, hence cov(q) ~= 2*sigma2*pinv(H).
h = cfg.hessianStep .* max(1, abs(q(:)));
H = nan(3);
f0 = profilePspThetaUnwhitened(q, pairData, cfg).rss;
for a = 1:3
    ea = zeros(3,1); ea(a) = h(a);
    H(a,a) = (profilePspThetaUnwhitened(q+ea, pairData, cfg).rss - 2*f0 + ...
        profilePspThetaUnwhitened(q-ea, pairData, cfg).rss) / h(a)^2;
    for b = a+1:3
        eb = zeros(3,1); eb(b) = h(b);
        H(a,b) = (profilePspThetaUnwhitened(q+ea+eb, pairData, cfg).rss - ...
            profilePspThetaUnwhitened(q+ea-eb, pairData, cfg).rss - ...
            profilePspThetaUnwhitened(q-ea+eb, pairData, cfg).rss + ...
            profilePspThetaUnwhitened(q-ea-eb, pairData, cfg).rss) / (4*h(a)*h(b));
        H(b,a) = H(a,b);
    end
end
H = (H + H.') / 2;
dof = max(1, numel(pairData.y) - size(pairData.G,2) - 2 - 3);
sigma2 = max(rss / dof, eps);
if any(~isfinite(H(:))) || rcond(H) < cfg.hessianRcond
    Vq = diag(repmat(cfg.fallbackQVariance, 3, 1));
    info = 'fallback: singular or nonfinite local Hessian';
else
    Vq = 2 * sigma2 * pinv(H);
    Vq = (Vq + Vq.') / 2;
    if any(~isfinite(Vq(:))) || any(eig(Vq) <= 0)
        Vq = diag(repmat(cfg.fallbackQVariance, 3, 1));
        info = 'fallback: non-positive local covariance';
    else
        info = 'finite-difference profile-RSS Hessian';
    end
end
end

function starts = makeStarts(q0, cfg)
starts = q0(:);
for s = 1:cfg.nAdditionalStarts
    candidate = q0(:) + cfg.startJitter(:) .* randn(3,1);
    starts(:,end+1) = min(cfg.qUpper, max(cfg.qLower, candidate)); %#ok<AGROW>
end
end

function q = thetaToQ(theta, cfg)
latencyFraction = (theta.latency - cfg.onsetBounds(1)) / diff(cfg.onsetBounds);
latencyFraction = min(1 - 1e-6, max(1e-6, latencyFraction));
q = [log(latencyFraction / (1-latencyFraction)); log(theta.tauRise); ...
    log(max(theta.tauDecay - theta.tauRise - cfg.minTauGap, cfg.minGapExcess))];
q = min(cfg.qUpper, max(cfg.qLower, q));
end

function flags = boundaryFlags(q, profile, cfg)
flags = struct('qLower', any(abs(q - cfg.qLower) <= cfg.boundaryTolerance), ...
    'qUpper', any(abs(q - cfg.qUpper) <= cfg.boundaryTolerance), ...
    'alphaHitCap', profile.alphaHitCap, 'amplitudeAtZero', profile.amplitudeAtZero);
end

function cfg = normaliseCfg(cfg)
required = {'onsetBounds','tauRiseBounds','tauDecayBounds','minTauGap','measurementKernel'};
for k = 1:numel(required)
    if ~isfield(cfg, required{k}), error('fitThetaIndependent:Cfg', 'cfg.%s is required.', required{k}); end
end
if ~isfield(cfg,'minGapExcess'), cfg.minGapExcess = 1e-6; end
if ~isfield(cfg,'qLower')
    cfg.qLower = [-12; log(cfg.tauRiseBounds(1)); log(cfg.minGapExcess)];
end
if ~isfield(cfg,'qUpper')
    cfg.qUpper = [12; log(cfg.tauRiseBounds(2)); log(max(cfg.minGapExcess, cfg.tauDecayBounds(2)-cfg.tauRiseBounds(1)-cfg.minTauGap))];
end
cfg.qLower = cfg.qLower(:); cfg.qUpper = cfg.qUpper(:);
if ~isfield(cfg,'nAdditionalStarts'), cfg.nAdditionalStarts = 4; end
if ~isfield(cfg,'startJitter'), cfg.startJitter = [1;0.35;0.35]; end
if ~isfield(cfg,'hessianStep'), cfg.hessianStep = 1e-3; end
if ~isfield(cfg,'hessianRcond'), cfg.hessianRcond = 1e-10; end
if ~isfield(cfg,'fallbackQVariance'), cfg.fallbackQVariance = 1; end
if ~isfield(cfg,'boundaryTolerance'), cfg.boundaryTolerance = 1e-5; end
if ~isfield(cfg,'fminconOptions') || isempty(cfg.fminconOptions)
    cfg.fminconOptions = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1500);
end
end
