% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function pool = fitThetaPoolTwoStage(independentFits, cfg)
%FITTHETAPOOLT WOSTAGE Empirical-Bayes partial pooling of transformed kinetics.
%   Fits a diagonal population covariance to independent q estimates and
%   their local covariance matrices, then returns the closed-form posterior
%   mean and covariance for every pair.  Only q is pooled; all linear terms
%   are re-profiled separately by the caller at qPooled.

if ~iscell(independentFits)
    error('fitThetaPoolTwoStage:Input', 'independentFits must be a cell array of fit structures.');
end
cfg = normaliseCfg(cfg);
nPair = numel(independentFits);
qPooled = nan(3,nPair); VPooled = nan(3,3,nPair);
valid = false(nPair,1);
for p = 1:nPair
    fit = independentFits{p};
    valid(p) = isstruct(fit) && isfield(fit,'qHat') && isfield(fit,'Vq') && ...
        numel(fit.qHat) == 3 && isequal(size(fit.Vq),[3 3]) && ...
        all(isfinite(fit.qHat)) && all(isfinite(fit.Vq(:)));
end
pool = struct('isAvailable',false,'message','','qPooled',qPooled,'VPooled',VPooled, ...
    'mu',nan(3,1),'Sigma',nan(3),'pooledDimensions',[],'validIndependent',valid, ...
    'populationExitFlag',NaN,'populationOutput',[]);
dims = find(cfg.poolOnset * [1 0 0] + [0 1 1]);
pool.pooledDimensions = dims;
validIndex = find(valid);
if numel(validIndex) < cfg.minimumPairs
    pool.message = sprintf('Need at least %d valid independent fits for partial pooling; found %d.', cfg.minimumPairs, numel(validIndex));
    return
end

Q = cell2mat(cellfun(@(x) x.qHat(:), independentFits(validIndex), 'UniformOutput', false));
V = independentFits(validIndex);
mu0 = mean(Q(dims,:), 2);
spread = std(Q(dims,:), 0, 2);
spread = min(cfg.populationStdBounds(2), max(cfg.populationStdBounds(1), spread));
x0 = [mu0; log(spread)];
lb = [-Inf(numel(dims),1); repmat(log(cfg.populationStdBounds(1)),numel(dims),1)];
ub = [ Inf(numel(dims),1); repmat(log(cfg.populationStdBounds(2)),numel(dims),1)];
objective = @(x) populationObjective(x, Q(dims,:), V, dims, cfg);
try
    [xHat, ~, exitFlag, output] = fmincon(objective, x0, [], [], [], [], lb, ub, [], cfg.populationOptions);
catch ME
    pool.message = ME.message;
    return
end
d = numel(dims);
mu = xHat(1:d);
Sigma = diag(exp(2*xHat(d+1:end)));
pool.mu(dims) = mu;
pool.Sigma(dims,dims) = Sigma;
pool.populationExitFlag = exitFlag;
pool.populationOutput = output;

for p = 1:nPair
    if ~valid(p), continue, end
    independent = independentFits{p};
    q = independent.qHat(:); Vp = stabiliseCovariance(independent.Vq(dims,dims), cfg);
    S = Sigma + Vp;
    delta = q(dims) - mu;
    qNew = q;
    qNew(dims) = mu + Sigma * (S \ delta);
    Vnew = independent.Vq;
    Vnew(dims,dims) = Sigma - Sigma * (S \ Sigma);
    qPooled(:,p) = qNew;
    VPooled(:,:,p) = (Vnew + Vnew.') / 2;
end
pool.qPooled = qPooled;
pool.VPooled = VPooled;
pool.isAvailable = true;
pool.message = 'Two-stage diagonal empirical-Bayes population fit.';
end

function value = populationObjective(x, Q, fits, dims, cfg)
d = numel(dims);
mu = x(1:d);
Sigma = diag(exp(2*x(d+1:end)));
value = 0;
for p = 1:size(Q,2)
    Vp = stabiliseCovariance(fits{p}.Vq(dims,dims), cfg);
    S = Sigma + Vp;
    [L, flag] = chol(S, 'lower');
    if flag ~= 0
        value = cfg.invalidPopulationObjective;
        return
    end
    delta = Q(:,p) - mu;
    value = value + 2*sum(log(diag(L))) + sum((L \ delta).^2);
end
end

function V = stabiliseCovariance(V, cfg)
V = (double(V) + double(V).') / 2;
if any(~isfinite(V(:)))
    V = cfg.fallbackVariance * eye(size(V));
    return
end
[U,D] = eig(V);
e = max(diag(D), cfg.minimumObservationVariance);
V = U * diag(e) * U.';
V = (V + V.') / 2;
end

function cfg = normaliseCfg(cfg)
if ~isfield(cfg,'poolOnset'), cfg.poolOnset = true; end
if ~isfield(cfg,'minimumPairs'), cfg.minimumPairs = 3; end
if ~isfield(cfg,'populationStdBounds'), cfg.populationStdBounds = [0.02 5]; end
if ~isfield(cfg,'minimumObservationVariance'), cfg.minimumObservationVariance = 1e-6; end
if ~isfield(cfg,'fallbackVariance'), cfg.fallbackVariance = 1; end
if ~isfield(cfg,'invalidPopulationObjective'), cfg.invalidPopulationObjective = 1e30; end
if ~isfield(cfg,'populationOptions') || isempty(cfg.populationOptions)
    cfg.populationOptions = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',3000);
end
end
