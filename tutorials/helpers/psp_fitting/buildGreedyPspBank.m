% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function bank = buildGreedyPspBank(t, stimTime, denseGrid, representativeConditions, cfg)
%BUILDGREEDYPSPBANK Greedily cover a physiological PSP-template family.
%   Coverage is computed after whitening and nuisance projection for each
%   representative condition, under both crosstalk active-set cases.
%   denseGrid is a table or struct with latency, tauRise, and tauDecay.
%   Each representative condition has G, c, and optional lower factor L.

if nargin < 5 || isempty(cfg)
    cfg = struct;
end
cfg = normaliseCfg(cfg);
params = normaliseGrid(denseGrid);
t = t(:);
if isempty(representativeConditions)
    error('buildGreedyPspBank:NoConditions', 'At least one representative condition is required.');
end
if ~isstruct(representativeConditions)
    error('buildGreedyPspBank:Conditions', 'representativeConditions must be a struct array.');
end

M = height(params);
PhiDense = zeros(numel(t), M);
for m = 1:M
    PhiDense(:,m) = makePspTemplate(t, stimTime, params.latency(m), params.tauRise(m), params.tauDecay(m), cfg.measurementKernel);
end

directions = cell(numel(representativeConditions), 2);
valid = true(M, 1);
for r = 1:numel(representativeConditions)
    condn = representativeConditions(r);
    if ~isfield(condn, 'G') || ~isfield(condn, 'c')
        error('buildGreedyPspBank:ConditionFields', 'Each representative condition needs G and c fields.');
    end
    G = condn.G; c = condn.c(:);
    if size(G,1) ~= numel(t) || numel(c) ~= numel(t)
        error('buildGreedyPspBank:ConditionSize', 'Representative G and c must match t.');
    end
    if isfield(condn, 'L') && ~isempty(condn.L)
        L = condn.L;
        GW = L \ G; cW = L \ c; PhiW = L \ PhiDense;
    else
        GW = G; cW = c; PhiW = PhiDense;
    end
    for active = 0:1
        if active == 1
            [Q, ~] = qr([GW, cW], 0);
        else
            [Q, ~] = qr(GW, 0);
        end
        projected = PhiW - Q * (Q' * PhiW);
        norms = sqrt(sum(projected.^2, 1));
        thisValid = norms > cfg.identifiabilityTolerance;
        valid = valid & thisValid(:);
        projected(:, thisValid) = projected(:, thisValid) ./ norms(thisValid);
        projected(:, ~thisValid) = 0;
        directions{r, active + 1} = projected;
    end
end
if ~any(valid)
    error('buildGreedyPspBank:Unidentifiable', 'No dense PSP template is identifiable after nuisance projection.');
end

transformed = [params.latency, log(params.tauRise), log(params.tauDecay - params.tauRise)];
centre = mean(transformed(valid,:), 1);
[~, seedLocal] = min(sum((transformed(valid,:) - centre).^2, 2));
validIdx = find(valid);
selected = validIdx(seedLocal);
coverageHistory = [];
coverage = zeros(M, 1);
while true
    coverage(:) = Inf;
    coverage(valid) = 1;
    for r = 1:size(directions, 1)
        for active = 1:2
            S = directions{r, active};
            matches = (S(:,valid).' * S(:,selected)).^2;
            coverage(valid) = min(coverage(valid), max(matches, [], 2));
        end
    end
    worstCoverage = min(coverage(valid));
    coverageHistory(end+1,1) = worstCoverage; %#ok<AGROW>
    if worstCoverage >= 1 - cfg.delta || numel(selected) >= cfg.maxBankSize
        break
    end
    candidateCoverage = coverage;
    candidateCoverage(~valid) = Inf;
    [~, next] = min(candidateCoverage);
    if any(selected == next)
        break
    end
    selected(end+1,1) = next; %#ok<AGROW>
end

bank.parameters = params(selected,:);
bank.Phi = PhiDense(:,selected);
bank.selectedDenseIndex = selected;
bank.coverageHistory = coverageHistory;
bank.finalCoverage = coverage;
bank.finalWorstCoverage = min(coverage(valid));
bank.identifiableDenseIndex = find(valid);
bank.delta = cfg.delta;
bank.measurementKernel = cfg.measurementKernel;
if bank.finalWorstCoverage < 1 - cfg.delta
    warning('buildGreedyPspBank:MaxBankSize', 'Stopped at maxBankSize before reaching requested coverage.');
end
end

function params = normaliseGrid(denseGrid)
if istable(denseGrid)
    params = denseGrid;
elseif isstruct(denseGrid)
    need = {'latency','tauRise','tauDecay'};
    if ~all(isfield(denseGrid, need))
        error('buildGreedyPspBank:GridFields', 'denseGrid needs latency, tauRise, and tauDecay.');
    end
    params = table(denseGrid.latency(:), denseGrid.tauRise(:), denseGrid.tauDecay(:), 'VariableNames', need);
else
    error('buildGreedyPspBank:GridType', 'denseGrid must be a table or struct.');
end
if height(params) == 0 || any(~isfinite(params{:,1:3}), 'all') || any(params.latency < 0) || any(params.tauRise <= 0) || any(params.tauDecay <= params.tauRise)
    error('buildGreedyPspBank:GridValues', 'Grid requires latency >= 0 and 0 < tauRise < tauDecay.');
end
end

function cfg = normaliseCfg(cfg)
cfg = defaultField(cfg, 'delta', 0.05);
cfg = defaultField(cfg, 'maxBankSize', Inf);
cfg = defaultField(cfg, 'identifiabilityTolerance', 1e-10);
cfg = defaultField(cfg, 'measurementKernel', []);
validateattributes(cfg.delta, {'numeric'}, {'scalar','>',0,'<',1});
validateattributes(cfg.maxBankSize, {'numeric'}, {'scalar','positive'});
end

function s = defaultField(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
