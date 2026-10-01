% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function results = fitDetectedPairThetasUnwhitened(detectionResults, pspTensor, crsTalkTemp, t, stimTime, mask, G, cfg)
%FITDETECTEDPAIRTHETASUNWHITENED Continue detection with two theta analyses.
%   All BH-detected pairs receive an independent continuous fit.  The same
%   complete set then supplies the two-stage theta-only partial-pooling fit.
%   This function deliberately uses no whitening, null waveforms, or new
%   delay/cap selection: it inherits each detected waveform's mask, selected
%   raw crosstalk delay, and alpha cap from the completed detector table.

required = {'isDetected','preCell','postCell','bestTemplateIndex','latencyHat', ...
    'tauRiseHat','tauDecayHat','crosstalkAlphaMax','crosstalkDelayIndex'};
if ~istable(detectionResults) || ~all(ismember(required, detectionResults.Properties.VariableNames))
    error('fitDetectedPairThetasUnwhitened:DetectionResults', 'The unwhitened detection table lacks required fitted-constraint columns.');
end
cfg = normaliseCfg(cfg);
t = t(:); mask = logical(mask(:));
if numel(t) ~= numel(mask) || size(G,1) ~= nnz(mask)
    error('fitDetectedPairThetasUnwhitened:Dimensions', 'G must already be restricted to the detector mask.');
end
detectedRows = find(detectionResults.isDetected & isfinite(detectionResults.bestTemplateIndex));
nPair = numel(detectedRows);
fits = cell(nPair,1); pairData = cell(nPair,1);
fprintf('Unwhitened theta fitting: %d BH-detected pair(s) to fit independently.\n', nPair);
fitTimer = tic;
lastProgressSeconds = 0;
for p = 1:nPair
    row = detectedRows(p);
    pre = detectionResults.preCell(row); post = detectionResults.postCell(row);
    yFull = squeeze(pspTensor(pre, post, :));
    cFull = crsTalkTemp(post,:).';
    if detectionResults.crosstalkDelayIndex(row) == 2
        cFull = [zeros(1,1,'like',cFull); cFull(1:end-1)];
    end
    pairData{p} = struct('pairID', detectionResults.pairIndex(row), 'detectionRow', row, ...
        'preCell',pre,'postCell',post,'y',yFull(mask),'G',G,'c',cFull(mask), ...
        'alphaMax',detectionResults.crosstalkAlphaMax(row),'t',t(mask), ...
        'stimTime',stimTime,'bankTheta',struct('latency',detectionResults.latencyHat(row), ...
        'tauRise',detectionResults.tauRiseHat(row),'tauDecay',detectionResults.tauDecayHat(row)));
    fits{p} = fitThetaIndependent(pairData{p}, cfg);
    elapsedSeconds = toc(fitTimer);
    if p == nPair || elapsedSeconds - lastProgressSeconds >= cfg.progressIntervalSeconds
        fprintf('[%6.1f s] Unwhitened theta: independent fit %d/%d (%d -> %d), RSS %.4g, converged %d.\n', ...
            elapsedSeconds, p, nPair, pre, post, fits{p}.rss, fits{p}.converged);
        lastProgressSeconds = elapsedSeconds;
    end
end

fprintf('Unwhitened theta fitting: estimating theta-only partial-pooling population model.\n');
pool = fitThetaPoolTwoStage(fits, cfg);
pooledFits = cell(nPair,1);
if pool.isAvailable
    poolTimer = tic;
    lastProgressSeconds = 0;
    for p = 1:nPair
        if all(isfinite(pool.qPooled(:,p)))
            pooledFits{p} = profilePspThetaUnwhitened(pool.qPooled(:,p), pairData{p}, cfg);
            pooledFits{p}.qHat = pool.qPooled(:,p);
            pooledFits{p}.Vq = pool.VPooled(:,:,p);
            pooledFits{p}.seTheta = pspThetaPhysicalSe(pool.qPooled(:,p), pool.VPooled(:,:,p), cfg);
        end
        elapsedSeconds = toc(poolTimer);
        if p == nPair || elapsedSeconds - lastProgressSeconds >= cfg.progressIntervalSeconds
            pair = pairData{p};
            fprintf('[%6.1f s] Unwhitened theta: pooled re-profile %d/%d (%d -> %d).\n', ...
                elapsedSeconds, p, nPair, pair.preCell, pair.postCell);
            lastProgressSeconds = elapsedSeconds;
        end
    end
else
    fprintf('Unwhitened theta fitting: partial pooling unavailable: %s\n', pool.message);
end
results = struct('detectedRows',detectedRows,'pairData',{pairData},'independentFits',{fits}, ...
    'pool',pool,'pooledFits',{pooledFits},'cfg',cfg);
results.table = makeSummaryTable(pairData, fits, pooledFits, pool);
end

function tableOut = makeSummaryTable(pairData, independentFits, pooledFits, pool)
n = numel(pairData);
pairID = nan(n,1); preCell = nan(n,1); postCell = nan(n,1);
indTheta = nan(n,3); pooledTheta = nan(n,3); indSe = nan(n,3); pooledSe = nan(n,3); indRss = nan(n,1); pooledRss = nan(n,1);
indAmp = nan(n,1); pooledAmp = nan(n,1); indAlpha = nan(n,1); pooledAlpha = nan(n,1);
shrinkageDistance = nan(n,1); independentConverged = false(n,1); alphaHitCap = false(n,1); amplitudeAtZero = false(n,1);
for p = 1:n
    d = pairData{p}; f = independentFits{p};
    pairID(p) = d.pairID; preCell(p) = d.preCell; postCell(p) = d.postCell;
    indTheta(p,:) = thetaRow(f.theta); indSe(p,:) = f.seTheta(:).'; indRss(p) = f.rss;
    indAmp(p) = f.amplitudeHat; indAlpha(p) = f.alphaHat; independentConverged(p) = f.converged;
    alphaHitCap(p) = f.alphaHitCap; amplitudeAtZero(p) = f.amplitudeAtZero;
    if pool.isAvailable && ~isempty(pooledFits{p})
        pf = pooledFits{p}; pooledTheta(p,:) = thetaRow(pf.theta); pooledSe(p,:) = pf.seTheta(:).'; pooledRss(p) = pf.rss;
        pooledAmp(p) = pf.amplitudeHat; pooledAlpha(p) = pf.alphaHat;
        shrinkageDistance(p) = norm(pf.qHat - f.qHat);
    end
end
tableOut = table(pairID,preCell,postCell,indTheta(:,1),indTheta(:,2),indTheta(:,3), ...
    pooledTheta(:,1),pooledTheta(:,2),pooledTheta(:,3),indSe(:,1),indSe(:,2),indSe(:,3),pooledSe(:,1),pooledSe(:,2),pooledSe(:,3), ...
    indRss,pooledRss,indAmp,pooledAmp,indAlpha,pooledAlpha,shrinkageDistance, ...
    independentConverged,alphaHitCap,amplitudeAtZero, ...
    'VariableNames',{'pairID','preCell','postCell','independentLatency','independentTauRise','independentTauDecay', ...
    'pooledLatency','pooledTauRise','pooledTauDecay','independentLatencySE','independentTauRiseSE','independentTauDecaySE','pooledLatencySE','pooledTauRiseSE','pooledTauDecaySE', ...
    'independentRSS','pooledRSS','independentAmplitude','pooledAmplitude','independentCrosstalk','pooledCrosstalk', ...
    'shrinkageDistance','independentConverged','independentAlphaHitCap','independentAmplitudeAtZero'});
end

function row = thetaRow(theta)
row = [theta.latency theta.tauRise theta.tauDecay];
end

function cfg = normaliseCfg(cfg)
if ~isfield(cfg,'onsetBounds'), cfg.onsetBounds = [0.001 0.010]; end
if ~isfield(cfg,'tauRiseBounds'), cfg.tauRiseBounds = [0.001 0.008]; end
if ~isfield(cfg,'tauDecayBounds'), cfg.tauDecayBounds = [0.006 0.100]; end
if ~isfield(cfg,'minTauGap'), cfg.minTauGap = 1e-5; end
if ~isfield(cfg,'measurementKernel'), cfg.measurementKernel = []; end
if ~isfield(cfg,'hessianStep'), cfg.hessianStep = 1e-3; end
if ~isfield(cfg,'progressIntervalSeconds'), cfg.progressIntervalSeconds = 10; end
if ~isfield(cfg,'lsqlinOptions'), cfg.lsqlinOptions = optimoptions('lsqlin','Display','off'); end
if ~isfield(cfg,'fminconOptions'), cfg.fminconOptions = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',1500); end
end
