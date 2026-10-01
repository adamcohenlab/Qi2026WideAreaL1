%% Tutorial 3d -- detecting a synapse against an empirical null
%
% We now have, per pair, a waveform (3b) and a model of the crosstalk that
% contaminates it (3c). The remaining problem is that the crosstalk is far
% larger than the IPSP and lives in the same time window, so "is there a
% negative deflection?" is not answerable by inspection.
%
% The method turns it into a model-comparison question with two sign
% constraints:
%
%   H0:  y = G*beta + alpha*c              0 <= alpha <= alphaMax
%   H1:  y = G*beta + alpha*c - A*phi_k    same, plus A >= 0
%
% G is a low-order baseline, c the crosstalk template, phi_k a candidate PSP
% shape. Crosstalk may only ever be POSITIVE and the synaptic term only ever
% NEGATIVE. That is what stops the two from absorbing each other: an
% unexplained positive residual cannot be soaked up by -A*phi.
%
% The statistic is how much the best PSP template improves the fit,
% T = max_k (RSS0 - RSS1k). This has no tabulated null distribution, so one is
% built empirically: the same pair's postsynaptic cell, averaged around pulses
% delivered >800 um away. Those carry the same crosstalk but cannot carry a
% synapse from this presynaptic cell. Fitting the identical model to each null
% waveform gives the distribution of T under "no connection", and the p-value
% is the add-one empirical tail.
%
% Sources: tutorial_detect_ipsp_unwhitened.m and the helpers in
% codex\psp fitting. This script calls those helpers unchanged.

clear; clc

%% 1. Configuration
% These are the published values, read back from the cfg stored in
% ipsp_detection_unwhitened_results.mat for both sessions. Do not change them
% if you intend to compare against the deposit.
cfg = struct;
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
% Optional published detection run, used only to compare. Not part of the
% deposit; leave '' to run from the NWB alone.
cfg.referenceMat = '';   % e.g. 'path/to/ipsp_detection_unwhitened_results.mat'

cfg.nPairsDetected    = 5;      % sample this many published-detected pairs
cfg.nPairsNotDetected = 5;      % and this many published-null pairs
cfg.parallelWorkers   = 0;      % 0 = serial; substreams make results identical
cfg.plotDiagnostics   = true;   % per-pair fit decomposition and null histogram
cfg.maxDiagnosticPlots = 2;     % how many detected pairs to show

cfg.nNullTotal   = 10000;
cfg.targetFdr    = 0.05;
cfg.baselineOrder = 2;
cfg.maskFirstSamples = 20;
cfg.rCutNear = 400/6.5;
cfg.rCutFar  = 800/6.5;
cfg.randomSeed = 0;
cfg.bankDelta = 0.05;
cfg.alphaCapStimWindow          = [0.001 0.020];
cfg.crosstalkEarlyWindow        = [0.001 0.020];
cfg.crosstalkNegligibleFraction = 0.05;
cfg.crosstalkCapUpperTailCount  = 20;
cfg.crosstalkCapMinimumPercentile = [];
cfg.latencies = 0.001:0.001:0.01;
cfg.tauRises  = 0.001:0.001:0.008;
cfg.tauDecays = 0.006:0.002:0.1;
cfg.measurementKernel = [];

% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

rng(cfg.randomSeed, 'twister');

%% 2. Load, and rebuild the two inputs from stages 3b/3c
data = loadSparsePulseHadForIpsp(cfg.nwbFile);
nCells = data.nCells; tBack = data.tBack; tFront = data.tFront;
nTime = tBack + tFront + 1;

% The deposited tensor IS the published detection input, differing only by the
% per-pair constant 1/(spkHgtNoBlue(post) * nSpikesPulse(pre)). T scales by the
% square of that constant and the null carries the same factor, so p-values are
% unchanged; amplitudes come out directly in the units of ipspAmpTestN2.
fprintf('\nReading pair-waveform tensor...\n');
pspTensor = double(h5read(cfg.nwbFile, data.deposited.wavePath));

selfWaveform = zeros(nCells, nTime);
for ii = 1:nCells, selfWaveform(ii,:) = pspTensor(ii,ii,:); end
fprintf('Building crosstalk templates (stage 3c)...\n');
crsTalkTemp = buildCrosstalkTemplate(selfWaveform, tBack, tFront);

%% 2b. The null set uses its own, different, trigger definition
%
% READ THIS BEFORE CHANGING ANYTHING HERE.
%
% Everything so far -- nSpikesPulse in 3a, the pair waveform and far-pulse
% baseline in 3b -- used data.t_blue, where tau = 0 is the LAST FRAME BEFORE
% the pulse (convention A). The published detection script does not reuse that.
% It rebuilds the trigger times from scratch, with a definition that lands ONE
% FRAME LATER, and uses those to draw the far-pulse null set:
%
%     tutorial_detect_ipsp_unwhitened.m, lines 95-100
%
% So inside the detector the observed waveform and its own null distribution
% are aligned one frame apart. The crosstalk template, built for convention A,
% therefore fits the nulls slightly worse than it fits the data, which makes
% them marginally easier to beat and mildly inflates significance.
%
% The tutorial reproduces this deliberately. Its job is to document the
% analysis that was actually run, and silently correcting the pipeline would
% produce numbers matching no existing result. But it IS an imperfection, and
% the larger of the two known ones:
%
%     correcting the trigger convention : 138 of ~950 connection calls change
%     correcting the annulus median     :  82 of ~950
%
% and the fix direction is the expected one -- it REMOVES detections, 955 to
% 919, consistent with the published nulls having been too easy to beat.
%
% TO CORRECT IT: pass data.t_blue instead of tBlueNull to generatePairNullSet
% below. That single substitution is the entire fix. The result will then
% disagree with the published pMat by those 138 calls.
%
% Full measurement and reasoning: ISSUE_tblue_trigger_convention.md.
%
% Note the source also reads the blue levels from cell 1 alone and applies no
% sparse-pulse epoch mask. Both are reproduced verbatim here. The mask turns
% out not to matter -- the analog threshold already rejects the Hadamard steps,
% so the two conventions select the same pulses and differ only in alignment.
tBlueNull = cell(nCells,1);
blueLevels = unique(data.blue_all(1,:));
assert(numel(blueLevels) >= 2, 'Tutorial3d:BlueLevels', ...
    'Cell 1 has no detectable blue-pulse level.');
for cellIndex = 1:nCells
    tBlueNull{cellIndex} = find([0, diff(data.blue_all(cellIndex,:))] > blueLevels(2));
end

nA = cellfun(@numel, data.t_blue); nB = cellfun(@numel, tBlueNull);
fprintf(['\ntrigger sets: convention A (waveform) median %d/cell, ' ...
         'convention B (null) median %d/cell\n'], median(nA), median(nB));
fprintf('same pulses, one frame apart -- see the comment above.\n');

%% 3. Time base, mask, baseline basis
t = (-tBack:tFront).' * data.dt;
stimTime = 0;
mask = true(nTime,1);
mask(1:min(cfg.maskFirstSamples, nTime)) = false;
Gfull = buildBaselineBasis(t, stimTime, cfg.baselineOrder);

%% 4. Enumerate pairs EXACTLY as the published run did
% The RNG substream is the pair's index in this list, so the enumeration must
% match or the null sets differ and p-values will not be comparable.
[preAll, postAll] = find(data.distMat <= cfg.rCutNear);
keep = preAll ~= postAll;
preAll = preAll(keep); postAll = postAll(keep);
nPairAll = numel(preAll);
fprintf('\n%d directed pairs within %.0f um.\n', nPairAll, cfg.rCutNear*data.dx);

%% 5. Greedy PSP bank
% A dense (latency, tauRise, tauDecay) grid is reduced to a bank that still
% captures any grid waveform to within bankDelta of its squared match. The
% representative condition is the first pair's crosstalk geometry.
[lg, rg, dg] = ndgrid(cfg.latencies, cfg.tauRises, cfg.tauDecays);
keepGrid = dg > rg;
denseGrid = struct('latency', lg(keepGrid), 'tauRise', rg(keepGrid), 'tauDecay', dg(keepGrid));

yBank = squeeze(pspTensor(preAll(1), postAll(1), :));
cBank = crsTalkTemp(postAll(1),:).';
xtBank = selectCrosstalkDelayAndCap(yBank, cBank, t, true(nTime,1), ~mask, cfg);
representative = struct('G', Gfull(mask,:), 'c', xtBank.cSelected(mask), 'L', []);
fprintf('Building greedy PSP bank from %d dense candidates...\n', numel(denseGrid.latency));
bank = buildGreedyPspBank(t(mask), stimTime, denseGrid, representative, ...
    struct('delta', cfg.bankDelta, 'measurementKernel', cfg.measurementKernel));
fprintf('  %d templates, worst coverage %.6f\n', size(bank.Phi,2), bank.finalWorstCoverage);

%% 6. Choose a subset of pairs to run
% Deterministic: the first N published-detected and first N published-null
% pairs, keeping each pair's GLOBAL index so its null stream is reproduced.
isDetPublished = false(nPairAll,1);
for k = 1:nPairAll
    isDetPublished(k) = data.deposited.pMat(preAll(k), postAll(k));
end
idxDet = find(isDetPublished, cfg.nPairsDetected);
idxNul = find(~isDetPublished, cfg.nPairsNotDetected);
runIdx = sort([idxDet; idxNul]);
nRun = numel(runIdx);
fprintf('\nRunning %d pairs (%d published-detected, %d published-null).\n', ...
    nRun, numel(idxDet), numel(idxNul));

%% 7. Detection
Tobs = nan(nRun,1); pValue = nan(nRun,1); bestTemplateIndex = nan(nRun,1);
amplitudeHat = nan(nRun,1); crosstalkScaleHat = nan(nRun,1);
crosstalkAlphaMax = nan(nRun,1); crosstalkHitCap = false(nRun,1);
crosstalkDelayIndex = nan(nRun,1); rss0 = nan(nRun,1); rss1 = nan(nRun,1);
nNullCal = nan(nRun,1); baselineBetaHat = nan(size(Gfull,2), nRun);
preCell = preAll(runIdx); postCell = postAll(runIdx); pairIndex = runIdx;

% The null statistics are kept so the diagnostic figure in section 11 can show
% the observed statistic against the distribution it was judged by, which is
% the one panel that actually explains the method. The source tutorial keeps
% them for the same reason. They are nNullTotal doubles per pair -- fine for a
% handful of pairs, but drop this if you ever run the full family.
nullStatistics = cell(nRun,1);

% MEMORY, if you raise cfg.parallelWorkers above 0. parfor ships a full copy of
% every variable the loop body names to every worker. Naming `data` therefore
% sends the whole loader struct, including blue_all -- a 3.78 GB matrix the loop
% never touches. Measured on the full-family run: 8.7 GB per worker instead of
% 4.8 GB. So slice out exactly what is needed and do not mention `data` inside
% the loop.
loopTraces       = data.traces_all;
loopDistMat      = data.distMat;
loopRemovedTimes = data.removedTimes;
loopNSpikesPulse = data.nSpikesPulse;

loopTimer = tic;
parfor (k = 1:nRun, cfg.parallelWorkers)
    g = runIdx(k);
    thisPre = preAll(g); thisPost = postAll(g);

    pairStream = RandStream('mrg32k3a', 'Seed', cfg.randomSeed);
    pairStream.Substream = g;

    % Far-pulse null set, in raw normalised-voltage units.
    % tBlueNull, NOT data.t_blue -- one frame later, deliberately. See the
    % block headed "2b. The null set uses its own, different, trigger
    % definition" above before changing this.
    Ynull = generatePairNullSet(thisPre, thisPost, tBlueNull, loopTraces, ...
        loopDistMat, tBack, tFront, cfg.rCutNear, cfg.rCutFar, cfg.nNullTotal, ...
        loopRemovedTimes, pairStream);

    % ...and onto the same scale as the observed waveform. The tensor carries a
    % 1/nSpikesPulse(pre) factor that an STA of the raw trace does not. Omitting
    % this leaves observed and null on different scales and every p-value is
    % wrong. See the loader's NORMALISATION note.
    Ynull = Ynull ./ loopNSpikesPulse(thisPre);

    y = squeeze(pspTensor(thisPre, thisPost, :));
    c = crsTalkTemp(thisPost,:).';

    obsXt  = selectCrosstalkDelayAndCap(y, c, t, true(nTime,1), ~mask, cfg);
    nullXt = selectCrosstalkDelayAndCapMatrix(Ynull, c, t, true(nTime,1), ~mask, cfg);

    obsFit = detectWithPspBank(y(mask), Gfull(mask,:), obsXt.cSelected(mask), [], ...
        bank.Phi, obsXt.alphaMax);
    nullFits = detectWithPspBankMatrix(Ynull(mask,:), Gfull(mask,:), ...
        nullXt.cSelected(mask,:), [], bank.Phi, nullXt.alphaMax);

    Tobs(k) = obsFit.T;
    bestTemplateIndex(k) = obsFit.bestTemplateIndex;
    amplitudeHat(k) = obsFit.amplitudeHat;
    crosstalkScaleHat(k) = obsFit.alphaHat;
    crosstalkAlphaMax(k) = obsFit.alphaMax;
    crosstalkHitCap(k) = obsFit.alphaHitCap;
    crosstalkDelayIndex(k) = obsXt.bestDelayIndex;
    rss0(k) = obsFit.rss0; rss1(k) = obsFit.rss1Best;
    baselineBetaHat(:,k) = obsFit.betaHat;
    nNullCal(k) = numel(nullFits.T);
    nullStatistics{k} = nullFits.T(:);
    pValue(k) = empiricalPValue(obsFit.T, nullFits.T);

    fprintf('[%6.1f s] pair %d/%d  (%d -> %d)  T=%.4g  p=%.4g\n', ...
        toc(loopTimer), k, nRun, thisPre, thisPost, Tobs(k), pValue(k));
end

%% 8. Benjamini-Hochberg
% CAVEAT: BH depends on the whole family of tests. The published run corrected
% across all %d pairs; correcting across this subset gives different q-values,
% so isDetected here is NOT comparable to the deposited pMat. The per-pair
% quantities -- T, amplitude, p-value -- are comparable, and section 10 uses
% those. Raise nPairsDetected/nPairsNotDetected to the full set to compare
% detection outcomes.
qValue = nan(nRun,1); isDetected = false(nRun,1);
valid = isfinite(pValue);
if any(valid)
    bh = benjaminiHochberg(pValue(valid), cfg.targetFdr);
    qValue(valid) = bh.qValues; isDetected(valid) = bh.isRejected;
end

ipspDetectionUnwhitened = table(pairIndex, preCell, postCell, Tobs, ...
    bestTemplateIndex, amplitudeHat, crosstalkScaleHat, crosstalkAlphaMax, ...
    crosstalkHitCap, crosstalkDelayIndex, rss0, rss1, nNullCal, pValue, qValue, ...
    isDetected);
validTemplate = isfinite(bestTemplateIndex);
ipspDetectionUnwhitened.latencyHat  = nan(nRun,1);
ipspDetectionUnwhitened.tauRiseHat  = nan(nRun,1);
ipspDetectionUnwhitened.tauDecayHat = nan(nRun,1);
ipspDetectionUnwhitened.latencyHat(validTemplate)  = bank.parameters.latency(bestTemplateIndex(validTemplate));
ipspDetectionUnwhitened.tauRiseHat(validTemplate)  = bank.parameters.tauRise(bestTemplateIndex(validTemplate));
ipspDetectionUnwhitened.tauDecayHat(validTemplate) = bank.parameters.tauDecay(bestTemplateIndex(validTemplate));
disp(ipspDetectionUnwhitened)

% Kept outside the displayed table: one nNullTotal-long vector per pair, for
% the diagnostic figure in section 11.
ipspDetectionNullStatistics = nullStatistics;

%% 9. Continuous kinetics for detected pairs -- INDEPENDENT fits
% fitDetectedPairThetasUnwhitened also computes a two-stage empirical-Bayes
% partial pooling. Those pooled estimates are NOT used: every published value
% comes from the independent fit (dataPrepsForMAT.m reads only independent*).
% Pooling is left in place because the helper computes it, not because it is
% part of the method.
thetaCfg = struct;
thetaCfg.onsetBounds    = [min(cfg.latencies), 0.02];
thetaCfg.tauRiseBounds  = [min(cfg.tauRises),  0.05];
thetaCfg.tauDecayBounds = [min(cfg.tauDecays), 0.4];
thetaCfg.minTauGap = min(diff(unique(sort([cfg.tauRises(:); cfg.tauDecays(:)]))));
thetaCfg.minTauGap = max(thetaCfg.minTauGap/10, eps);
thetaCfg.measurementKernel = cfg.measurementKernel;
thetaCfg.poolOnset = true;
thetaCfg.nAdditionalStarts = 4;
thetaCfg.minimumPairs = 3;
thetaCfg.progressIntervalSeconds = 10;
thetaCfg.lsqlinOptions = optimoptions('lsqlin','Display','off');
thetaCfg.fminconOptions = optimoptions('fmincon','Display','off', ...
    'Algorithm','interior-point','MaxFunctionEvaluations',1500);

ipspThetaFitsUnwhitened = [];
if any(ipspDetectionUnwhitened.isDetected)
    ipspThetaFitsUnwhitened = fitDetectedPairThetasUnwhitened( ...
        ipspDetectionUnwhitened, pspTensor, crsTalkTemp, t, stimTime, mask, ...
        Gfull(mask,:), thetaCfg);
    disp(ipspThetaFitsUnwhitened.table)
end

%% 10. Check against the published detection run
if ~isempty(cfg.referenceMat) && exist(cfg.referenceMat,'file') == 2
    fprintf('\nLoading published detection table (this is a few GB)...\n');
    R = load(cfg.referenceMat, 'ipspDetectionUnwhitened', 'bank');
    refT = R.ipspDetectionUnwhitened;

    fprintf('\n--- greedy bank ---\n');
    fprintf('templates        : %d here, %d published\n', size(bank.Phi,2), size(R.bank.Phi,2));
    fprintf('worst coverage   : %.6f here, %.6f published\n', ...
        bank.finalWorstCoverage, R.bank.finalWorstCoverage);

    % --- (a) T and amplitude are expected to differ by a per-pair CONSTANT.
    % The published run fitted the unnormalised waveform; this fitted the same
    % waveform divided by s = spkHgtNoBlue(post)*nSpikesPulse(pre). Amplitude is
    % linear in the data and T is quadratic, so if scale is the ONLY difference
    %       T_mine/T_pub == (A_mine/A_pub)^2
    % must hold exactly. That is a sharp test with no free parameters.
    fprintf('\n--- (a) is the difference purely the normalisation constant? ---\n');
    fprintf('%6s %6s %14s %14s %12s\n', 'pre','post','A_mine/A_pub','sqrt(T ratio)','rel. gap');
    scaleGap = nan(nRun,1); impliedS = nan(nRun,1);
    for k = 1:nRun
        r = find(refT.preCell == preCell(k) & refT.postCell == postCell(k), 1);
        if isempty(r), continue; end
        aR = amplitudeHat(k) / refT.amplitudeHat(r);
        tR = sqrt(Tobs(k) / refT.Tobs(r));
        scaleGap(k) = abs(aR - tR) / aR;
        impliedS(k) = 1/aR;
        fprintf('%6d %6d %14.6f %14.6f %12.3e\n', preCell(k), postCell(k), aR, tR, scaleGap(k));
    end
    fprintf('worst disagreement between the two ratios: %.3e\n', max(scaleGap));

    % --- (b) the implied scale must factor as h(post)*n(pre). n(pre) is known,
    % so h(post) = s/n(pre) must come out the SAME for every pair sharing a
    % post cell, even though those pairs have different presynaptic cells.
    fprintf('\n--- (b) does the implied scale factor as h(post)*nSpikesPulse(pre)? ---\n');
    impliedH = impliedS ./ data.nSpikesPulse(preCell);
    for pc = unique(postCell)'
        sel = postCell == pc & isfinite(impliedH);
        if nnz(sel) < 2, continue; end
        v = impliedH(sel);
        fprintf('post cell %3d: %d pairs, implied spkHgtNoBlue = %.4f, spread %.2e\n', ...
            pc, nnz(sel), mean(v), (max(v)-min(v))/mean(v));
    end

    % --- (c) the comparison that needs no reference file: normalised amplitude
    % and decay against the deposit. ipspAmpTestN2 holds the detection amplitude
    % for undetected pairs and the theta-fit amplitude for detected ones.
    fprintf('\n--- (c) normalised amplitude and decay vs the deposit ---\n');
    fprintf('%6s %6s %5s %14s %14s %10s %12s %12s\n', ...
        'pre','post','det','ampHere','ampDeposit','relAmp','tauHere','tauDeposit');
    for k = 1:nRun
        pc = postCell(k); pr = preCell(k);
        depA = data.deposited.ipspAmpTestN2(pr,pc);
        if ipspDetectionUnwhitened.isDetected(k) && ~isempty(ipspThetaFitsUnwhitened)
            q = find(ipspThetaFitsUnwhitened.table.preCell == pr & ...
                     ipspThetaFitsUnwhitened.table.postCell == pc, 1);
            % A is the MAGNITUDE of inhibition: the model is y = ... - A*phi
            % with A >= 0, and the deposit stores it positive. Figure code
            % negates it only for display.
            myA = ipspThetaFitsUnwhitened.table.independentAmplitude(q);
            myTau = ipspThetaFitsUnwhitened.table.independentTauDecay(q);
        else
            myA = amplitudeHat(k); myTau = nan;
        end
        fprintf('%6d %6d %5d %14.6f %14.6f %10.2e %12.5f %12.5f\n', pr, pc, ...
            ipspDetectionUnwhitened.isDetected(k), myA, depA, ...
            abs(myA-depA)/max(abs(depA),eps), myTau, data.deposited.ipspDecay(pr,pc));
    end

    % --- (d) p-values. The tutorial uses the same convention-B nulls as the
    % published run, so these should agree exactly, not merely closely: the
    % mrg32k3a substream is keyed to the pair index, so the same pairs draw the
    % same null events. Any nonzero difference here means the pair enumeration
    % or the trigger set has drifted out of step with the published run.
    fprintf('\n--- (d) p-value and detection agreement with the published run ---\n');
    fprintf('%6s %6s %12s %12s %12s %8s %8s\n', ...
        'pre','post','pHere','pPublished','delta','detHere','detPub');
    dP = nan(nRun,1); agree = 0; nCmp = 0;
    for k = 1:nRun
        r = find(refT.preCell == preCell(k) & refT.postCell == postCell(k), 1);
        if isempty(r), continue; end
        dP(k) = abs(pValue(k) - refT.pValue(r));
        dHere = ipspDetectionUnwhitened.isDetected(k);
        dPub = data.deposited.pMat(preCell(k), postCell(k));
        agree = agree + (dHere == dPub); nCmp = nCmp + 1;
        fprintf('%6d %6d %12.5f %12.5f %12.5f %8d %8d\n', ...
            preCell(k), postCell(k), pValue(k), refT.pValue(r), dP(k), dHere, dPub);
    end
    fprintf('\nmax |p-value| difference : %.3e\n', max(dP));
    fprintf('Monte Carlo SE at p=0.5, nNull=%d : %.3e\n', cfg.nNullTotal, ...
        sqrt(0.25/cfg.nNullTotal));
    fprintf('detection labels agreeing with deposited pMat : %d/%d\n', agree, nCmp);
    fprintf(['\nNote: BH here corrects across %d pairs, the published run across %d,\n' ...
             'so agreement of the labels is informative but not a like-for-like test.\n'], ...
             nRun, nPairAll);
    clear R refT
end

%% 11. Look at what was detected
% Numbers in a table are not the point. For each pair shown, the figure pulls
% the fit apart into baseline, crosstalk and IPSP, shows the fitted IPSP on its
% own with its kinetics, and -- the panel that carries the argument -- puts the
% observed statistic against the empirical null it was judged by.
if cfg.plotDiagnostics
    show = find(ipspDetectionUnwhitened.isDetected);
    if isempty(show)
        fprintf('\nNothing detected in this subset; showing the strongest pair instead.\n');
        [~, show] = max(ipspDetectionUnwhitened.Tobs);
    end
    show = show(1:min(numel(show), cfg.maxDiagnosticPlots));
    for s = show(:)'
        plotIpspDetectionDiagnostics(ipspDetectionUnwhitened, s, pspTensor, ...
            crsTalkTemp, t, mask, Gfull(mask,:), bank, ...
            ipspDetectionNullStatistics{s}, 'Tutorial 3d', ipspThetaFitsUnwhitened);
    end
    % A pair that was NOT detected, for contrast: same machinery, statistic
    % sitting inside the null rather than outside it.
    notDet = find(~ipspDetectionUnwhitened.isDetected, 1);
    if ~isempty(notDet)
        plotIpspDetectionDiagnostics(ipspDetectionUnwhitened, notDet, pspTensor, ...
            crsTalkTemp, t, mask, Gfull(mask,:), bank, ...
            ipspDetectionNullStatistics{notDet}, 'Tutorial 3d (not detected)');
    end
end

%% 12. What to take away
% The empirical null is the load-bearing idea. Nothing about the fitted model is
% assumed to be correct -- not the crosstalk shape, not the PSP family, not the
% noise. The same possibly-wrong model is applied to waveforms that cannot
% contain a synapse, and only an improvement larger than that reference counts.
% Model error common to both cancels.
%
% This is also why the trigger convention had to be unified: if the null
% waveforms are aligned one frame differently from the observed waveform, the
% crosstalk template no longer fits them equally well, the cancellation breaks,
% and the p-values are biased. See ISSUE_tblue_trigger_convention.md.
