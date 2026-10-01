function result = runFullSpikeletDetection(nwbFile, outDir, smokeTest)
%RUNFULLSPIKELETDETECTION Whole-family spikelet detection, the published run.
%
%   result = RUNFULLSPIKELETDETECTION(nwbFile, outDir) corrects every
%   presynaptic row, builds a null for every directed pair closer than
%   400 um, and applies Benjamini-Hochberg across that whole family. This is
%   what Tutorial 4d demonstrates on a handful of pairs.
%
%   result = RUNFULLSPIKELETDETECTION(nwbFile, outDir, true) runs a smoke
%   test: 3 presynaptic cells and 200 nulls, about 2 minutes, exercising
%   every branch. Run this before committing to the full job.
%
%   BH IS WHY THIS CANNOT BE SUBSETTED. The cutoff depends on the ranks of
%   every p-value in the family, so a subset gives different decisions even
%   for pairs it shares with the full run. Detection OUTCOMES can only be
%   compared between complete runs.
%
%   Cost: about 17 s per pair at nRand = 5000, over 9,124 pairs, so roughly
%   43 hours. The loop is SERIAL on purpose. A parfor over pairs has to
%   broadcast the working set — the 4 GB trace matrix and the two 1 GB
%   masks inside `trig` — to every worker, and eight workers exhaust memory
%   and abort the pool. If you want the wall-clock back, split `preList`
%   across separate MATLAB processes, each writing its own output, and
%   concatenate the tables before running BH. BH must see the whole family,
%   so it has to happen after the merge, not inside each process.
%
%   Writes <outDir>/spikelet_detection_full.mat with the per-pair table, the
%   BH outcome, and cfg.
%
%   See also TUTORIAL4D_SPIKELETDETECTION, BENJAMINIHOCHBERG.

if nargin < 3 || isempty(smokeTest)
    smokeTest = false;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

cfg = struct();
cfg.nwbFile   = nwbFile;
cfg.nBack     = 20;
cfg.nFront    = 20;
cfg.nBack_c   = 10;
cfg.nRand     = 5000;      % source literal for crossHadOrthoSpklet_BH.mat
cfg.rLim      = 400;       % um, both the null radius and the BH family radius
cfg.bhQ       = 0.05;
cfg.seed      = 0;
cfg.smokeTest = smokeTest;
if smokeTest
    cfg.nRand = 200;
end

totalTimer = tic;

%% ---------------------------------------------------------------- load
data = loadSparsePulseHadForSpikelet(nwbFile);
trig = buildOrthogonalTriggers(data);
tracesSub = subtractSpikes(data.traces_all, data.spk_t);
data.traces_all = [];
traceMat = tracesSub';          % T x nCells, kept single: 4 GB, not 8
clear tracesSub
nCells = data.nCells;

%% ------------------------------------------------- the family of pairs
distMask = data.distMat < cfg.rLim / data.dx & ~logical(eye(nCells));
preList = unique(find(any(distMask, 2))');
if smokeTest
    preList = preList(round(linspace(1, numel(preList), 3)));
    distMask = distMask & ismember((1:nCells)', preList);
end
[it, jt] = find(distMask);
nPair = numel(it);
fprintf('\nfamily: %d directed pairs under %d um, %d presynaptic cells, nRand = %d\n', ...
    nPair, cfg.rLim, numel(preList), cfg.nRand);

%% ------------------------------- corrected amplitudes, one row at a time
ampCor = nan(nCells);
fprintf('correcting %d presynaptic rows:\n', numel(preList));
rowTimer = tic;
for k = 1:numel(preList)
    ii = preList(k);
    staRow = spikeletStaRow(ii, trig, traceMat, 'nBack', cfg.nBack, 'nFront', cfg.nFront);
    out = correctCommonMode(staRow, data.distMat(ii,:), ...
        'nBack', cfg.nBack, 'dx', data.dx, 'selfIndex', ii);
    ampCor(ii,:) = spikeletAmplitude(out.trRes, cfg.nBack)';
    if mod(k, 25) == 0 || k == numel(preList)
        fprintf('  %3d/%d  (%.0f s)\n', k, numel(preList), toc(rowTimer));
    end
end

%% ----------------------------------------------------- nulls, pair by pair
% Each pair draws its own surrogate sets, seeded from the pair index rather
% than from a shared stream, so a split-and-merge run reproduces a
% single-process run exactly.
pVal = nan(nPair,1);
nSpk = nan(nPair,1);
nCtrl = nan(nPair,1);
ampPair = ampCor(sub2ind([nCells nCells], it, jt));
distMatLocal = data.distMat;
tracesByCell = traceMat.';       % nCells x T view for the null helper

fprintf('null sets for %d pairs:\n', nPair);
pairTimer = tic;
for ip = 1:nPair
    ii = it(ip); jj = jt(ip);
    ns = spikeletNullSet(ii, jj, trig, tracesByCell, distMatLocal, ...
        'nRand', cfg.nRand, 'rLim', cfg.rLim, ...
        'nBack_c', cfg.nBack_c, 'dx', data.dx, 'Seed', cfg.seed + ip);
    pVal(ip) = sum(ampPair(ip) < ns.nullAmp) / (cfg.nRand + 1);
    nSpk(ip) = ns.nSpk;
    nCtrl(ip) = ns.nCtrl;
    if mod(ip, 100) == 0 || ip == nPair
        el = toc(pairTimer);
        fprintf('  %5d/%d  (%.0f s, %.1f s/pair, %.1f h remaining)\n', ...
            ip, nPair, el, el/ip, (nPair-ip)*el/ip/3600);
    end
end
fprintf('  %d pairs in %.0f s (%.1f s/pair)\n', nPair, toc(pairTimer), toc(pairTimer)/nPair);

%% ------------------------------------------------------------------ BH
bh = benjaminiHochberg(pVal, cfg.bhQ);
isDetected = bh.isRejected;

pMat_spklet = false(nCells);
pMat_spklet(sub2ind([nCells nCells], it(isDetected), jt(isDetected))) = true;

%% ------------------------------------------------------------- assemble
result = struct();
result.cfg = cfg;
result.table = table(it, jt, distMatLocal(sub2ind([nCells nCells], it, jt))*data.dx, ...
    nSpk, nCtrl, ampPair, pVal, bh.qValues, isDetected, ...
    'VariableNames', {'preCell','postCell','distanceUm','nSpk','nCtrl', ...
                      'amplitudeCorrected','pValue','qValue','isDetected'});
result.pMat_spklet = pMat_spklet;
result.ampCor = ampCor;
result.bhCutoff = bh.threshold;
result.elapsedSeconds = toc(totalTimer);

%% ------------------------------------------ comparison with the deposit
dep = data.deposited.pMat_spklet;
inFamily = false(nCells); inFamily(sub2ind([nCells nCells], it, jt)) = true;
agree = sum(pMat_spklet(inFamily) == dep(inFamily));
fprintf('\n--- versus the deposit, over the %d pairs in this run ---\n', nPair);
fprintf('  BH cutoff        : %.6f (published full-family value 0.024795)\n', bh.threshold);
fprintf('  detected         : %d (deposit: %d)\n', sum(isDetected), sum(dep(inFamily)));
fprintf('  matching decisions: %d (%.2f%%)\n', agree, 100*agree/nPair);
if smokeTest
    fprintf('  NOTE: smoke test. nRand = %d and the family is a subset, so BH\n', cfg.nRand);
    fprintf('  differs from the published run by construction.\n');
end

outFile = fullfile(outDir, 'spikelet_detection_full.mat');
save(outFile, '-struct', 'result', '-v7.3');
fprintf('\nwrote %s  (%.0f s total)\n', outFile, result.elapsedSeconds);
end
