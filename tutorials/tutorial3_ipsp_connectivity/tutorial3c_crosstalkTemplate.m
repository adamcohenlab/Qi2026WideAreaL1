%% Tutorial 3c -- modelling the blue-light crosstalk
%
% After both corrections in tutorial3b, the pair waveform still contains a large
% positive transient. The corrections removed the field-wide and mid-range
% artifact; what is left is the LOCAL crosstalk, and nearby pairs are exactly
% where scattered light is strongest. So the residual artifact is largest
% precisely where the synaptic signal lives, and it cannot be baselined away.
%
% The alternative is to model it. A pulse aimed at a cell hits that cell's own
% Voltron directly, so the cell's response to its OWN pulse -- the diagonal of
% the tensor -- is essentially pure crosstalk plus its own depolarisation, with
% no synaptic self-input. Its shape, rescaled, is the contaminant in every pair
% where that cell is the RECORDED one.
%
% Tutorial 3d then fits that template jointly with a strictly negative synaptic
% term, with the crosstalk scale constrained non-negative. Modelling rather than
% subtracting is what lets the two be separated.
%
% Source: crosstalkTemplate_nullSet.m.

clear; clc

%% 1. Configuration
cfg = struct;
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
% Pipeline ground truth, used only to validate. Not part of the deposit; set to
% '' to skip the comparison and run from the NWB alone.
cfg.referenceMat = '';   % e.g. 'path/to/ipspWaveform_detect_fit.mat'
cfg.useMat2gray = false;   % true reproduces the source's Image Processing call
cfg.plotCells = [3 17 42 151];

%% 2. Load the tensor diagonal
% Chunks are [320 320 25] in MATLAB order, so a single-pair read pulls a whole
% 320x320x25 slab. Reading the 110 MB tensor once and taking the diagonal is
% much cheaper than 320 strided reads.
data = loadSparsePulseHadForIpsp(cfg.nwbFile, 'LoadTraces', false);
nCells = data.nCells; tBack = data.tBack; tFront = data.tFront;
nTau = tBack + tFront + 1;
tau = (-tBack:tFront) * data.dt;

fprintf('\nReading pair-waveform tensor...\n');
timer = tic;
W = h5read(cfg.nwbFile, data.deposited.wavePath);
selfWaveform = zeros(nCells, nTau);
for ii = 1:nCells, selfWaveform(ii,:) = squeeze(W(ii,ii,:)); end
clear W
fprintf('  %.1f s; %d of %d diagonals are finite\n', toc(timer), ...
    sum(all(isfinite(selfWaveform),2)), nCells);

%% 3. Build the templates
fprintf('\nFitting %d cells...\n', nCells);
timer = tic;
[crsTalkTemp, fitInfo] = buildCrosstalkTemplate(selfWaveform, tBack, tFront, cfg.useMat2gray);
fprintf('  %.1f s; %d cells fitted\n', toc(timer), sum(fitInfo.ok));

fprintf('\nfitted decay rate (1/frames): median %.4f [%.4f, %.4f]\n', ...
    median(fitInfo.tauDecay(fitInfo.ok)), min(fitInfo.tauDecay(fitInfo.ok)), ...
    max(fitInfo.tauDecay(fitInfo.ok)));
fprintf('equivalent decay time constant: median %.1f ms\n', ...
    median(1./fitInfo.tauDecay(fitInfo.ok)) * data.dt * 1e3);
fprintf('\nThe fitted rise rate is computed but intentionally unused; see\n');
fprintf('buildCrosstalkTemplate.m. For comparison, its median is %.4f 1/frames.\n', ...
    median(fitInfo.tauRise(fitInfo.ok)));

%% 4. Check against the pipeline template
if ~isempty(cfg.referenceMat) && exist(cfg.referenceMat, 'file') == 2
    ref = load(cfg.referenceMat, 'crsTalkTemp');
    a = crsTalkTemp; b = double(ref.crsTalkTemp);
    good = all(isfinite(a),2) & all(isfinite(b),2);
    perCell = max(abs(a(good,:) - b(good,:)), [], 2);

    % tauDecay is the only free parameter, recoverable from the first decay
    % sample since that equals exp(-tauDecay).
    tauRef  = -log(b(good, tBack+1+18));
    tauMine = -log(a(good, tBack+1+18));
    tauRel = abs(tauMine - tauRef) ./ tauRef;

    % The metric that actually matters: the detector uses matched-filter
    % projections, so the meaningful error is the squared match lost, not the
    % pointwise difference.
    cs = sum(a(good,:).*b(good,:), 2) ./ ...
         sqrt(sum(a(good,:).^2, 2) .* sum(b(good,:).^2, 2));
    matchLoss = max(1 - cs.^2);

    fprintf('\n--- reproduction of crsTalkTemp ---\n');
    fprintf('cells compared              : %d\n', sum(good));
    fprintf('max abs diff (unit peak)    : %.3e\n', max(perCell));
    fprintf('tauDecay relative diff      : median %.3e, max %.3e\n', ...
        median(tauRel), max(tauRel));
    fprintf('worst squared-match loss    : %.3e\n', matchLoss);
    fprintf('greedy bank tolerance in 3d : %.3e\n', 0.05);

    if matchLoss < 1e-6
        fprintf('\nPASS: recovered from the deposit alone.\n');
        fprintf('The source fitted the UNNORMALISED self waveform; this fitted the\n');
        fprintf('normalised one. Only the fitted rate enters the template, and a rate\n');
        fprintf('is unchanged by a positive rescaling. The residual is nonlinear-solver\n');
        fprintf('tolerance: fit stops on a residual that scales with the data, so a few\n');
        fprintf('cells converge a digit differently. In the matched-filter metric the\n');
        fprintf('detector uses, that is ~6 orders of magnitude below the tolerance the\n');
        fprintf('PSP bank is itself built to in 3d, so it cannot change a detection.\n');
    else
        fprintf('\nMISMATCH -- see TUTORIAL3_VALIDATION.md.\n');
    end
else
    fprintf('\nReference MAT not available; skipping the comparison.\n');
end

%% 5. Optional: template against the waveform it models
if ~isempty(cfg.plotCells)
    plotCrosstalkTemplateDiagnostics(selfWaveform, crsTalkTemp, fitInfo, ...
        tau, tBack, cfg.plotCells);
end

%% 6. What to take away
% Every template is zero until frame +2, rises to a unit peak at +17, then
% decays. Only the decay rate differs between cells. That single degree of
% freedom is what tutorial3d scales by a non-negative alpha, per pair, while
% simultaneously fitting a strictly negative IPSP -- the sign constraint is what
% keeps the two terms from absorbing each other.
%
% Note the template is indexed by the RECORDED cell, not the stimulated one:
% crsTalkTemp(postCell,:).
