function [crsTalkTemp, fitInfo] = buildCrosstalkTemplate(selfWaveform, tBack, tFront, useMat2gray)
%BUILDCROSSTALKTEMPLATE Per-cell optical-crosstalk template from the self waveform.
%
%   [crsTalkTemp, fitInfo] = BUILDCROSSTALKTEMPLATE(selfWaveform, tBack, tFront)
%
%   selfWaveform is nCells-by-nTau: each cell's own response to its own blue
%   pulse, i.e. the diagonal of the pair-waveform tensor. crsTalkTemp is
%   nCells-by-nTau with unit positive peak, one template per cell.
%
%   WHAT THIS IS. A pulse aimed at a cell drives that cell optogenetically and
%   also hits its Voltron directly. The diagonal waveform is therefore almost
%   pure crosstalk plus the cell's own depolarisation, with no synaptic input
%   from itself. Its shape, rescaled, is the artifact that contaminates every
%   pair involving that cell as the RECORDED cell -- which is why the detector
%   in 3d indexes the template by the postsynaptic cell.
%
%   SHAPE. Two exponentials are fitted to the self waveform: a rise over frames
%   +2..+17 and a decay over +18..+68. The template is then built as a
%   normalised rise followed by a pure decay, zero before frame +2.
%
%   WHY THE RISE USES THE DECAY CONSTANT. The template's rising phase is built
%   from tau_decay, not from the fitted tau_rise, which is computed and
%   discarded. This is deliberate, not an oversight: the cell is spiking hard
%   during the pulse, and those spikes corrupt the rising phase badly enough
%   that its fitted time constant is not trustworthy. The decay is measured
%   after the pulse ends, when spiking has stopped, so it is the reliable one
%   and is reused to shape the rise. This is the same constraint that forces the
%   spectral spike-count estimator in tutorial3a: nothing spike-shaped is
%   measurable during the pulse.
%
%   WHY FRAME +2. Under the convention-A trigger used throughout, tau = 0 is the
%   last frame BEFORE the pulse, so the blue is not commanded on until +1, and
%   the first camera frame materially affected is +2. Measured, not assumed --
%   see ISSUE_tblue_trigger_convention.md.
%
%   SCALE INVARIANCE. Only the fitted rate a3 enters the template, and a3 is
%   unchanged when the data are multiplied by a positive constant (the start
%   point scales the same way). So this returns the same template from the
%   deposited normalised waveform as from the source unnormalised one.
%
%   useMat2gray (default false) uses the Image Processing Toolbox mat2gray, as
%   the source does. The default expands it inline as (x-min)/(max-min), which
%   is what mat2gray computes for this monotone input, and drops the dependency.
%
%   Source: crosstalkTemplate_nullSet.m. Requires Curve Fitting Toolbox.

if nargin < 4, useMat2gray = false; end

nCells = size(selfWaveform,1);
nTau = tBack + tFront + 1;
assert(size(selfWaveform,2) == nTau, 'selfWaveform must be nCells-by-(tBack+tFront+1).');

ftRise  = fittype('a1+a2*(1-exp(-a3*x))');
ftDecay = fittype('a1+a2*exp(-a3*x)');

xRise  = (2:17)';
xDecay = (0:50)';

crsTalkTemp = zeros(nCells, nTau);
fitInfo = struct('tauRise', nan(nCells,1), 'tauDecay', nan(nCells,1), ...
                 'ok', false(nCells,1));

for ii = 1:nCells
    yRise  = double(selfWaveform(ii, tBack+1+xRise)).';
    yDecay = double(selfWaveform(ii, tBack+1+18+xDecay)).';

    if any(~isfinite(yRise)) || any(~isfinite(yDecay))
        crsTalkTemp(ii,:) = nan;      % cell has no usable self waveform
        continue
    end

    % Start points exactly as in the source. tau0 is the half-maximum crossing,
    % which is scale-free, so the whole start point scales with the data.
    [~, tau0] = min(abs(yRise - max(yRise)/2));
    fRise = fit(xRise, yRise, ftRise, 'startPoint', [yRise(1) max(yRise) 1/tau0]);

    [~, tau0] = min(abs(yDecay - max(yDecay)/2));
    fDecay = fit(xDecay, yDecay, ftDecay, ...
        'startPoint', [mean(yDecay(end-5:end)) max(yDecay) 1/tau0]);

    c = coeffvalues(fRise);  fitInfo.tauRise(ii)  = c(3);   % kept, not used below
    c = coeffvalues(fDecay); tauDecay = c(3);
    fitInfo.tauDecay(ii) = tauDecay;
    fitInfo.ok(ii) = true;

    riseShape = 1 - exp(-tauDecay*(0:15));
    if useMat2gray
        riseShape = mat2gray(riseShape);
    else
        riseShape = (riseShape - min(riseShape)) ./ (max(riseShape) - min(riseShape));
    end

    crsTalkTemp(ii, tBack+1+(2:17))  = riseShape;
    crsTalkTemp(ii, tBack+1+18:end)  = exp(-tauDecay*(1:tFront-17));
end
end
