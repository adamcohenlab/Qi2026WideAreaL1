function selection = selectMotionFreeFrames(motionX, motionY, options)
% SELECTMOTIONFREEFRAMES  Frames judged free of large motion, from the flow field.
%
% selection = SELECTMOTIONFREEFRAMES(motionX, motionY, options)
%
% Registration and a frozen spatial footprint do not remove all motion artifact,
% so every later stage of the pipeline estimates its thresholds and templates only
% from frames where the tissue was reasonably still. This computes that frame set
% from the per-cell displacement traces the optical-flow registration produced.
%
% The statistic
% -------------
% Average the displacement across cells, remove line noise, high-pass with a
% running median, and take the product of the absolute x and y residuals. Using the
% product rather than the norm means a frame counts as moving only when it moves in
% both axes, which suppresses single-axis noise in the flow estimate.
%
% Choosing the cut
% ----------------
% The published pipeline uses a multiple of the statistic's mean, chosen by eye
% from plots of the motion trace. It uses a *different* multiple at each stage:
%
%     4   demixing            large.voltage_FOV_NMFDemix_v8
%     2   spike finding       large_voltage_FOV_MF_negThresh_behavior
%     6   spike-height fit    get_reliable_spikeHeight_behavior_v2
%
% So `multiplier` has no single correct value; pass the one belonging to the stage
% being reproduced. The default here is 2, the spike-finding value.
%
% A more robust alternative, which the author has endorsed but which is NOT the
% published behaviour, is to model the statistic's distribution as a mixture: still
% and moving frames should come from components with different means, by analogy
% with the whisking-motion mixture fit used elsewhere in the manuscript. Set
% options.reportMixtureCut to compute that cut alongside the published one, for
% comparison. It is reported, never applied -- substituting it would change
% published results.
%
% Inputs
% ------
% motionX, motionY : [nCells x nSamples] or [1 x nSamples] double
%     Per-cell x and y displacement traces. If a cell axis is present it is
%     averaged, matching the source pipeline.
%
% options : struct, all optional
%     .multiplier        cut = multiplier * mean(statistic). Default 2.
%     .medianWindow      high-pass running-median window, samples. Default 400.
%     .tailGuard         samples excluded at the end of the recording. Default 10.
%     .reportMixtureCut  also fit a two-component Gaussian mixture to the log
%                        statistic and report the implied cut. Default false.
%
% Outputs
% -------
% selection : struct
%     .frames             indices of motion-free samples
%     .statistic          the per-sample motion statistic
%     .threshold          the applied cut
%     .fractionKept       proportion of samples retained
%     .mixtureThreshold   the alternative cut, or NaN if not requested
%     .mixtureFraction    proportion that alternative would retain, or NaN
%
% Dependency note
% ---------------
% Line-noise removal uses fft_clean, whose licensing is unresolved (see
% PIPELINE_SCRIPTS_AUDIT.md section H). When it is not on the path this function
% falls back to a local notch implementation and warns, so the tutorial still runs;
% results then differ slightly from the published run. Use the original when
% reproducing published numbers.
%
% Requires Statistics and Machine Learning Toolbox for the optional mixture fit.
% MATLAB R2019b compatible.

if nargin < 3, options = struct(); end
options = setDefault(options, 'multiplier', 2);
options = setDefault(options, 'medianWindow', 400);
options = setDefault(options, 'tailGuard', 10);
options = setDefault(options, 'reportMixtureCut', false);

if size(motionX, 1) > 1, motionX = mean(motionX, 1); end
if size(motionY, 1) > 1, motionY = mean(motionY, 1); end
motion = double([motionX(:), motionY(:)]);
nSample = size(motion, 1);

% Line noise sits in the flow estimate and would otherwise inflate the statistic.
if ~isempty(which('fft_clean'))
    motion = fft_clean(motion);   % double, matching the source call despite fft_clean's header
else
    warning('Motion:NoFftClean', ...
        ['fft_clean is not on the path; using a local notch fallback. Results ' ...
         'will differ slightly from the published run.']);
    motion = localNotch(motion);
end

residual = abs(double(motion) - movmedian(double(motion), options.medianWindow, 1));
statistic = residual(:,1) .* residual(:,2);

threshold = mean(statistic) * options.multiplier;
frames = find(statistic < threshold);
frames(frames > (nSample - options.tailGuard)) = [];

selection.frames        = frames(:)';
selection.statistic     = statistic(:)';
selection.threshold     = threshold;
selection.fractionKept  = numel(frames) / nSample;
selection.mixtureThreshold = NaN;
selection.mixtureFraction  = NaN;
selection.options       = options;

if options.reportMixtureCut
    [selection.mixtureThreshold, selection.mixtureFraction] = ...
        mixtureCut(statistic, nSample);
end
end

%% ========================================================= local functions

function [cut, fraction] = mixtureCut(statistic, nSample)
% Two-component Gaussian mixture on the log statistic, which is far closer to
% Gaussian than the raw product. The cut is the crossing point between the
% components, i.e. where a sample becomes more likely to have come from the
% moving component than the still one.
positive = statistic(statistic > 0);
if numel(positive) < 100
    cut = NaN; fraction = NaN; return
end
x = log(positive);
try
    model = fitgmdist(x, 2, 'Replicates', 5, 'Options', statset('MaxIter', 500));
catch fitError
    warning('Motion:MixtureFailed', 'Mixture fit did not converge: %s', fitError.message);
    cut = NaN; fraction = NaN; return
end
[~, stillComponent] = min(model.mu);
grid = linspace(min(x), max(x), 4000);
posteriorStill = posterior(model, grid');
crossing = find(diff(posteriorStill(:, stillComponent) > 0.5) ~= 0, 1, 'last');
if isempty(crossing)
    cut = NaN; fraction = NaN; return
end
cut = exp(grid(crossing));
fraction = sum(statistic < cut) / nSample;
end

function out = localNotch(signal)
% Minimal stand-in for fft_clean: zero frequency bins whose magnitude exceeds a
% robust multiple of the local spectral median, which removes narrow line-noise
% peaks while leaving the broadband motion spectrum alone.
% The input is real, so |X(k)| is already symmetric about the Nyquist bin. A
% threshold on magnitude therefore selects a Hermitian-symmetric set of bins on
% its own, and zeroing them leaves the inverse transform real. No explicit
% conjugate mirroring is needed, and attempting it by hand is where this kind of
% helper usually goes wrong.
out = signal;
for c = 1:size(signal, 2)
    offset = mean(signal(:,c));
    X = fft(signal(:,c) - offset);
    magnitude = abs(X);
    baseline = movmedian(magnitude, 201);
    spread = 1.4826 * median(abs(magnitude - baseline));
    isLine = magnitude > baseline + 6*spread;
    isLine(1) = false;                      % never remove DC
    X(isLine) = 0;
    reconstructed = ifft(X);
    assert(max(abs(imag(reconstructed))) < 1e-9 * max(abs(real(reconstructed))), ...
        'Motion:NotchNotReal', 'Notch filtering produced a complex signal.');
    out(:,c) = real(reconstructed) + offset;
end
end

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = value; end
end
