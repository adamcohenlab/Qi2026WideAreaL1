function fit = fitSpikeHeightMixture(detection, options)
% FITSPIKEHEIGHTMIXTURE  Recover the spike-height distribution, and both error rates.
%
% fit = FITSPIKEHEIGHTMIXTURE(detection, options)
%
% The detected spike heights are a contaminated, truncated sample: contaminated by
% the false positives the threshold admits by construction, and truncated because
% everything below threshold was never detected. This models them as a two-component
% mixture of truncated Gaussians and recovers the underlying distribution, which
% yields three things the raw detection cannot give:
%
%   1. a per-spike posterior probability of being a false positive
%   2. the false-negative rate, as the mass of the true-spike component below
%      threshold
%   3. an unbiased mean spike height, corrected for both truncation and
%      contamination, which is the per-cell normalization used downstream
%
% The model
% ---------
%   pdf(v)  =  p * p1(v)  +  (1 - p) * p2(v)
%
%   p1  false positives. Not free: it is fitted directly to the retained negative
%       peak heights, which are a measured sample of exactly this distribution.
%   p2  true spike heights. Mean and variance unknown; this is what we want.
%   p   mixing weight. Also not free: it is the expected number of false positives
%       in the analysis window divided by the number of spikes detected in it.
%
% Only p2's two parameters are actually free. The fit proceeds outward from what is
% measured: first the noise scale from the negative excursions, then p1 seeded from
% and bounded near it, then p2 with p1 and p held fixed.
%
% Inputs
% ------
% detection : struct from detectSpikesMatchedFilter
%
% options : struct
%   Required:
%     .framePeriod         seconds per sample
%   Optional, defaulted to the published values:
%     .analysisWindow      seconds from the start of the recording over which
%                          spike heights are fitted. Default 20*60. This limits
%                          photobleaching bias on spike height; it is a heuristic,
%                          not a measured bleach constant.
%     .falsePositiveRate   expected false positives per second per cell. Defaults
%                          to detection.options.falsePositiveRate so the two stages
%                          cannot silently disagree.
%     .maxMixingWeight     cap on p. Default 0.99.
%     .minFitPercentile    the retry loop trims the upper tail in 0.1 percentile
%                          steps until mle converges; below this it gives up.
%                          Default 90.
%     .boundaryGuard       samples around each chunk boundary to exclude.
%                          Default [-100 100].
%     .chunkLength         samples per chunk. Defaults to detection.options.chunkLength.
%
% Outputs
% -------
% fit : struct
%     .mixingWeight            p
%     .falsePositiveParams     [mu sigma] of p1
%     .trueSpikeParams         [mu sigma] of p2
%     .noiseParams             [mu sigma] of the zero-truncated noise fit
%     .threshold               the detection threshold the components are truncated at
%     .meanSpikeHeight         trueSpikeParams(1), the unbiased mean height
%     .falseNegativeRate       normcdf(threshold, mu2, sigma2): the fraction of true
%                              spikes falling below threshold. The source scripts
%                              never evaluate this although the fit supplies it.
%     .falsePositiveProbability per-spike pFP, aligned with .analysedSpikes
%     .analysedSpikes          spike indices inside the analysis window
%     .analysedValues          their filtered heights
%     .fitPercentile           [p1 p2] percentile actually accepted by the retry
%                              loops. The source does not record this, so the fit
%                              window used is not recoverable from saved output.
%     .expectedFalsePositivesInWindow
%
% Densities are returned as PARAMETERS, not function handles. results.mat stores
% these components as anonymous handles whose parent workspace was a Live Editor
% temp file; those handles no longer resolve on load, so the saved densities are
% unusable while the parameters beside them are fine. Use truncatedNormalPdf below.
%
% Requires Statistics and Machine Learning Toolbox (mle, normpdf, normcdf, prctile).
% MATLAB R2019b compatible.

%% ---------------------------------------------------------------- options
if nargin < 2, options = struct(); end
options = setDefault(options, 'analysisWindow', 20*60);
options = setDefault(options, 'maxMixingWeight', 0.99);
options = setDefault(options, 'minFitPercentile', 90);
options = setDefault(options, 'boundaryGuard', [-100 100]);
options = setDefault(options, 'falsePositiveRate', detection.options.falsePositiveRate);
options = setDefault(options, 'chunkLength', detection.options.chunkLength);
assert(isfield(options, 'framePeriod'), 'Mixture:MissingOption', ...
    'options.framePeriod is required.');

dt        = options.framePeriod;
filtered  = detection.filtered;
threshold = detection.threshold;
nSample   = numel(filtered);

%% ------------------------------------------------- restrict to the window
% Spikes late in the recording sit on a more bleached trace, so the height
% distribution is fitted only over the leading window, with chunk boundaries
% excluded.
lastSample = min(round(options.analysisWindow / dt), nSample);
boundary = boundaryIndices(nSample, options.chunkLength, options.boundaryGuard);
isBoundary = false(1, nSample);
isBoundary(boundary) = true;

spikes = detection.spikes(:)';
spikes = spikes(spikes <= lastSample & ~isBoundary(spikes));
values = filtered(spikes);
assert(~isempty(values), 'Mixture:NoSpikesInWindow', ...
    'No spikes fall inside the analysis window; cannot fit.');

% The mixing weight is a count ratio, so the expected false positives must be
% counted over the same window the detected spikes were counted over.
expectedFP = options.falsePositiveRate * lastSample * dt;
p = min(expectedFP / numel(values), options.maxMixingWeight);

%% ------------------------------------- step 1: noise scale, truncated at 0
% Negative excursions of the filtered trace are pure noise. Fitting a Gaussian
% truncated at zero to their magnitudes fixes the noise scale, which then seeds
% and bounds the false-positive component.
negativeMagnitude = -double(filtered(filtered < 0));
noiseParams = mle(negativeMagnitude(negativeMagnitude < prctile(negativeMagnitude, 99)), ...
    'pdf', @(x, mu, sigma) truncatedNormalPdf(x, mu, sigma, 0), ...
    'Start', double([0 std(negativeMagnitude)]), 'LowerBound', [-inf -inf]);

%% ------------------- step 2: false-positive component, truncated at threshold
% Fitted to the retained negative peaks, which are a measured sample of the
% false-positive height distribution rather than an assumption about it.
[fpParams, fpPercentile] = fitWithTailTrim( ...
    double(detection.negativePeakValues(:)), ...
    @(x, mu, sigma) truncatedNormalPdf(x, mu, sigma, threshold), ...
    double(noiseParams), noiseParams - abs(noiseParams)/5, [inf inf], ...
    options.minFitPercentile);

%% ----------------- step 3: true-spike component, with p1 and p held fixed
mixturePdf = @(x, mu, sigma) ...
    p .* truncatedNormalPdf(x, fpParams(1), fpParams(2), threshold) + ...
    (1 - p) .* truncatedNormalPdf(x, mu, sigma, threshold);
[spikeParams, spikePercentile] = fitWithTailTrim(double(values(:)), mixturePdf, ...
    double([mean(values) std(values)]), [threshold 0], [max(values) inf], ...
    options.minFitPercentile);

%% -------------------------------------------------- posteriors and rates
pFP = p .* truncatedNormalPdf(values, fpParams(1), fpParams(2), threshold);
pTP = (1 - p) .* truncatedNormalPdf(values, spikeParams(1), spikeParams(2), threshold);
falsePositiveProbability = pFP ./ (pFP + pTP);

% The false-negative rate is the mass of the fitted true-spike component that lies
% below the detection threshold: real spikes too small to have been detected. The
% truncated fit makes this available for free, and the source scripts never take it.
falseNegativeRate = normcdf(threshold, spikeParams(1), spikeParams(2));

% The true-spike mean is bounded below by the threshold, because a component whose
% mean sits under the threshold is not identifiable from data that was truncated
% there. When the optimiser pins the mean against that bound the fit has failed to
% separate the components, and the false-negative rate degenerates to exactly 50%
% by construction -- a number that looks like a measurement and is not one. This
% happens for cells whose spikes are barely above the noise, which are exactly the
% cells where the error rates matter most, so it is flagged rather than returned
% silently.
meanAtLowerBound = abs(spikeParams(1) - threshold) < 1e-6 * max(1, abs(threshold));
if meanAtLowerBound
    warning('Mixture:MeanAtBound', ...
        ['The true-spike mean pinned at the detection threshold (%.4g), so the ' ...
         'two components are not separable for this cell. falseNegativeRate is ' ...
         '%.3f only because of that bound and must not be interpreted. Treat the ' ...
         'spike-height fit for this cell as failed.'], threshold, falseNegativeRate);
end

%% ------------------------------------------------------------------ pack
fit.mixingWeight             = p;
fit.falsePositiveParams      = fpParams;
fit.trueSpikeParams          = spikeParams;
fit.noiseParams              = noiseParams;
fit.threshold                = threshold;
fit.meanSpikeHeight          = spikeParams(1);
fit.falseNegativeRate        = falseNegativeRate;
fit.meanAtLowerBound         = meanAtLowerBound;
fit.falsePositiveProbability = falsePositiveProbability;
fit.analysedSpikes           = spikes;
fit.analysedValues           = values;
fit.fitPercentile            = [fpPercentile spikePercentile];
fit.expectedFalsePositivesInWindow = expectedFP;
fit.options                  = options;
end

%% ========================================================= local functions

function y = truncatedNormalPdf(x, mu, sigma, lowerLimit)
% Normal density conditioned on x > lowerLimit.
y = normpdf(x, mu, sigma) ./ (1 - normcdf(lowerLimit, mu, sigma));
end

function [params, percentileUsed] = fitWithTailTrim(x, pdfFun, start, lower, upper, minPercentile)
% mle can fail to converge when a few extreme values dominate the likelihood. The
% source pipeline retries with a progressively trimmed upper tail; that is
% preserved here, and the accepted percentile is returned so the fit window is
% recoverable from the saved output, which it is not in the source.
% An all-infinite bound is omitted rather than passed. Supplying it would be
% mathematically equivalent but can send mle down its constrained-optimisation
% path, which need not land on the same numbers as the published unconstrained fit.
extraArgs = {};
if any(isfinite(lower)), extraArgs = [extraArgs, {'LowerBound', lower}]; end
if any(isfinite(upper)), extraArgs = [extraArgs, {'UpperBound', upper}]; end

percentileUsed = 100;
while true
    try
        params = mle(x(x <= prctile(x, percentileUsed)), 'pdf', pdfFun, ...
            'Start', start, extraArgs{:});
        return
    catch fitError
        percentileUsed = percentileUsed - 0.1;
        if percentileUsed < minPercentile
            error('Mixture:FitFailed', ...
                ['Maximum likelihood fit did not converge even after trimming to ' ...
                 'the %.1fth percentile. Last error: %s'], ...
                minPercentile, fitError.message);
        end
    end
end
end

function idx = boundaryIndices(nTotal, chunkLength, guard)
offsets = guard(1):guard(2);
idx = (0:chunkLength:nTotal)' + offsets;
idx = unique(idx(:));
idx = idx(idx >= 1 & idx <= nTotal)';
end

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = value; end
end
