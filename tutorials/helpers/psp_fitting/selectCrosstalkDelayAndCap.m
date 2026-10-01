% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function out = selectCrosstalkDelayAndCap(y, c, t, validMask, artifactMask, cfg)
%SELECTCROSSTALKDELAYANDCAP Select a capped nominal/delayed crosstalk fit.
%   Builds nominal C0 and one-frame delayed C1 raw templates, derives a
%   conservative upper-tail baseline from samples where both are negligible,
%   calculates ALPHAMAX, and selects the delay with the smallest capped RSS
%   in a prespecified early window. Selection is intentionally performed
%   before masking for the final fit or applying a whitening transform.

validateattributes(y, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'y');
validateattributes(c, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'c');
validateattributes(t, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 't');
if nargin < 4 || isempty(validMask), validMask = true(numel(y), 1); end
if nargin < 5 || isempty(artifactMask), artifactMask = false(numel(y), 1); end
if nargin < 6 || isempty(cfg), cfg = struct; end
cfg = normaliseCfg(cfg);

y = y(:); c0 = c(:); t = t(:);
n = numel(y);
if numel(c0) ~= n || numel(t) ~= n
    error('selectCrosstalkDelayAndCap:DimensionMismatch', 'y, c, and t must have the same number of samples.');
end
validMask = logical(validMask(:)); artifactMask = logical(artifactMask(:));
if numel(validMask) ~= n || numel(artifactMask) ~= n
    error('selectCrosstalkDelayAndCap:Mask', 'validMask and artifactMask must match y.');
end

c1 = [zeros(1, 1, 'like', c0); c0(1:end-1)];
CBank = [c0, c1];
usableMask = validMask & ~artifactMask;
cEnvelope = max(CBank, [], 2);
envelopePeak = max(cEnvelope);
if ~isfinite(envelopePeak) || envelopePeak <= eps(max(1, max(abs(CBank(:)))))
    error('selectCrosstalkDelayAndCap:InvalidTemplate', 'Crosstalk template must have a positive peak.');
end
cMask = cEnvelope < cfg.crosstalkNegligibleFraction * envelopePeak & usableMask;
nSample = nnz(cMask);
if nSample <= cfg.crosstalkCapUpperTailCount
    error('selectCrosstalkDelayAndCap:TooFewLowCrosstalkSamples', ...
        'Need more than %d valid low-crosstalk samples; found %d.', cfg.crosstalkCapUpperTailCount, nSample);
end
percentileLevel = 100 * (nSample - cfg.crosstalkCapUpperTailCount) / nSample;
if ~isempty(cfg.crosstalkCapMinimumPercentile)
    percentileLevel = max(cfg.crosstalkCapMinimumPercentile, percentileLevel);
end
percentileLevel = min(percentileLevel, 100);
baselineUpper = prctile(double(y(cMask)), percentileLevel);

stimMask = inWindow(t, cfg.alphaCapStimWindow) & usableMask;
earlyMask = inWindow(t, cfg.crosstalkEarlyWindow) & usableMask;
if ~any(stimMask) || ~any(earlyMask)
    error('selectCrosstalkDelayAndCap:EmptyWindow', 'The stimulation and early-delay windows must contain usable samples.');
end

% The supplied crosstalk template is normally unit peak. Dividing by the
% selected raw-template peak preserves the cap's meaning if that convention
% was not met exactly.
templatePeak = max(CBank, [], 1);
if any(templatePeak <= eps(max(1, envelopePeak)))
    error('selectCrosstalkDelayAndCap:InvalidTemplate', 'Each candidate crosstalk template must have a positive peak.');
end
stimulusPeak = max(double(y(stimMask)));
alphaMaxByDelay = max(0, (stimulusPeak - baselineUpper) ./ double(templatePeak));

yEarly = double(y(earlyMask)) - baselineUpper;
rssDelay = Inf(1, 2);
alphaDelay = zeros(1, 2);
for s = 1:2
    cEarly = double(CBank(earlyMask, s));
    denom = cEarly' * cEarly;
    if denom <= eps(max(1, denom))
        alphaUnconstrained = 0;
    else
        alphaUnconstrained = (cEarly' * yEarly) / denom;
    end
    alphaDelay(s) = min(alphaMaxByDelay(s), max(0, alphaUnconstrained));
    residual = yEarly - alphaDelay(s) * cEarly;
    rssDelay(s) = residual' * residual;
end
[~, bestDelayIndex] = min(rssDelay);

out.cSelected = CBank(:, bestDelayIndex);
out.alphaMax = alphaMaxByDelay(bestDelayIndex);
out.baselineUpper = baselineUpper;
out.bestDelayIndex = bestDelayIndex;
out.delayRss = rssDelay;
out.delayAlpha = alphaDelay;
out.alphaMaxByDelay = alphaMaxByDelay;
out.cMask = cMask;
out.stimMask = stimMask;
out.earlyMask = earlyMask;
out.percentileLevel = percentileLevel;
out.templatePeak = templatePeak(bestDelayIndex);
end

function mask = inWindow(t, window)
mask = t >= window(1) & t <= window(2);
end

function cfg = normaliseCfg(cfg)
cfg = defaultField(cfg, 'alphaCapStimWindow', [0 0.020]);
cfg = defaultField(cfg, 'crosstalkEarlyWindow', cfg.alphaCapStimWindow);
cfg = defaultField(cfg, 'crosstalkNegligibleFraction', 0.05);
cfg = defaultField(cfg, 'crosstalkCapUpperTailCount', 20);
cfg = defaultField(cfg, 'crosstalkCapMinimumPercentile', []);
validateattributes(cfg.alphaCapStimWindow, {'numeric'}, {'vector','numel',2,'real','finite'}, mfilename, 'cfg.alphaCapStimWindow');
validateattributes(cfg.crosstalkEarlyWindow, {'numeric'}, {'vector','numel',2,'real','finite'}, mfilename, 'cfg.crosstalkEarlyWindow');
validateattributes(cfg.crosstalkNegligibleFraction, {'numeric'}, {'scalar','>',0,'<',1}, mfilename, 'cfg.crosstalkNegligibleFraction');
validateattributes(cfg.crosstalkCapUpperTailCount, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'cfg.crosstalkCapUpperTailCount');
if ~isempty(cfg.crosstalkCapMinimumPercentile)
    validateattributes(cfg.crosstalkCapMinimumPercentile, {'numeric'}, {'scalar','>=',0,'<=',100}, mfilename, 'cfg.crosstalkCapMinimumPercentile');
end
if cfg.alphaCapStimWindow(1) > cfg.alphaCapStimWindow(2) || cfg.crosstalkEarlyWindow(1) > cfg.crosstalkEarlyWindow(2)
    error('selectCrosstalkDelayAndCap:WindowOrder', 'Window endpoints must be ordered [start end].');
end
end

function s = defaultField(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
