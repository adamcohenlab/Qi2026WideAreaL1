% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function out = selectCrosstalkDelayAndCapMatrix(Y, c, t, validMask, artifactMask, cfg)
%SELECTCROSSTALKDELAYANDCAPMATRIX Batched capped crosstalk-delay selection.
%   Y is nTime-by-nWaveforms.  Each output column/row corresponds to one
%   waveform, and is algebraically equivalent to calling
%   selectCrosstalkDelayAndCap for each column of Y.

validateattributes(Y, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'Y');
validateattributes(c, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'c');
validateattributes(t, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 't');
if nargin < 4 || isempty(validMask), validMask = true(size(Y,1), 1); end
if nargin < 5 || isempty(artifactMask), artifactMask = false(size(Y,1), 1); end
if nargin < 6 || isempty(cfg), cfg = struct; end
cfg = normaliseCfg(cfg);

c0 = c(:); t = t(:);
n = size(Y, 1); R = size(Y, 2);
if numel(c0) ~= n || numel(t) ~= n
    error('selectCrosstalkDelayAndCapMatrix:DimensionMismatch', ...
        'Y, c, and t must have the same number of rows.');
end
validMask = logical(validMask(:)); artifactMask = logical(artifactMask(:));
if numel(validMask) ~= n || numel(artifactMask) ~= n
    error('selectCrosstalkDelayAndCapMatrix:Mask', ...
        'validMask and artifactMask must match the rows of Y.');
end

c1 = [zeros(1, 1, 'like', c0); c0(1:end-1)];
CBank = [c0, c1];
usableMask = validMask & ~artifactMask;
cEnvelope = max(CBank, [], 2);
envelopePeak = max(cEnvelope);
if ~isfinite(envelopePeak) || envelopePeak <= eps(max(1, max(abs(CBank(:)))))
    error('selectCrosstalkDelayAndCapMatrix:InvalidTemplate', ...
        'Crosstalk template must have a positive peak.');
end
cMask = cEnvelope < cfg.crosstalkNegligibleFraction * envelopePeak & usableMask;
nSample = nnz(cMask);
if nSample <= cfg.crosstalkCapUpperTailCount
    error('selectCrosstalkDelayAndCapMatrix:TooFewLowCrosstalkSamples', ...
        'Need more than %d valid low-crosstalk samples; found %d.', ...
        cfg.crosstalkCapUpperTailCount, nSample);
end
percentileLevel = 100 * (nSample - cfg.crosstalkCapUpperTailCount) / nSample;
if ~isempty(cfg.crosstalkCapMinimumPercentile)
    percentileLevel = max(cfg.crosstalkCapMinimumPercentile, percentileLevel);
end
percentileLevel = min(percentileLevel, 100);
baselineUpper = prctile(double(Y(cMask,:)), percentileLevel, 1);

stimMask = inWindow(t, cfg.alphaCapStimWindow) & usableMask;
earlyMask = inWindow(t, cfg.crosstalkEarlyWindow) & usableMask;
if ~any(stimMask) || ~any(earlyMask)
    error('selectCrosstalkDelayAndCapMatrix:EmptyWindow', ...
        'The stimulation and early-delay windows must contain usable samples.');
end

templatePeak = max(CBank, [], 1);
if any(templatePeak <= eps(max(1, envelopePeak)))
    error('selectCrosstalkDelayAndCapMatrix:InvalidTemplate', ...
        'Each candidate crosstalk template must have a positive peak.');
end
stimulusPeak = max(double(Y(stimMask,:)), [], 1);
alphaMaxByDelay = max(0, (stimulusPeak - baselineUpper).' ./ templatePeak);

yEarly = double(Y(earlyMask,:)) - baselineUpper;
cEarly = double(CBank(earlyMask,:));
denom = sum(cEarly.^2, 1);
alphaUnconstrained = zeros(2, R);
validDenom = denom > eps(max(1, denom));
alphaUnconstrained(validDenom,:) = ...
    (cEarly(:,validDenom).' * yEarly) ./ denom(validDenom).';
alphaDelay = min(alphaMaxByDelay.', max(0, alphaUnconstrained));
residualNominal = yEarly - cEarly(:,1) * alphaDelay(1,:);
residualDelayed = yEarly - cEarly(:,2) * alphaDelay(2,:);
rssDelay = [sum(residualNominal.^2, 1); sum(residualDelayed.^2, 1)];
[~, bestDelayIndex] = min(rssDelay, [], 1);

linearIndex = sub2ind([R, 2], (1:R).', bestDelayIndex.');
out.cSelected = CBank(:, bestDelayIndex);
out.alphaMax = alphaMaxByDelay(linearIndex).';
out.baselineUpper = baselineUpper;
out.bestDelayIndex = bestDelayIndex;
out.delayRss = rssDelay.';
out.delayAlpha = alphaDelay.';
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
    error('selectCrosstalkDelayAndCapMatrix:WindowOrder', 'Window endpoints must be ordered [start end].');
end
end

function s = defaultField(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
