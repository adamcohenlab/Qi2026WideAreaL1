% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function out = detectWithPspBankMatrix(Y, G, C, L, PhiBank, alphaMax)
%DETECTWITHPSPBANKMATRIX Batched exact constrained IPSP detection.
%   Each column of Y and C is one waveform/crosstalk-template pair.  The
%   baseline is profiled once and every calibration-null waveform is fitted
%   simultaneously.  ALPHAMAX is one nonnegative upper bound per waveform
%   (or a scalar shared bound) and defaults to Inf. The implementation
%   evaluates the relevant active sets of the two-coefficient bounded
%   least-squares problem, which is algebraically equivalent to calling
%   LSQLIN once per waveform/template but avoids that inner loop.
%
%   L=[] performs the unwhitened calculation.  The outputs are R-by-1,
%   where R is the number of waveform columns.

if nargin < 4 || isempty(L)
    L = [];
end
validateattributes(Y, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'Y');
validateattributes(G, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'G');
validateattributes(C, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'C');
validateattributes(PhiBank, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'PhiBank');
n = size(Y,1); R = size(Y,2);
if size(G,1) ~= n || size(PhiBank,1) ~= n
    error('detectWithPspBankMatrix:DimensionMismatch', 'Y, G, and PhiBank must have the same number of rows.');
end
if isvector(C)
    C = repmat(C(:), 1, R);
end
if ~isequal(size(C), size(Y))
    error('detectWithPspBankMatrix:DimensionMismatch', 'C must be nTime-by-nWaveforms, or one nTime-vector.');
end
if rank(G) < size(G,2)
    error('detectWithPspBankMatrix:RankDeficientBaseline', 'The baseline matrix G must have full column rank.');
end
if nargin < 6 || isempty(alphaMax)
    alphaMax = Inf(1, R);
elseif isscalar(alphaMax)
    alphaMax = repmat(alphaMax, 1, R);
else
    alphaMax = alphaMax(:).';
end
if numel(alphaMax) ~= R || any(isnan(alphaMax) | ~isreal(alphaMax) | alphaMax < 0)
    error('detectWithPspBankMatrix:AlphaMax', 'alphaMax must be a nonnegative scalar or one value per waveform.');
end

% Do the small per-pair regression in double precision. This also keeps
% behavior stable when the recorded tensor is stored as single.
Y = double(Y); G = double(G); C = double(C); PhiBank = double(PhiBank);
if isempty(L)
    YW = Y; GW = G; CW = C; PhiW = PhiBank;
else
    validateattributes(L, {'numeric'}, {'2d','real','finite','size',[n n]}, mfilename, 'L');
    L = double(L);
    YW = L \ Y; GW = L \ G; CW = L \ C; PhiW = L \ PhiBank;
end

[Q, ~] = qr(GW, 0);
YR = YW - Q * (Q' * YW);
CR = CW - Q * (Q' * CW);
PhiR = PhiW - Q * (Q' * PhiW);

cc = sum(CR.^2, 1);
cy = sum(CR .* YR, 1);
alpha0 = zeros(1, R);
validC = cc > eps(max(1, max(cc)));
alpha0(validC) = min(alphaMax(validC), max(0, cy(validC) ./ cc(validC)));
R0 = YR - CR .* alpha0;
rss0 = sum(R0.^2, 1);

K = size(PhiR, 2);
bestT = zeros(1, R);
bestIndex = ones(1, R);
bestRss1 = rss0;
bestAlpha = alpha0;
bestAmplitude = zeros(1, R);

for k = 1:K
    w = -PhiR(:,k);
    ww = w' * w;
    alphaK = alpha0;
    amplitudeK = zeros(1, R);
    rssK = rss0;

    if ww > eps(max(1, ww))
        % Face with crosstalk inactive and PSP active.
        wy = w' * YR;
        amplitudeOnly = max(0, wy ./ ww);
        rssPspOnly = sum((YR - w .* amplitudeOnly).^2, 1);
        usePspOnly = rssPspOnly < rssK;
        alphaK(usePspOnly) = 0;
        amplitudeK(usePspOnly) = amplitudeOnly(usePspOnly);
        rssK(usePspOnly) = rssPspOnly(usePspOnly);

        % Interior face with both nonnegative coefficients active.
        cw = sum(CR .* w, 1);
        determinant = cc .* ww - cw.^2;
        determinantTolerance = 1e-12 * max(1, cc .* ww);
        validBoth = determinant > determinantTolerance;
        alphaBoth = zeros(1, R); amplitudeBoth = zeros(1, R);
        alphaBoth(validBoth) = (ww .* cy(validBoth) - cw(validBoth) .* wy(validBoth)) ./ determinant(validBoth);
        amplitudeBoth(validBoth) = (cc(validBoth) .* wy(validBoth) - cw(validBoth) .* cy(validBoth)) ./ determinant(validBoth);
        validBoth = validBoth & alphaBoth >= 0 & alphaBoth <= alphaMax & amplitudeBoth >= 0;
        if any(validBoth)
            rssBoth = sum((YR - CR .* alphaBoth - w .* amplitudeBoth).^2, 1);
            useBoth = validBoth & rssBoth < rssK;
            alphaK(useBoth) = alphaBoth(useBoth);
            amplitudeK(useBoth) = amplitudeBoth(useBoth);
            rssK(useBoth) = rssBoth(useBoth);
        end

        % Face where the crosstalk coefficient reaches its finite upper
        % bound. The PSP amplitude is then solved exactly under A >= 0.
        atCap = isfinite(alphaMax);
        if any(atCap)
            alphaCap = zeros(1, R);
            alphaCap(atCap) = alphaMax(atCap);
            residualAtCap = YR - CR .* alphaCap;
            amplitudeAtCap = max(0, (w' * residualAtCap) ./ ww);
            rssAtCap = sum((residualAtCap - w .* amplitudeAtCap).^2, 1);
            useCap = atCap & rssAtCap < rssK;
            alphaK(useCap) = alphaCap(useCap);
            amplitudeK(useCap) = amplitudeAtCap(useCap);
            rssK(useCap) = rssAtCap(useCap);
        end
    end

    Tk = max(0, rss0 - rssK);
    if k == 1
        bestT = Tk; bestIndex(:) = 1; bestRss1 = rssK;
        bestAlpha = alphaK; bestAmplitude = amplitudeK;
    else
        useTemplate = Tk > bestT;
        bestT(useTemplate) = Tk(useTemplate);
        bestIndex(useTemplate) = k;
        bestRss1(useTemplate) = rssK(useTemplate);
        bestAlpha(useTemplate) = alphaK(useTemplate);
        bestAmplitude(useTemplate) = amplitudeK(useTemplate);
    end
end

out.T = bestT(:);
out.bestTemplateIndex = bestIndex(:);
out.rss0 = rss0(:);
out.rss1Best = bestRss1(:);
out.alpha0 = alpha0(:);
out.alphaHat = bestAlpha(:);
out.amplitudeHat = bestAmplitude(:);
out.alphaMax = alphaMax(:);
out.alphaHitCap = isfinite(alphaMax(:)) & abs(bestAlpha(:) - alphaMax(:)) <= 1e-8 .* max(1, alphaMax(:));
end
