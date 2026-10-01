% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function out = detectWithPspBank(y, G, c, L, PhiBank, alphaMax)
%DETECTWITHPSPBANK Exact constrained RSS-improvement IPSP detector.
%   Fits 0 <= alpha <= ALPHAMAX under the null and
%   0 <= alpha <= ALPHAMAX, A >= 0 for every template under the
%   alternative. ALPHAMAX defaults to Inf, preserving the original
%   nonnegative-only fit. PhiBank columns are positive unit-peak PSPs.
%   L=[] performs the otherwise identical unwhitened calculation.

if nargin < 4 || isempty(L)
    L = [];
end
if nargin < 6 || isempty(alphaMax)
    alphaMax = Inf;
end
[y, G, c] = validateInputs(y, G, c);
validateattributes(alphaMax, {'numeric'}, {'scalar','real','nonnegative'}, mfilename, 'alphaMax');
validateattributes(PhiBank, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'PhiBank');
n = numel(y);
if size(PhiBank, 1) ~= n
    error('detectWithPspBank:DimensionMismatch', 'PhiBank must have one row per sample in y.');
end
if isempty(L)
    yW = y; GW = G; cW = c; PhiW = PhiBank;
else
    validateattributes(L, {'numeric'}, {'2d','real','finite','size',[n n]}, mfilename, 'L');
    yW = L \ y; GW = L \ G; cW = L \ c; PhiW = L \ PhiBank;
end

[Q, ~] = qr(GW, 0);
yR = yW - Q * (Q' * yW);
cR = cW - Q * (Q' * cW);
PhiR = PhiW - Q * (Q' * PhiW);
denom = cR' * cR;
if denom <= eps(max(1, denom))
    alpha0 = 0;
else
    alpha0 = min(alphaMax, max(0, (cR' * yR) / denom));
end
r0 = yR - alpha0 * cR;
rss0 = r0' * r0;

K = size(PhiR, 2);
Tk = zeros(1, K);
rss1 = zeros(1, K);
alpha1 = zeros(1, K);
amplitude = zeros(1, K);
options = optimoptions('lsqlin', 'Display', 'off');
for k = 1:K
    D = [cR, -PhiR(:, k)];
    % Cast only the small constrained problem because MATLAB R2019b's
    % bounded least-squares solvers require double inputs.
    z = lsqlin(double(D), double(yR), [], [], [], [], ...
        double([0; 0]), double([alphaMax; Inf]), [], options);
    r1 = yR - D * z;
    rss1(k) = r1' * r1;
    Tk(k) = max(0, rss0 - rss1(k));
    alpha1(k) = z(1);
    amplitude(k) = z(2);
end
[T, bestTemplateIndex] = max(Tk);
alphaHat = alpha1(bestTemplateIndex);
amplitudeHat = amplitude(bestTemplateIndex);
betaHat = GW \ (yW - alphaHat * cW + amplitudeHat * PhiW(:,bestTemplateIndex));

out.T = T;
out.Tk = Tk;
out.bestTemplateIndex = bestTemplateIndex;
out.rss0 = rss0;
out.rss1 = rss1;
out.alpha0 = alpha0;
out.alpha1 = alpha1;
out.amplitude = amplitude;
out.alphaHat = alphaHat;
out.amplitudeHat = amplitudeHat;
out.alphaMax = alphaMax;
out.alphaHitCap = isfinite(alphaMax) && abs(alphaHat - alphaMax) <= 1e-8 * max(1, alphaMax);
out.rss1Best = rss1(bestTemplateIndex);
out.yResidual = yR;
out.cResidual = cR;
out.PhiResidual = PhiR;
% Raw-scale components for diagnostic plots and downstream inspection.
out.betaHat = betaHat;
out.fittedBaseline = G * betaHat;
out.fittedCrosstalk = alphaHat * c;
out.fittedPsp = -amplitudeHat * PhiBank(:,bestTemplateIndex);
out.fittedTotal = out.fittedBaseline + out.fittedCrosstalk + out.fittedPsp;
out.residualRaw = y - out.fittedTotal;
end

function [y, G, c] = validateInputs(y, G, c)
validateattributes(y, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'y');
validateattributes(G, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'G');
validateattributes(c, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'c');
y = y(:); c = c(:);
if size(G, 1) ~= numel(y) || numel(c) ~= numel(y)
    error('detectWithPspBank:DimensionMismatch', 'y, G, and c must have the same number of rows.');
end
if rank(G) < size(G, 2)
    error('detectWithPspBank:RankDeficientBaseline', 'The baseline matrix G must have full column rank.');
end
end
