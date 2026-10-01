function out = correctCommonMode(staRow, distRow, varargin)
%CORRECTCOMMONMODE Remove the distance trend and the far-field common mode.
%
%   out = CORRECTCOMMONMODE(staRow, distRow) takes one presynaptic cell's row
%   of spikelet waveforms — nCells-by-nTau, in normalized units, offset zero
%   at index nBack+1 — together with that cell's distances to every other
%   cell in pixels, and returns the corrected row. This is
%   `sparsePulseHad_spikelet_cm_correction.m` for a single presynaptic cell.
%
%   Options (defaults are the source literals):
%     'nBack'        20             offset-zero index is nBack+1
%     'dx'           6.5            um per pixel
%     'distEdges'    linspace(1,2000,11)/dx   distance bins, in pixels
%     'polyOrder'    3
%     'polyExclude'  -2:4           offsets left out of the polynomial fit
%     'cmExclude'    -2:2           offsets left out of the common-mode fit
%     'rBaseline'    800            um; cells beyond this define the common mode
%     'selfIndex'    []             presynaptic index, excluded from the fits
%
%   FIELDS
%     trRaw     nCells x nTau input, mean subtracted; the starting point
%     trCor     nCells x nTau after the distance-trend subtraction
%     trRes     nCells x nTau after both corrections; this is the output
%     trBsl     1 x nTau far-field common-mode waveform
%     c         nCells x 1 non-negative common-mode coefficients
%     yFit      nTau x nBin fitted distance-bin trends
%     binOf     1 x nCells bin assignment, 0 = outside every bin
%     varTot, varRes   per-post-cell variance before and after the second step
%
%   THE TWO THINGS BEING REMOVED, AND WHY THEY ARE DIFFERENT
%
%   Averaging a post cell's trace on a pre cell's spikes does not return a
%   flat line for unconnected pairs. Two artifacts survive, with different
%   shapes and different distance dependence.
%
%   1. A smooth, slow, distance-dependent trend. Optical crosstalk and
%      scattered light from the pre cell's own stimulation reach the post
%      cell, and how much depends on how far away it is. This is the term
%      that is worst exactly where real coupling is most likely, so it
%      cannot be left in. It is estimated by binning post cells by distance
%      from the pre cell, averaging within each bin, fitting a cubic in
%      time to each bin average, and interpolating that fit across distance
%      for every post cell. A cubic because the trend is slow relative to a
%      3-frame spikelet: it can bend, but it cannot manufacture a peak.
%
%   2. A common-mode waveform shared by every post cell regardless of
%      distance. The pre cell's spike is not an isolated event — it is
%      correlated with network activity, and the imaging shares optics and
%      illumination across the field. Cells beyond rBaseline are too far for
%      a gap junction, so whatever they show at offset zero is this shared
%      term. It is removed by regressing each post cell's waveform onto that
%      far-field average.
%
%   The coefficient is clamped at zero, `max(c, 0)`. A negative coefficient
%   would mean the post cell carried an inverted copy of the common mode,
%   which is not a physical possibility here, and letting it go negative
%   would ADD a positive bump at offset zero and manufacture spikelets.
%
%   WHY BOTH FITS EXCLUDE THE PEAK
%
%   Both corrections are estimated from data that includes the signal being
%   looked for, so both leave the peak region out: offsets -2:4 for the
%   polynomial, -2:2 for the common-mode regression. Without that, a genuine
%   spikelet in the nearest distance bin would be fitted as part of the
%   trend and subtracted away, and the test would be self-defeating. The
%   asymmetry of the polynomial exclusion, -2:4 rather than -2:2, matches the
%   slower falling side of the spikelet.
%
%   The exclusions are why the correction cannot be judged by how flat it
%   makes things. It is designed to leave the peak alone.
%
%   THE DIAGONAL
%   With the source's default distEdges the first edge is 1 um, so the pre
%   cell's own zero distance falls below every bin and is assigned bin 0,
%   which drops it from the trend fit; rBaseline drops it from the common
%   mode too. That is load-bearing. Pass 'selfIndex' to exclude it
%   explicitly as well, which matters if the row came from the deposit,
%   where the diagonal holds a self spike waveform about one spike height
%   tall rather than a cross average.
%
%   See also SPIKELETAMPLITUDE, TUTORIAL4C_COMMONMODECORRECTION.

%% ---------------------------------------------------------------- options
opt = struct('nBack', 20, 'dx', 6.5, 'distEdges', [], 'polyOrder', 3, ...
    'polyExclude', -2:4, 'cmExclude', -2:2, 'rBaseline', 800, 'selfIndex', []);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
if isempty(opt.distEdges)
    opt.distEdges = linspace(1, 2000, 11) / opt.dx;
end

[nCells, nTau] = size(staRow);
assert(numel(distRow) == nCells, 'Tutorial4:DistanceLength', ...
    'staRow has %d rows but distRow has %d entries.', nCells, numel(distRow));
distRow = reshape(double(distRow), 1, []);
nBack = opt.nBack;
nFront = nTau - nBack - 1;
xFit = -nBack:nFront;

%% ------------------------------------ step 1: distance-dependent trend
% Bin post cells by distance. Bin 0 means "outside every bin", which is
% where the pre cell itself lands.
[~, ~, binOf] = histcounts(distRow, opt.distEdges);
distCtr = mean([opt.distEdges(1:end-1); opt.distEdges(2:end)], 1);
nBin = numel(distCtr);

if ~isempty(opt.selfIndex)
    binOf(opt.selfIndex) = 0;
end

idxPoly = true(1, nTau);
idxPoly(nBack + 1 + opt.polyExclude) = false;

yFit = zeros(nTau, nBin);
for ib = 1:nBin
    inBin = (binOf == ib);
    if ~any(inBin)
        continue
    end
    y = mean(staRow(inBin,:), 1);
    y = y - mean(y);
    p = polyfit(xFit(idxPoly), y(idxPoly), opt.polyOrder);
    yFit(:,ib) = polyval(p, xFit);
end

tr = staRow - mean(staRow, 2);
binUse = unique(binOf(binOf > 0));
assert(numel(binUse) >= 2, 'Tutorial4:TooFewBins', ...
    'Only %d distance bins are populated; cannot interpolate the trend.', numel(binUse));
trend = interp1(distCtr(binUse), yFit(:,binUse)', distRow, 'linear', 'extrap');
trCor = tr - trend;

%% ------------------------------------ step 2: far-field common mode
isFar = distRow > opt.rBaseline / opt.dx;
if ~isempty(opt.selfIndex)
    isFar(opt.selfIndex) = false;
end
assert(any(isFar), 'Tutorial4:NoFarCells', ...
    'No cell is beyond %g um; the common mode cannot be estimated.', opt.rBaseline);
trBsl = mean(trCor(isFar,:), 1);

idxCm = true(1, nTau);
idxCm(nBack + 1 + opt.cmExclude) = false;

% SeeResiduals_vec(..., remOffset = 0, outChoice = 1) returns the regression
% coefficient only. Clamped at zero, see above.
c = max(SeeResiduals_vec(trCor(:,idxCm), trBsl(idxCm), 0, 1), 0);
trRes = trCor - c * trBsl;

%% ------------------------------------------------------------- assemble
out = struct();
out.trRaw   = tr;
out.trCor   = trCor;
out.trRes   = trRes;
out.trBsl   = trBsl;
out.c       = c;
out.yFit    = yFit;
out.binOf   = binOf;
out.distCtr = distCtr;
out.isFar   = isFar;
out.varTot  = var(trCor, [], 2);
out.varRes  = var(trRes, [], 2);
out.nBack   = nBack;
out.nFront  = nFront;
end
