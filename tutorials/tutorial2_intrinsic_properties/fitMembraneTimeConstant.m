function out = fitMembraneTimeConstant(blueOffSta, nBack, dt, varargin)
%FITMEMBRANETIMECONSTANT Decay constant of the voltage after blue light off.
%
%   out = FITMEMBRANETIMECONSTANT(blueOffSta, nBack, dt) fits
%
%       y = A1 * exp(-x/tau1) + C
%
%   to the first 41 samples after each cell's blue-off transition and returns
%
%     membraneC   nCells x 1, tau1 in milliseconds.
%     coefAll     nCells x 3, [A1 C tau1] as coeffvalues orders them.
%     rSqMem      nCells x 1, adjusted R squared.
%     atBound     nCells x 1 logical, tau1 resting on a fit bound.
%     skipped     nCells x 1 logical, rows the fit could not be run on.
%
%   FITMEMBRANETIMECONSTANT(..., 'Name', value) accepts:
%
%     'nSample'  samples after the transition to fit, default 41 (t = 0:40).
%     'TauUpper' upper bound on tau1 in frames, default 60.
%     'Verbose'  true (default).
%
%   COEFFICIENT ORDER
%   MATLAB's fittype orders coefficients alphabetically, so for
%   'A1*exp(-x/tau1) + C' coeffvalues returns [A1 C tau1] and tau1 is the
%   THIRD element, which is what the source's coefAll(:,3) picks up. The
%   source's doPlot title prints those three against the labels A1, B, C and
%   tau1, a leftover from the four-parameter variant just above it; the fit
%   itself is correct and the label is not.
%
%   THE BOUNDS BIND
%   The source constrains A1 >= 0, tau1 in [0 60] frames, which is 0 to
%   76.2 ms. Both ends are reached in the deposited data: two QC-passing cells
%   in M-YQ0201-27 sit at tau1 = 0 and one in M-YQ0201-29 sits at 76.2 ms. A
%   railed fit is not an estimate, so `atBound` flags them instead of letting
%   them enter a population histogram as if they were measurements.
%
%   WHAT IS BEING MEASURED
%   Turning the light off removes the CheRiff conductance as a step, and the
%   membrane relaxes to rest with its own time constant. The blue-off average
%   pools all five down-steps in the protocol, the four square-pulse offsets
%   and the ramp offset, exactly as the source does, so the estimate is an
%   average over a range of preceding depolarisations rather than a
%   small-signal constant.
%
%   Requires the Curve Fitting Toolbox.
%
%   See also SELFTRIGGEREDAVERAGE, TUTORIAL2D_ADPANDMEMBRANECONSTANT.

opt = struct('nSample', 41, 'TauUpper', 60, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end

nCells = size(blueOffSta, 1);
dbExp = fittype('A1*exp(-x/tau1) + C');
coefAll = zeros(nCells, 3);
rSqMem = zeros(nCells, 1);
skipped = false(nCells, 1);

x = (0:opt.nSample-1)';
for ii = 1:nCells
    y = blueOffSta(ii, nBack + 1 + (0:opt.nSample-1))';
    % The source guards with isnan(std(row)); a row left at zero by
    % get_sta_mat_self has zero variance and no decay to fit, so both cases
    % are skipped here and recorded.
    if any(~isfinite(y)) || std(y) == 0
        skipped(ii) = true;
        continue;
    end
    [f, gof] = fit(x, y, dbExp, 'startPoint', [1, -1, 30], ...
        'lower', [0 -inf 0], 'upper', [inf inf opt.TauUpper]);
    coefAll(ii, :) = coeffvalues(f);
    rSqMem(ii) = gof.adjrsquare;
end

tau1 = coefAll(:, 3);
out = struct();
out.membraneC = tau1 * dt / 1e-3;          % frames -> ms
out.coefAll = coefAll;
out.rSqMem = rSqMem;
out.skipped = skipped;
out.atBound = ~skipped & (tau1 <= 0 | tau1 >= opt.TauUpper - 1e-9);
out.tauUpperMs = opt.TauUpper * dt / 1e-3;

if opt.Verbose
    good = ~skipped & ~out.atBound;
    fprintf('Membrane time constant over %d cells\n', nCells);
    fprintf('  fitted %d, skipped %d, at a bound %d\n', ...
        sum(~skipped), sum(skipped), sum(out.atBound));
    fprintf('  tau (ms): median %.2f [%.2f %.2f] over unrailed fits\n', ...
        median(out.membraneC(good)), min(out.membraneC(good)), max(out.membraneC(good)));
    fprintf('  adjusted R^2: median %.3f\n', median(rSqMem(~skipped)));
end
end
