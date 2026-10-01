function plotIpspDetectionDiagnostics(results, rowIndex, pspTensor, crsTalkTemp, ...
    t, mask, Gmasked, bank, nullStatistics, figureName, thetaFits)
%PLOTIPSPDETECTIONDIAGNOSTICS Show why one pair was, or was not, called a synapse.
%
%   plotIpspDetectionDiagnostics(results, rowIndex, pspTensor, crsTalkTemp, ...
%       t, mask, Gmasked, bank, nullStatistics, figureName)
%
%   Draws three panels for the pair in results(rowIndex,:):
%
%     1. the observed waveform with the fitted model pulled apart into its
%        three parts -- low-order baseline, positive crosstalk, negative IPSP;
%     2. the fitted IPSP on its own, with its kinetics;
%     3. the empirical null distribution of the detection statistic, with the
%        observed value marked.
%
%   Panel 3 is the one that carries the argument. The statistic T is how much
%   the best PSP template improves the fit, and it has no tabulated null
%   distribution -- so the null is built from far-away stimulation, which
%   produces the same crosstalk but cannot produce a synapse. A detection means
%   the observed improvement sits outside that distribution. Seeing the observed
%   value against the histogram is the whole method in one picture.
%
%   This function draws only. It recomputes the baseline coefficients by least
%   squares from the stored crosstalk and amplitude estimates rather than
%   requiring them to be carried in the table, so nothing it does can change a
%   result.
%
%   nullStatistics may be [] to omit panel 3.
%
%   See also TUTORIAL3D_DETECTANDFITIPSP.

if nargin < 10 || isempty(figureName), figureName = 'IPSP detection'; end

r    = results(rowIndex,:);
pre  = r.preCell;  post = r.postCell;
y    = squeeze(double(pspTensor(pre, post, :)));

% crosstalk template for the RECORDED cell, with the delay the detector chose
c = crsTalkTemp(post,:).';
if r.crosstalkDelayIndex == 2
    c = [0; c(1:end-1)];
end

% the PSP shape the bank selected
phi = zeros(size(y));
if isfinite(r.bestTemplateIndex)
    phi(mask) = bank.Phi(:, r.bestTemplateIndex);
end

% re-profile the baseline at the stored crosstalk scale and IPSP amplitude:
%   y = G*beta + alpha*c - A*phi
alpha = r.crosstalkScaleHat;
A     = r.amplitudeHat;
resid = y(mask) - alpha*c(mask) + A*phi(mask);
beta  = Gmasked \ resid;

baselineFit = nan(size(y)); baselineFit(mask) = Gmasked*beta;
crossFit    = alpha * c;
ipspFit     = -A * phi;
totalFit    = baselineFit + crossFit + ipspFit;

tms = t*1e3;
nPanel = 2 + ~isempty(nullStatistics);
figure('Name', sprintf('%s: %d -> %d', figureName, pre, post), 'Color', 'w');

% ---------------------------------------------------------------- panel 1
subplot(1, nPanel, 1)
plot(tms, y, 'k', 'LineWidth', 1.4); hold on
plot(tms, crossFit, 'Color', [0 .45 .85], 'LineWidth', 1.1);
plot(tms, baselineFit, 'Color', [.6 .6 .6], 'LineWidth', 1.1);
plot(tms, ipspFit, 'r', 'LineWidth', 1.4);
plot(tms, totalFit, '--', 'Color', [0 .6 .2], 'LineWidth', 1.2);
plot([0 0], ylim, 'b:', 'HandleVisibility', 'off');
xlim([-20 150]); xlabel('time from pulse (ms)'); ylabel('normalised \DeltaV')
title(sprintf('%d \\rightarrow %d   p = %.4g%s', pre, post, r.pValue, ...
    ternary(r.isDetected, '  (detected)', '')))
legend({'observed','crosstalk \alpha c','baseline G\beta','IPSP -A\phi','total fit'}, ...
    'Location','southeast'); legend boxoff

% ---------------------------------------------------------------- panel 2
% IMPORTANT: the red curve is the BANK TEMPLATE that won the detection search,
% not the reported kinetics. The bank is a coarse grid whose only job is to
% detect; its parameters frequently rail against the grid edges (here
% latency <= 10 ms, tauRise <= 8 ms, tauDecay <= 100 ms). The kinetics that get
% reported come from the continuous fit in section 9 of the driver, which is
% unconstrained by the grid. Where that fit is available it is overlaid, and
% the two are usually quite different -- which is the point of doing it.
subplot(1, nPanel, 2)
plot(tms, ipspFit, 'r', 'LineWidth', 1.6); hold on

refined = [];
if nargin >= 11 && ~isempty(thetaFits) && isfield(thetaFits,'table') && ~isempty(thetaFits.table)
    q = find(thetaFits.table.preCell == pre & thetaFits.table.postCell == post, 1);
    if ~isempty(q), refined = thetaFits.table(q,:); end
end
if ~isempty(refined)
    u = t - refined.independentLatency;
    tr = refined.independentTauRise; td = refined.independentTauDecay;
    g = (exp(-u/td) - exp(-u/tr)) .* (u >= 0);
    upk = tr*td/(td-tr)*log(td/tr);
    g = g / ((exp(-upk/td) - exp(-upk/tr)));
    plot(tms, -refined.independentAmplitude * g, 'Color', [0 .5 0], 'LineWidth', 1.6);
    legend({'bank template (detection)','continuous fit (reported)'}, ...
        'Location','southeast'); legend boxoff
end
% guide lines: kept out of the legend
plot([0 0], ylim, 'b:', 'HandleVisibility', 'off');
plot(xlim, [0 0], 'k:', 'HandleVisibility', 'off');
xlim([-20 150]); xlabel('time from pulse (ms)'); ylabel('normalised \DeltaV')

if isfinite(r.bestTemplateIndex)
    ttl = sprintf(['bank template: A %.4f, lat %.1f, \\tau_r %.1f, \\tau_d %.1f ms\n' ...
                   '(grid-limited -- detection only)'], ...
        A, r.latencyHat*1e3, r.tauRiseHat*1e3, r.tauDecayHat*1e3);
    if ~isempty(refined)
        ttl = sprintf('%s\ncontinuous fit: A %.4f, lat %.1f, \\tau_r %.1f, \\tau_d %.1f ms', ...
            ttl, refined.independentAmplitude, refined.independentLatency*1e3, ...
            refined.independentTauRise*1e3, refined.independentTauDecay*1e3);
    end
    title(ttl)
else
    title('no template selected')
end

% ---------------------------------------------------------------- panel 3
if ~isempty(nullStatistics)
    subplot(1, nPanel, 3)
    Tn = nullStatistics(:);
    histogram(Tn, 60, 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none'); hold on
    % set the log scale BEFORE reading ylim: on a linear axis the lower limit
    % is 0, and a line drawn to 0 does not render once the axis becomes log.
    set(gca, 'YScale', 'log');
    yl = ylim; yl(1) = max(yl(1), 0.7);
    plot([r.Tobs r.Tobs], yl, 'r', 'LineWidth', 2);
    ylim(yl); xlim([0 max([Tn(:); r.Tobs])*1.05]);
    xlabel('T = RSS_0 - RSS_1'); ylabel('null waveforms')
    nGE = sum(Tn >= r.Tobs);
    title(sprintf(['empirical null, n = %d\n%d null(s) reached T_{obs}\n' ...
                   'p = (1+%d)/(1+%d) = %.4g'], numel(Tn), nGE, nGE, numel(Tn), r.pValue))
    legend({'null statistics','observed T'}, 'Location','northeast'); legend boxoff
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
