function plotCrosstalkTemplateDiagnostics(selfWaveform, crsTalkTemp, fitInfo, ...
    tau, tBack, showCells)
%PLOTCROSSTALKTEMPLATEDIAGNOSTICS Show the crosstalk template against the data.
%
%   plotCrosstalkTemplateDiagnostics(selfWaveform, crsTalkTemp, fitInfo, ...
%       tau, tBack, showCells)
%
%   Top row: for each cell in showCells, that cell's response to its own blue
%   pulse (the tensor diagonal, essentially pure crosstalk) with the fitted
%   template over it, both scaled to unit early peak.
%
%   Bottom: every fitted template together. They are one shape family with a
%   single degree of freedom, the decay rate. That is what tutorial3d scales by
%   a non-negative alpha per pair.
%
%   Two details are worth seeing rather than reading about. Every template is
%   flat until frame +2, which is not arbitrary: under the convention-A trigger
%   tau = 0 is the last frame BEFORE the pulse, so +1 is the first commanded
%   frame and +2 the first the camera actually registers. And the rising phase
%   is built from the fitted DECAY constant, because the cell is spiking hard
%   during the pulse and its rise cannot be measured reliably.
%
%   Draws only.
%
%   See also BUILDCROSSTALKTEMPLATE, TUTORIAL3C_CROSSTALKTEMPLATE.

if nargin < 6 || isempty(showCells), showCells = 1:min(4, size(selfWaveform,1)); end
nShow = numel(showCells);
early = tBack + 1 + (0:20);

figure('Name', 'Tutorial 3c -- crosstalk templates', 'Color', 'w');

for k = 1:nShow
    ii = showCells(k);
    subplot(2, nShow, k)
    y = selfWaveform(ii,:);
    plot(tau, y / max(y(early)), 'k'); hold on
    plot(tau, crsTalkTemp(ii,:), 'r', 'LineWidth', 1.3);
    plot([0 0], ylim, 'b:');
    xlim([-.02 .1]); xlabel('time from trigger (s)')
    if isfinite(fitInfo.tauDecay(ii))
        title(sprintf('cell %d  (\\tau_{decay} rate %.3f/frame)', ii, fitInfo.tauDecay(ii)))
    else
        title(sprintf('cell %d', ii))
    end
    if k == 1
        ylabel('unit early peak')
        legend({'self waveform','template'}, 'Location','best'); legend boxoff
    end
end

subplot(2,1,2)
ok = fitInfo.ok;
plot(tau, crsTalkTemp(ok,:)', 'Color', [.3 .3 .3 .15]);
xlim([-.02 .15]); xlabel('time from trigger (s)'); ylabel('template')
title(sprintf(['all %d templates: one shape family, differing only in decay rate ' ...
               '(median %.1f ms)'], sum(ok), ...
               median(1./fitInfo.tauDecay(ok)) * (tau(2)-tau(1)) * 1e3))
end
