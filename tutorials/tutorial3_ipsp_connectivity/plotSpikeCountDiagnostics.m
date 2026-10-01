function plotSpikeCountDiagnostics(out, cellIndex, prm)
%PLOTSPIKECOUNTDIAGNOSTICS Show how spikes per pulse are read off the spectrum.
%
%   plotSpikeCountDiagnostics(out, cellIndex, prm)
%
%   out is the struct returned by ESTIMATESPIKESPERPULSE for one cell. Four
%   panels follow the estimator step by step:
%
%     1. the pulse-triggered trials, their median, and the fitted stereotyped
%        response -- what is about to be removed;
%     2. the residual after removing it, which is mostly spiking;
%     3. the per-trial autocovariance of that residual;
%     4. its power spectrum, with the centroid marked. Centroid times pulse
%        duration is the spike count.
%
%   The logic of panel 4 is the whole method: spikes are brief, so more of them
%   broadens the residual's spectrum. The estimator never identifies an
%   individual spike, which is the point -- during the pulse they cannot be
%   detected reliably.
%
%   Draws only; computes nothing that affects results.
%
%   See also ESTIMATESPIKESPERPULSE, TUTORIAL3A_BLUEPULSESPIKECOUNT.

tau = (-prm.nBack:prm.nFront) * prm.dt;
lag = (-(prm.tOff-prm.tOn+1):(prm.tOff-prm.tOn+1)) * prm.dt;

figure('Name', sprintf('Tutorial 3a -- cell %d', cellIndex), 'Color', 'w');

subplot(2,2,1)
% plot() of an nTrials-by-nTau matrix returns nTrials line handles, so the
% legend must be given explicit handles -- otherwise it labels the first three
% trial lines and every swatch comes out grey.
hTrials = plot(tau, out.stMat', 'Color', [.75 .75 .75]); hold on
hMedian = plot(tau, median(out.stMat,1), 'k', 'LineWidth', 1.5);
hTemplate = plot(tau, out.template, 'r', 'LineWidth', 1.5);
xlim([-.05 .1]); xlabel('time from trigger (s)'); ylabel('normalised V')
title(sprintf('cell %d: %d pulse-triggered trials', cellIndex, out.nTrials))
legend([hTrials(1) hMedian hTemplate], {'trials','median','fitted template'}, ...
    'Location','best'); legend boxoff

subplot(2,2,2)
plot(tau, out.stMatCor', 'Color', [.75 .75 .75]); hold on
plot(tau, mean(out.stMatCor,1), 'k', 'LineWidth', 1.5);
xlim([-.05 .1]); xlabel('time from trigger (s)'); ylabel('residual')
title('template removed: what is left is mostly spikes')

subplot(2,2,3)
imagesc(lag, 1:out.nTrials, out.autoMat);
xlabel('lag (s)'); ylabel('trial'); title('residual autocovariance')

subplot(2,2,4)
plot(out.freq, out.P, 'k', 'LineWidth', 1.2); hold on
centroid = sum(out.freq .* out.P) / sum(out.P);
plot([centroid centroid], ylim, 'r--', 'LineWidth', 1.5);
xlabel('frequency (Hz)'); ylabel('power (floor removed)')
title(sprintf('centroid %.1f Hz x %g s = %.3f spikes', ...
    centroid, prm.tPulse, out.nSpikesPulse))
legend({'spectrum','centroid'}, 'Location','best'); legend boxoff
end
