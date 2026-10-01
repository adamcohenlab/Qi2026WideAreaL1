function fig = plotSpikeletDiagnostics(mode, varargin)
%PLOTSPIKELETDIAGNOSTICS Figures for Tutorial 4. Computes nothing.
%
%   fig = PLOTSPIKELETDIAGNOSTICS(mode, ...) draws one figure. Every number
%   it draws is passed in; this function must not be where an analysis
%   decision lives.
%
%   MODES
%     'orthogonality'  (trig, pre, post, chunkIdx, nFrame)
%         The two cells' Hadamard patterns over one chunk, and the frames
%         that survive as triggers.
%     'staCheck'       (tau, ours, deposited, pre, post, dt)
%         Reconstructed versus deposited waveform for one pair.
%     'waveformRow'    (tau, staRow, distRow, dx, pre, dt, nBack)
%         Every post cell's raw waveform for one pre cell, coloured by
%         distance, with the amplitude windows marked.
%     'correction'     (out, distRow, dx, pre, dt)
%         The five-panel view of the two correction steps, the same layout
%         as the doPlot block of the source script.
%     'amplitudeVsDistance' (dist, ampRaw, ampCor, dx)
%         What the correction does to the distance profile.
%     'nullPair'       (nullAmp, ampCor, p, pre, post)
%         One pair's null distribution and where the observation sits.
%
%   See also TUTORIAL4A_ORTHOGONALTRIGGERS.

switch lower(mode)
    case 'orthogonality'
        fig = plotOrthogonality(varargin{:});
    case 'stacheck'
        fig = plotStaCheck(varargin{:});
    case 'waveformrow'
        fig = plotWaveformRow(varargin{:});
    case 'correction'
        fig = plotCorrection(varargin{:});
    case 'amplitudevsdistance'
        fig = plotAmplitudeVsDistance(varargin{:});
    case 'nullpair'
        fig = plotNullPair(varargin{:});
    otherwise
        error('Tutorial4:UnknownPlotMode', 'Unknown plot mode "%s".', mode);
end
end

%% ------------------------------------------------------------------------
function fig = plotOrthogonality(trig, pre, post, chunkIdx, nFrame)
f0 = (chunkIdx-1)*nFrame + 1;
idx = f0 : f0 + nFrame - 1;
t = (0:nFrame-1) * 1.27e-3;

trigFrames = trig.pairTriggers(pre, post);
inChunk = trigFrames(trigFrames >= idx(1) & trigFrames <= idx(end));
allHad = trig.tHad{pre};
allHad = allHad(allHad >= idx(1) & allHad <= idx(end));

fig = figure('Name', sprintf('orthogonality %d -> %d', pre, post), 'Color', 'w');
ax(1) = subplot(3,1,1);
plot(t, double(trig.blueHadOn(pre, idx)), 'b', 'LineWidth', 1); ylim([-0.2 1.2])
ylabel(sprintf('cell %d\nblue on', pre)); title(sprintf('Hadamard chunk %d', chunkIdx))
ax(2) = subplot(3,1,2);
plot(t, double(trig.blueHadOn(post, idx)), 'Color', [0 .5 0], 'LineWidth', 1); ylim([-0.2 1.2])
ylabel(sprintf('cell %d\nblue on', post))
ax(3) = subplot(3,1,3);
ortho = trig.blueHadOn(pre, idx) & ~trig.blueHadOn(post, idx);
area(t, double(ortho), 'FaceColor', [.85 .85 .85], 'EdgeColor', 'none'); hold on
plot((allHad - f0)*1.27e-3,  0.5*ones(size(allHad)),  '.', 'Color', [.6 .6 .6])
plot((inChunk - f0)*1.27e-3, 0.5*ones(size(inChunk)), 'r.', 'MarkerSize', 8)
ylim([-0.2 1.9]); ylabel('usable'); xlabel('Time in chunk (s)')
set(gca, 'YTick', [0 1])
legend({'orthogonal frames', 'pre spikes in Hadamard steps', 'surviving triggers'}, ...
    'Location', 'north', 'Box', 'off', 'Orientation', 'horizontal', 'FontSize', 8)
linkaxes(ax, 'x'); xlim([0 t(end)])
end

%% ------------------------------------------------------------------------
function fig = plotStaCheck(tau, ours, deposited, pre, post, dt)
fig = figure('Name', sprintf('STA check %d -> %d', pre, post), 'Color', 'w');
subplot(1,2,1)
plot(tau*dt*1e3, deposited, 'k', 'LineWidth', 2); hold on
plot(tau*dt*1e3, ours, 'r--', 'LineWidth', 1)
xlabel('Peri-spike time (ms)'); ylabel('Normalized spike height')
legend({'deposited', 'recomputed'}, 'Box', 'off')
title(sprintf('pair %d \\rightarrow %d', pre, post))
subplot(1,2,2)
plot(tau*dt*1e3, ours(:) - deposited(:), 'k')
xlabel('Peri-spike time (ms)'); ylabel('Difference')
title(sprintf('max |diff| = %.2e', max(abs(ours(:)-deposited(:)))))
end

%% ------------------------------------------------------------------------
function fig = plotWaveformRow(tau, staRow, distRow, dx, pre, dt, nBack)
d = distRow(:) * dx;
keep = (1:numel(d))' ~= pre;
[~, order] = sort(d);
order = order(ismember(order, find(keep)));

fig = figure('Name', sprintf('raw waveforms, pre %d', pre), 'Color', 'w');
subplot(1,3,1:2)
y = staRow(order,:) - mean(staRow(order,:), 2);
plot(tau*dt*1e3, y', 'LineWidth', 0.5); colororder(gca, jet(numel(order)))
hold on
plot(tau*dt*1e3, median(y,1), 'k', 'LineWidth', 2)
yl = ylim;
patch([tau(nBack+1) tau(nBack+2) tau(nBack+2) tau(nBack+1)]*dt*1e3, ...
    [yl(1) yl(1) yl(2) yl(2)], [1 .9 .9], 'FaceAlpha', .5, 'EdgeColor', 'none')
for f = [-4 4]
    patch(([f f+2 f+2 f]+0)*dt*1e3, [yl(1) yl(1) yl(2) yl(2)], ...
        [.9 .9 1], 'FaceAlpha', .5, 'EdgeColor', 'none')
end
set(gca,'Children',flipud(get(gca,'Children')))
xlabel('Peri-spike time (ms)'); ylabel('Normalized spike height')
title(sprintf('pre cell %d: every post cell, coloured near \\rightarrow far', pre))

subplot(1,3,3)
imagesc(1, sort(d(order)), reshape(jet(numel(order)), [], 1, 3))
set(gca, 'XTick', []); ylabel('Distance (\mum)'); axis on
end

%% ------------------------------------------------------------------------
function fig = plotCorrection(out, distRow, dx, pre, dt)
tau = (-out.nBack:out.nFront) * dt * 1e3;
d = distRow(:) * dx;
[~, order] = sort(d);
order = order(order ~= pre);
cmap = jet(numel(order));

% The artifacts being removed are ~3e-3 and the per-pair noise is ~8e-3, so
% all 319 individual waveforms drawn together hide exactly what these panels
% are about. Panels 1, 3 and 4 therefore show DISTANCE-BIN AVERAGES, using
% the same bins the trend is fitted in, which cuts the noise by sqrt(n) and
% leaves the trend, the common mode and the surviving spikelet all legible
% on one scale. Stage 4b is where the individual waveforms are shown.
binUse = reshape(unique(out.binOf(out.binOf > 0)), 1, []);
% -> nTau x nBin, so each column is one bin's average waveform.
binAvg = @(M) cell2mat(arrayfun(@(b) mean(M(out.binOf == b, :), 1)', ...
    binUse, 'UniformOutput', false));
binCmap = jet(numel(binUse));
binLabel = arrayfun(@(b) sprintf('%.0f \\mum', out.distCtr(b)*dx), binUse, ...
    'UniformOutput', false);

fig = figure('Name', sprintf('common-mode correction, pre %d', pre), ...
    'Color', 'w', 'Position', [80 200 1500 400]);

ax(1) = subplot(1,4,1);
plot(tau, binAvg(out.trRaw), 'LineWidth', 1); colororder(gca, binCmap)
title({'1. raw, mean subtracted', 'distance-bin averages'})
ylabel('Normalized spike height')
legend(binLabel, 'Box', 'off', 'FontSize', 7, 'Location', 'northwest')

ax2 = subplot(1,4,2);
plot(tau, out.yFit(:,binUse), 'LineWidth', 1); colororder(gca, binCmap)
ax2.YAxis.Exponent = 0;      % else the x10^-3 label lands on the title
title('2. the fitted distance trend')
xlabel('Peri-spike time (ms)')

ax(2) = subplot(1,4,3);
plot(tau, binAvg(out.trCor), 'LineWidth', 1); colororder(gca, binCmap); hold on
plot(tau, out.trBsl, 'k', 'LineWidth', 2)
title({'3. trend removed,', 'far-field common mode in black'})

ax(3) = subplot(1,4,4);
plot(tau, binAvg(out.trRes), 'LineWidth', 1); colororder(gca, binCmap); hold on
plot(tau, mean(out.trRes(out.isFar,:), 1), 'k', 'LineWidth', 2)
title({'4. corrected', 'nearest bin keeps its spikelet'})

linkaxes(ax, 'y')
sgtitle(sprintf('pre cell %d: panels 1, 3 and 4 share a scale', pre))
end

%% ------------------------------------------------------------------------
function fig = plotAmplitudeVsDistance(dist, ampRaw, ampCor, dx)
d = dist(:) * dx;
fig = figure('Name', 'amplitude vs distance', 'Color', 'w');
plot(d, ampRaw, '.', 'Color', [.7 .7 .7]); hold on
plot(d, ampCor, '.', 'Color', [.85 .2 .2])
yline(0, 'k:');
xlabel('Distance (\mum)'); ylabel('Spikelet amplitude (normalized)')
legend({'raw', 'corrected'}, 'Box', 'off')
end

%% ------------------------------------------------------------------------
function fig = plotNullPair(nullAmp, ampCor, p, pre, post)
fig = figure('Name', sprintf('null %d -> %d', pre, post), 'Color', 'w');
histogram(nullAmp, 60, 'Normalization', 'pdf', 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none')
hold on
yl = ylim;
plot([ampCor ampCor], yl, 'r', 'LineWidth', 2)
xlabel('Spikelet amplitude (normalized)'); ylabel('Null density')
title(sprintf('%d \\rightarrow %d:  p = %.4f  (%d nulls)', pre, post, p, numel(nullAmp)))
legend({'far-cell null', 'observed, corrected'}, 'Box', 'off')
end
