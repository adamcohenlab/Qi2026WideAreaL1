function figureHandle = plotSpikeFindingDiagnostics(result, session)
% PLOTSPIKEFINDINGDIAGNOSTICS  Diagnostic figure for one cell of tutorial 2.
%
% figureHandle = PLOTSPIKEFINDINGDIAGNOSTICS(result, session)
%
% Four panels, arranged so the argument of the tutorial reads left to right:
%
%   A  a stretch of the filtered trace with the threshold and detected spikes,
%      so the detection can be judged by eye
%   B  positive against negative peak heights. This is the panel that carries the
%      idea: the negative peaks are the noise, measured from this same trace, and
%      the threshold is simply a chosen quantile of them
%   C  the mixture fit over the detected heights, with both components drawn
%      separately, showing which part of the distribution is contamination
%   D  per-spike false-positive probability against height, with the false-
%      negative mass shaded below threshold
%
% Plotting only. All quantities come from detectSpikesMatchedFilter and
% fitSpikeHeightMixture; nothing is computed here that affects results.
%
% MATLAB R2019b compatible.

detection = result.detection;
fit       = result.fit;
dt        = session.framePeriod;

figureHandle = figure('Name', sprintf('Tutorial 2, cell %d', result.cellIndex), ...
    'Color', 'w', 'Position', [80 80 1180 760]);

%% --- A: filtered trace with threshold ------------------------------------
subplot(2, 2, 1);
% A window with some spikes in it, rather than the first samples, which may be quiet.
if isempty(detection.spikes)
    windowStart = 1;
else
    middleSpike = detection.spikes(max(1, round(numel(detection.spikes)/2)));
    windowStart = max(1, middleSpike - round(2/dt));
end
windowIdx = windowStart:min(windowStart + round(4/dt), numel(detection.filtered));
plot(windowIdx*dt, detection.filtered(windowIdx), 'Color', [.35 .35 .35]); hold on
inWindow = detection.spikes(detection.spikes >= windowIdx(1) & detection.spikes <= windowIdx(end));
plot(inWindow*dt, detection.filtered(inWindow), 'r.', 'MarkerSize', 9);
yline_compat(detection.threshold, 'r-', 'threshold');
yline_compat(-detection.threshold, 'b--', '-threshold');
xlabel('time (s)'); ylabel('matched-filter output');
title(sprintf('A  cell %d: %d spikes detected', result.cellIndex, numel(detection.spikes)));
axis tight; box off

%% --- B: the calibration idea ---------------------------------------------
subplot(2, 2, 2);
negative = sort(detection.negativePeakValues, 'descend');
positive = sort(detection.filtered(detection.spikes), 'descend');
plot(1:numel(positive), positive, 'k-', 'LineWidth', 1.4); hold on
plot(1:numel(negative), negative, 'b-', 'LineWidth', 1.4);
plot(numel(negative), detection.threshold, 'ro', 'MarkerFaceColor', 'r');
set(gca, 'XScale', 'log');
xlabel('rank'); ylabel('peak height');
legend({'positive peaks (spikes + noise)', 'negative peaks (noise only)', ...
    sprintf('threshold at rank %d', detection.expectedFalsePositives)}, ...
    'Location', 'northeast', 'Box', 'off');
title('B  threshold is a chosen quantile of the cell''s own noise');
box off

%% --- C: the mixture fit ---------------------------------------------------
subplot(2, 2, 3);
values = fit.analysedValues;
edges = linspace(fit.threshold, prctile(values, 99.5), 60);
histogram(values, edges, 'Normalization', 'pdf', ...
    'FaceColor', [.75 .75 .75], 'EdgeColor', 'none'); hold on
grid = linspace(fit.threshold, edges(end), 400);
componentFP = fit.mixingWeight * truncNorm(grid, fit.falsePositiveParams, fit.threshold);
componentTP = (1-fit.mixingWeight) * truncNorm(grid, fit.trueSpikeParams, fit.threshold);
plot(grid, componentFP, 'b-', 'LineWidth', 1.6);
plot(grid, componentTP, 'g-', 'LineWidth', 1.6);
plot(grid, componentFP + componentTP, 'k--', 'LineWidth', 1.4);
xlabel('matched-filter spike height'); ylabel('density');
legend({'detected heights', ...
    sprintf('false positives (p = %.3f)', fit.mixingWeight), ...
    'true spikes', 'mixture'}, 'Location', 'northeast', 'Box', 'off');
title('C  detected heights are contamination plus signal');
box off

%% --- D: error rates -------------------------------------------------------
subplot(2, 2, 4);
% The false-negative mass: the fitted true-spike component continued below the
% threshold, which is the part of the real distribution detection cannot see.
below = linspace(max(0, fit.trueSpikeParams(1) - 4*fit.trueSpikeParams(2)), ...
    fit.threshold, 200);
untruncated = normpdf(below, fit.trueSpikeParams(1), fit.trueSpikeParams(2));
area(below, untruncated, 'FaceColor', [1 .8 .8], 'EdgeColor', 'none'); hold on
above = linspace(fit.threshold, prctile(values, 99.5), 200);
plot(above, normpdf(above, fit.trueSpikeParams(1), fit.trueSpikeParams(2)), ...
    'g-', 'LineWidth', 1.6);
yyaxis right
plot(values, fit.falsePositiveProbability, '.', 'Color', [.2 .2 .8], 'MarkerSize', 6);
ylabel('P(false positive)'); ylim([0 1]);
yyaxis left
ylabel('true-spike density');
yline_compat(NaN, 'k-', '');   % keep the left axis active for the xline below
xline_compat(fit.threshold, 'r-', 'threshold');
xlabel('matched-filter spike height');
title(sprintf('D  false negatives %.1f%% (shaded), false positives %d expected', ...
    100*fit.falseNegativeRate, detection.expectedFalsePositives));
box off
end

%% ========================================================= local functions

function y = truncNorm(x, params, lowerLimit)
y = normpdf(x, params(1), params(2)) ./ (1 - normcdf(lowerLimit, params(1), params(2)));
end

function yline_compat(value, style, label)
% yline was introduced in R2018b, but its text placement changed later; drawing
% the line explicitly keeps the figure identical across releases.
if isnan(value), return; end
limits = get(gca, 'XLim');
plot(limits, [value value], style, 'HandleVisibility', 'off');
if ~isempty(label)
    text(limits(1), value, [' ' label], 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'Color', style(1));
end
end

function xline_compat(value, style, label)
limits = get(gca, 'YLim');
plot([value value], limits, style, 'HandleVisibility', 'off');
if ~isempty(label)
    text(value, limits(2), [' ' label], 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'Color', style(1));
end
end
