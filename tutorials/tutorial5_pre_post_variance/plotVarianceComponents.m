function fig = plotVarianceComponents(summary, boot, titleText, yLabelText)
%PLOTVARIANCECOMPONENTS Bar chart of variance components and the pre-post bootstrap.
%
%   fig = PLOTVARIANCECOMPONENTS(summary, boot, titleText, yLabelText)
%   draws two panels:
%     left   presynaptic, postsynaptic and residual variance with their 95%
%            bootstrap intervals, and the pre-vs-post p-value as a bracket;
%     right  the bootstrap distribution of sigmaPre2 - sigmaPost2, with zero
%            and the basic 95% interval marked.
%
%   Error bars are drawn as explicit lines because a bootstrap interval need
%   not be symmetric about its estimate.

c = summary.components;
fig = figure('Name', titleText, 'Color', 'w', 'Position', [100 100 900 380]);

subplot(1, 2, 1);
bar(1:3, c.Estimate, 0.6, 'FaceColor', [0.55 0.65 0.8], 'EdgeColor', 'none');
hold on
for k = 1:3
    line([k k], [c.Lower95(k) c.Upper95(k)], 'Color', 'k', 'LineWidth', 1.5);
    line(k + [-0.1 0.1], [c.Lower95(k) c.Lower95(k)], 'Color', 'k', 'LineWidth', 1.5);
    line(k + [-0.1 0.1], [c.Upper95(k) c.Upper95(k)], 'Color', 'k', 'LineWidth', 1.5);
end
top = max(c.Upper95(1:2)) * 1.08;
line([1 1 2 2], [top * 0.97 top top top * 0.97], 'Color', 'k');
text(1.5, top, sprintf('p %s %.3g', pComparator(summary), summary.pTwoSided), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off
set(gca, 'XTick', 1:3, 'XTickLabel', {'Pre','Post','Residual'}, ...
    'Box', 'off', 'TickDir', 'out');
ylim([0 max(c.Upper95) * 1.25]);
ylabel(yLabelText);
title(titleText);

subplot(1, 2, 2);
valid = boot.replicates.deltaVariance(boot.replicates.isValid);
histogram(valid, 40, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
hold on
xline(0, 'k-', 'LineWidth', 1);
xline(summary.delta.Estimate, 'r-', 'LineWidth', 1.5);
xline(summary.delta.BasicLower95, 'r--');
xline(summary.delta.BasicUpper95, 'r--');
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('\sigma^2_{pre} - \sigma^2_{post}');
ylabel('Bootstrap replicates');
title(sprintf('%d of %d replicates beyond zero', summary.nBeyondZero, summary.nValid));
end

function s = pComparator(summary)
% When no replicate crosses zero the p-value is a bound, not a measurement.
if summary.nBeyondZero == 0
    s = '\leq';
else
    s = '=';
end
end
