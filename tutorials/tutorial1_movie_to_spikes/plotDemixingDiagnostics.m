function figureHandle = plotDemixingDiagnostics(result, session)
% PLOTDEMIXINGDIAGNOSTICS  Diagnostic figure for one cell of tutorial 1.
%
% figureHandle = PLOTDEMIXINGDIAGNOSTICS(result, session)
%
% Four panels:
%
%   A  the initial footprints of the target and its neighbours, summed, showing
%      that they physically overlap. This is the problem statement.
%   B  the target's initial footprint beside the footprint the NMF learned, which
%      is what separates it from its neighbours
%   C  the naive and demixed traces averaged around every spike of the most
%      contaminating neighbour. A bump at lag zero in the naive trace and none in
%      the demixed one is the leakage being removed, over hundreds of events rather
%      than one hand-picked window
%   D  spike-band leakage per neighbour, before and after
%
% Plotting only. Nothing here affects results.
% MATLAB R2019b compatible.

demix = result.demix;
dt    = session.framePeriod;
pixelIdx  = demix.options.pixelIndices;
frameSize = demix.options.frameSize;

% A tight bounding box around the patch, so the footprints are visible rather than
% a few pixels in a 448 x 328 frame.
[rows, cols] = ind2sub(frameSize, pixelIdx);
rowRange = max(min(rows)-4, 1):min(max(rows)+4, frameSize(1));
colRange = max(min(cols)-4, 1):min(max(cols)+4, frameSize(2));

figureHandle = figure('Name', sprintf('Tutorial 1, cell %d', result.cellIndex), ...
    'Color', 'w', 'Position', [70 70 1200 780]);

%% --- A: the overlap problem ---------------------------------------------
subplot(2, 3, 1);
composite = zeros(frameSize);
composite(pixelIdx) = sum(demix.initialFootprints, 2);
imagesc(colRange, rowRange, composite(rowRange, colRange)); axis image; colormap(gca, hot);
hold on
plot(session.cellColumn(result.cellIndex), session.cellRow(result.cellIndex), ...
    'co', 'MarkerSize', 9, 'LineWidth', 1.4);
if ~isempty(result.neighbours)
    plot(session.cellColumn(result.neighbours), session.cellRow(result.neighbours), ...
        'c+', 'MarkerSize', 8, 'LineWidth', 1.2);
end
xlabel('column (px)'); ylabel('row (px)');
title(sprintf('A  cell %d and %d overlapping neighbours', ...
    result.cellIndex, numel(result.neighbours)));

%% --- B: initial versus learned footprint ---------------------------------
subplot(2, 3, 2);
initialImage = zeros(frameSize);
initialImage(pixelIdx) = demix.initialFootprints(:, 1);
imagesc(colRange, rowRange, initialImage(rowRange, colRange)); axis image
colormap(gca, hot); xlabel('column (px)');
title('B1  initial footprint (seed)');

subplot(2, 3, 3);
learnedImage = zeros(frameSize);
learnedImage(pixelIdx) = demix.footprint;
imagesc(colRange, rowRange, learnedImage(rowRange, colRange)); axis image
colormap(gca, hot); xlabel('column (px)');
title('B2  footprint learned by NMF');

%% --- C: spike-triggered average on the neighbour's spikes ----------------
% Averaging beats eyeballing here. A single window shows one or two of the
% neighbour's spikes and the reader has to take the rest on faith; triggering on
% every one of its spikes and averaging makes leakage unmissable and quantitative.
% If the neighbour's spikes leak into this cell, the naive trace has a
% spike-shaped bump at lag zero and the demixed trace does not.
%
% This is also the same construction the manuscript uses to measure cross-talk
% between cells, so the panel teaches a method the reader will meet again.
subplot(2, 3, [4 5]);
% A 20-sample (25 ms) high-pass, matching the spike band the leakage is measured
% in. A slower filter leaves subthreshold voltage in and the spikes disappear.
hp = @(x) x(:)' - movmedian(x(:)', 20);
% Each trace is divided by its OWN robust noise level before averaging. This is
% not cosmetic. The naive trace is a footprint-weighted sum of raw counts while the
% demixed trace comes out of the factorisation with its own scaling, so their
% absolute amplitudes are not comparable and plotting them on one axis unnormalised
% would let a scale difference masquerade as a result. In noise units the question
% the panel asks is well posed: relative to what each trace can resolve, how large
% is the deflection at the neighbour's spike?
robustSD = @(x) 1.4826 * median(abs(x - median(x)));
normalise = @(x) hp(x) / robustSD(hp(x));
naive = normalise(demix.naiveTrace);
mixed = normalise(demix.trace);
if ~isempty(result.leakage)
    [~, worst] = max(abs([result.leakage.naive]));
    nb = result.leakage(worst);
    nbTrace = nb.trace / (1.4826*median(abs(nb.trace - median(nb.trace))));
    [~, trigger] = findpeaks(nbTrace, 'MinPeakHeight', 4);

    halfWindow = round(0.05/dt);                 % +/- 50 ms
    trigger = trigger(trigger > halfWindow & trigger <= numel(naive) - halfWindow);
    lags = (-halfWindow:halfWindow) * dt * 1e3;  % ms
    window = trigger(:) + (-halfWindow:halfWindow);

    naiveSTA = mean(naive(window), 1);
    mixedSTA = mean(mixed(window), 1);
    naiveSTA = naiveSTA - mean(naiveSTA(1:20));  % baseline to the window start
    mixedSTA = mixedSTA - mean(mixedSTA(1:20));

    plot(lags, naiveSTA, 'Color', [.85 .35 .1], 'LineWidth', 1.6); hold on
    plot(lags, mixedSTA, 'Color', [.1 .4 .75], 'LineWidth', 1.6);
    plot([0 0], get(gca, 'YLim'), 'k:', 'HandleVisibility', 'off');
    legend({'naive mask average', 'demixed'}, 'Location', 'northwest', 'Box', 'off');
    xlabel(sprintf('lag from neighbour %d''s spikes (ms)', nb.neighbour));
    ylabel('robust SD');
    naivePeak = max(naiveSTA);
    mixedPeak = max(mixedSTA);
    title(sprintf(['C  lag-0 deflection %.1f -> %.1f noise SD after demixing ' ...
        '(n = %d spikes)'], naivePeak, mixedPeak, numel(trigger)));
    text(0.02, -0.30, ...
        {sprintf(['The residual %.1f SD is NOT necessarily leakage. These cells ' ...
                  'are electrically coupled, and a'], mixedPeak); ...
         'gap-junction spikelet also appears at lag 0 with a spike-like shape. ' ...
         'This measurement'; ...
         'cannot separate the two; that is what the spikelet analysis is for.'}, ...
        'Units', 'normalized', 'FontSize', 8, 'Color', [.4 .4 .4]);
    axis tight; box off
else
    text(0.5, 0.5, 'no overlapping neighbours', 'HorizontalAlignment', 'center');
    axis off
end

%% --- D: leakage before and after -----------------------------------------
subplot(2, 3, 6);
if isempty(result.leakage)
    text(0.5, 0.5, 'no overlapping neighbours', 'HorizontalAlignment', 'center');
    axis off
else
    values = [abs([result.leakage.naive]); abs([result.leakage.demixed])]';
    bar(values);
    set(gca, 'XTickLabel', arrayfun(@(n) num2str(n), [result.leakage.neighbour], ...
        'UniformOutput', false));
    xlabel('neighbour cell'); ylabel('|correlation|, spike band');
    legend({'naive', 'demixed'}, 'Location', 'northeast', 'Box', 'off');
    title('D  spike-band shared signal reduced');
    box off
    % The caveat belongs on the figure, not only in the text, because this panel is
    % what gets screenshotted out of context.
    text(0.02, -0.24, ...
        {'Spike band. Shared signal here is mostly optical leakage, but not'; ...
         'only: these cells are electrically coupled and spikelets live in'; ...
         'this band too. The change is the result; the residual is not a rate.'}, ...
        'Units', 'normalized', 'FontSize', 8, 'Color', [.4 .4 .4]);
end
end
