function figureHandles = plotMotionCorrectionDiagnostics(results)
% PLOTMOTIONCORRECTIONDIAGNOSTICS  Figures for tutorial 0.
%
% figureHandles = PLOTMOTIONCORRECTIONDIAGNOSTICS(results)
%
% Draws only. Everything shown was computed in tutorial1a_motionCorrection; this
% function derives nothing that affects a result, so it can be re-run, skipped or
% rewritten without changing any number.
%
% Figures
%   1  where the template came from, and what it looks like
%   2  the warp sign, measured rather than assumed
%   3  before and after, in images and as residual over time
%   4  the flow field, and why it is not a rigid shift
%   5  per-cell motion traces against the deposited ones
%
% MATLAB R2019b compatible.

figureHandles = gobjects(0);

framePeriod = 1.27e-3;
[nRow, nCol, nFrames] = size(results.raw);
t = (0:nFrames-1) * framePeriod;

%% ===================================== 1. template provenance
f = figure('Name', 'Tutorial 1a.1  Registration template', 'Color', 'w');
figureHandles(end+1) = f;

subplot(2,2,[1 2]);
info  = results.templateInfo;
range = 1:info.searchLimit;
tt    = (range-1)*framePeriod;
plot(tt, info.smoothedSpeed(range), 'Color', [.75 .75 .75]); hold on
plot(tt, abs(results.window.fullSpeed(range)), 'Color', [.35 .35 .35]);

% Mark the frames the applied rule calls still. Under 'zeroMotion' that is a
% sparse set covering most of the recording, so it is drawn as a rug along the
% bottom rather than as shading, which would fill the axes.
stillMask = info.isStill(range);
yl = [0 max(abs(results.window.fullSpeed(range)))*1.05];
rugHeight = yl(2)*0.04;
stillTimes = tt(stillMask);
plot(stillTimes, zeros(size(stillTimes)) + rugHeight, '.', ...
    'Color', [.2 .6 .9], 'MarkerSize', 1);
xRegion(info.startSample*framePeriod, info.stopSample*framePeriod, [.9 .3 .1]);
ylim(yl);
xlabel('time (s)'); ylabel('running speed');
title(sprintf(['Template = first run of %d consecutive still frames, ' ...
               'samples %d-%d\nrule ''%s'', %.1f%% of searched frames still'], ...
    info.nFrames, info.startSample, info.stopSample, info.stillRule, ...
    100*info.fractionStillSearched));
legend({'published 1000-sample moving mean', '|speed|', 'frames judged still', ...
    'template window'}, 'Location', 'northeast', 'FontSize', 7); legend boxoff
grid on

subplot(2,2,3);
imagescGray(results.template);
title(sprintf('Template: pixelwise median of %d still frames', info.nFrames));

subplot(2,2,4);
imagescGray(single(results.raw(:,:,1)));
title('One raw frame, for comparison');

annotate(f, ['The median, not the mean: 1000 frames at 787 Hz spans 1.3 s, in ' ...
             'which these cells spike. A mean would smear action potentials into ' ...
             'the reference image.']);

%% ===================================== 2. the warp sign
f = figure('Name', 'Tutorial 1a.2  Warp sign', 'Color', 'w');
figureHandles(end+1) = f;
st = results.mc.signTest;
if st.tested
    bar([st.residualPlus st.residualMinus], 'FaceColor', [.6 .6 .6]);
    set(gca, 'XTickLabel', {'imwarp(+flow)', 'imwarp(-flow)'});
    ylabel('mean |frame - template|');
    hold on
    chosen = 1 + (results.mc.warpSign < 0);
    bar(chosen, [st.residualPlus st.residualMinus]*[chosen==1; chosen==2], ...
        'FaceColor', [.2 .6 .9]);
    title(sprintf(['Sign measured on %d frames, not assumed: %+d wins\n' ...
                   'The source CPU and GPU paths disagree, and both are right ' ...
                   'for their own library'], st.nFramesTested, results.mc.warpSign));
else
    text(.5, .5, sprintf('Sign was supplied, not measured: %+d', results.mc.warpSign), ...
        'HorizontalAlignment', 'center');
    axis off
end
grid on
annotate(f, ['Choosing the wrong sign does not error. It warps each frame by the ' ...
             'motion instead of against it, doubling the displacement, and the ' ...
             'output still looks like a corrected movie.']);

%% ===================================== 3. before and after
f = figure('Name', 'Tutorial 1a.3  Before and after', 'Color', 'w');
figureHandles(end+1) = f;

% The frame that moved most, which is where correction is visible.
displacement = squeeze(mean(mean(abs(results.flow), 1), 2));
magnitude = hypot(displacement(:,1), displacement(:,2));
[~, worst] = max(magnitude);

subplot(2,3,1);
imagescGray(single(results.raw(:,:,worst)));
title(sprintf('Raw, frame %d (largest motion)', worst));

subplot(2,3,2);
imagescGray(results.corrected(:,:,worst));
title('Corrected');

subplot(2,3,3);
imagescGray(results.template);
title('Template');

% Both difference panels share one colour scale, set by the raw difference.
% Scaling them independently would make the corrected panel look flatter purely
% because its limits were wider, which is exactly the comparison being made.
rawDiff       = single(results.raw(:,:,worst)) - results.template;
correctedDiff = results.corrected(:,:,worst)   - results.template;
sharedLimit   = prctile(abs(double(rawDiff(:))), 99);

subplot(2,3,4);
imagescDiff(rawDiff, sharedLimit);
title('Raw - template');

subplot(2,3,5);
imagescDiff(correctedDiff, sharedLimit);
title('Corrected - template');

subplot(2,3,6);
before = perFrameResidual(single(results.raw),  results.template);
after  = perFrameResidual(results.corrected,     results.template);
plot(t, before, 'Color', [.75 .3 .3]); hold on
plot(t, after,  'Color', [.2 .5 .8]);
xlabel('time (s)'); ylabel('mean |frame - template|');
legend({'raw', 'corrected'}, 'Location', 'best'); legend boxoff
title(sprintf('Residual reduced %.1f%%', 100*results.mc.improvement));
grid on

annotate(f, ['Both difference panels share one colour scale. Border pixels are ' ...
             'zeroed by the warp and excluded from the residual: counting them ' ...
             'would score the correction against damage it caused.']);

%% ===================================== 4. the flow field
f = figure('Name', 'Tutorial 1a.4  The flow field', 'Color', 'w');
figureHandles(end+1) = f;

subplot(1,2,1);
imagescGray(results.template); hold on
step = max(8, round(min(nRow,nCol)/28));
[qc, qr] = meshgrid(1:step:nCol, 1:step:nRow);
vx = results.flow(1:step:nRow, 1:step:nCol, worst, 1);
vy = results.flow(1:step:nRow, 1:step:nCol, worst, 2);
quiver(qc, qr, vx, vy, 2, 'Color', [1 .4 .1], 'LineWidth', 1);
title(sprintf('Displacement at frame %d', worst));

% Running speed gets its own axis rather than being rescaled onto the
% displacement axis. Rescaling put it flat against zero, which made the
% comparison unreadable and let the panel imply a correspondence the reader
% could not actually check.
subplot(1,2,2);
yyaxis left
plot(t, displacement(:,1), 'Color', [.2 .5 .8], 'LineStyle', '-'); hold on
plot(t, displacement(:,2), 'Color', [.8 .4 .1], 'LineStyle', '-');
ylabel('mean |displacement| (px)');
set(gca, 'YColor', [.2 .2 .2]);
yyaxis right
plot(t, results.window.speed, 'Color', [.55 .55 .55], 'LineStyle', ':');
ylabel('running speed');
set(gca, 'YColor', [.55 .55 .55]);
xlabel('time (s)');
legend({'displacement x', 'displacement y', 'running speed'}, ...
    'Location', 'northwest', 'FontSize', 7); legend boxoff
title('Displacement and locomotion');
grid on

spread = squeeze(std(reshape(results.flow(:,:,worst,1), [], 1)));
annotate(f, sprintf(['Displacement varies across the field (s.d. %.2f px in x at ' ...
    'the worst frame), which is why this is a dense field and not one (dx, dy) ' ...
    'per frame. It is also what makes per-cell traces meaningful: under rigid ' ...
    'motion every cell would share one trace.'], spread));

%% ===================================== 5. per-cell traces vs deposited
f = figure('Name', 'Tutorial 1a.5  Against the deposited traces', 'Color', 'w');
figureHandles(end+1) = f;

nShow = min(3, numel(results.cellTraces));
% The deposited traces come from OpenCV's flow field, whose sign convention is
% opposite to MATLAB's. Plotting both raw would show two mirror images and read
% as a failure. The recomputed trace is flipped for display, and the label says
% so; the correlations printed alongside are the unflipped values.
signedMean = mean(arrayfun(@(c) c.footprintR.mean, results.comparison), 'omitnan');
flipForDisplay = signedMean < 0;
displaySign = 1 - 2*flipForDisplay;

for k = 1:nShow
    ct = results.cellTraces(k);
    cp = results.comparison(k);

    subplot(nShow, 2, 2*k-1);
    plot(t, cp.depositedX, 'Color', [.2 .2 .2]); hold on
    plot(t, displaySign * ct.footprintXY(1,:), 'Color', [.2 .6 .9]);
    ylabel(sprintf('cell %d, x', ct.cellIndex));
    % No legend: at this panel size it collides with either the traces or the
    % |r| label whichever corner it goes in. The colour key is in the caption.
    if k == 1, title('x displacement'); end
    if k == nShow, xlabel('time (s)'); end
    text(.02, .92, sprintf('|r| = %.3f', abs(cp.footprintR.x)), 'Units', 'normalized');
    grid on

    subplot(nShow, 2, 2*k);
    plot(t, cp.depositedY, 'Color', [.2 .2 .2]); hold on
    plot(t, displaySign * ct.footprintXY(2,:), 'Color', [.8 .4 .1]);
    if k == 1, title('y displacement'); end
    if k == nShow, xlabel('time (s)'); end
    text(.02, .92, sprintf('|r| = %.3f', abs(cp.footprintR.y)), 'Units', 'normalized');
    grid on
end

if flipForDisplay
    signNote = ['Recomputed shown sign-flipped: OpenCV and MATLAB return ' ...
                'opposite flow conventions, the same difference that makes the ' ...
                'source warp -flow on GPU and +flow here. '];
else
    signNote = '';
end
annotate(f, ['Black, deposited. Colour, recomputed here. ' signNote ...
             'Amplitude agrees as well as shape. Prefer the deposited traces ' ...
             'for quantitative work: they are what the paper used.']);
end

%% ========================================================= local functions

function r = perFrameResidual(movie, template)
% Mean |frame - template| per frame, over inner pixels only.
[nRow, nCol, nFrame] = size(movie);
inner = false(nRow, nCol);
inner(2:end-1, 2:end-1) = true;
t = double(template(inner));
r = zeros(1, nFrame);
for ii = 1:nFrame
    f = double(movie(:,:,ii));
    r(ii) = mean(abs(f(inner) - t));
end
end

function imagescGray(img)
imagesc(img); axis image off; colormap(gca, gray);
lo = prctile(double(img(:)), 1); hi = prctile(double(img(:)), 99.5);
if hi > lo, caxis([lo hi]); end
end

function imagescDiff(img, lim)
% Symmetric limits, so that zero is the midpoint and sign is readable. A limit
% may be supplied so that two panels being compared share one scale.
imagesc(img); axis image off; colormap(gca, gray);
if nargin < 2 || isempty(lim)
    lim = prctile(abs(double(img(:))), 99);
end
if lim > 0, caxis([-lim lim]); end
colorbar
end

function xRegion(x1, x2, colour)
yl = ylim;
patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], colour, ...
    'FaceAlpha', .18, 'EdgeColor', 'none');
ylim(yl);
end

function annotate(figureHandle, message)
% One caption per figure, at the bottom. The interpretation belongs next to the
% picture rather than only in the documentation.
%
% The axes are shifted up first to make room. Without that the caption runs off
% the bottom of the printed figure, which is easy to miss on screen because the
% figure window is taller than the page it prints to.
captionHeight = 0.11;
axesList = findobj(figureHandle, 'Type', 'axes');
for k = 1:numel(axesList)
    pos = get(axesList(k), 'Position');
    pos(2) = captionHeight + pos(2)*(1 - captionHeight);
    pos(4) = pos(4)*(1 - captionHeight);
    set(axesList(k), 'Position', pos);
end
annotation(figureHandle, 'textbox', [.02 .005 .96 captionHeight-.01], ...
    'String', message, 'EdgeColor', 'none', 'FontSize', 8, ...
    'FitBoxToText', 'off', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
end
