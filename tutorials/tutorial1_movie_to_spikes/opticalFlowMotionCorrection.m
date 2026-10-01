function [corrected, flow, info] = opticalFlowMotionCorrection(movie, template, options)
% OPTICALFLOWMOTIONCORRECTION  Register frames to a template by dense optical flow.
%
% [corrected, flow, info] = OPTICALFLOWMOTIONCORRECTION(movie, template, options)
%
% Estimates a dense, per-pixel displacement field between a reference template and
% every frame, then warps each frame back onto the template.
%
% Why dense flow rather than rigid translation
% --------------------------------------------
% The field of view is 2.9 x 2.1 mm of cortex imaged through a cranial window.
% Brain motion under that window is not a rigid shift: the tissue deforms, so a
% cell at one edge can move in a different direction from a cell at the other.
% A single (dx, dy) per frame cannot describe that. Dense flow gives a
% displacement at every pixel, which is also what makes the per-cell motion trace
% possible -- averaging the field over one cell's region gives that cell's own
% displacement, and those traces are what the later stages regress against.
%
% Relationship to the published implementation
% --------------------------------------------
% The published data was produced by `optical_flow_motion_correction_gpu`, which
% calls a compiled MEX over OpenCV's CUDA Farneback. That is not distributable
% (it needs a MEX build, mexOpenCV and CUDA), so this uses MATLAB's own
% `opticalFlowFarneback`, following the CPU sibling
% `optical_flow_motion_correction`. Same algorithm family, different
% implementation, so results are close but not identical.
%
% Two deliberate departures from the CPU sibling, both documented in
% TUTORIAL1_VALIDATION.md:
%
%   1. Template and movie are scaled together, into [0, 1].
%
%      MATLAB's opticalFlowFarneback requires its input in [0, 1] (or uint8).
%      Give it anything larger and it returns EXACTLY ZERO FLOW, with no error
%      and no warning. Measured on this data: input at raw scale or divided by
%      256 both yield max|Vx| = 0.00000, while the same frames divided by 65535
%      yield 0.664. A synthetic 3-pixel shift is likewise recovered as zero at
%      the wrong scale.
%
%      This makes the CPU sibling `optical_flow_motion_correction` a silent
%      no-op when used the way Stage 1 uses it. That function normalises the
%      movie by its own maximum -- correct on its own -- but passes a supplied
%      template through unscaled, and Stage 1 supplies a template in raw camera
%      counts. Template out of range, flow identically zero, corrected movie
%      identical to the input. Do not copy that pattern.
%
%      The GPU path's /2^8 is not a counterexample: OpenCV's CUDA Farneback
%      does its own normalisation, so the same constant that works there
%      disables the MATLAB estimator.
%   2. The warp sign is determined empirically by default, rather than hardcoded.
%      The CPU and GPU siblings use OPPOSITE signs, because MATLAB's and OpenCV's
%      flow conventions differ. Getting it wrong doubles the motion instead of
%      removing it, and the failure is silent, so `warpSign = 'auto'` measures it
%      instead of assuming.
%
% Inputs
% ------
% movie : [nRow x nCol x nFrame] numeric
%     Raw frames. Not modified; the returned movie is warped from these, while
%     flow is estimated from a smoothed, rescaled copy.
% template : [nRow x nCol] numeric
%     Reference image the frames are warped onto. Same intensity scale as movie.
%
% options : struct, all optional
%     .intensityScale     divide movie and template by this before flow
%                         estimation. Default 'auto': intmax of the input class
%                         for integer data, otherwise the joint maximum. The
%                         result must land in [0, 1] or the estimator silently
%                         does nothing, so a numeric value you supply is
%                         checked and rejected if it does not.
%     .temporalSmoothing  'sgolay' (default, matching the GPU path), 'movmean'
%                         (the CPU sibling) or 'none'. Flow is estimated on the
%                         smoothed copy; the warp is applied to the raw frames.
%     .sgolayOrder        default 3
%     .sgolayFrames       default 9, must be odd
%     .movmeanFrames      default 10
%     .warpSign           'auto' (default), +1 or -1.
%     .interpolation      passed to imwarp. Default 'cubic'.
%     .zeroBorder         zero the one-pixel frame border after warping, which
%                         the source does because the warp leaves it undefined.
%                         Default true.
%     .verbose            default true.
%
% Outputs
% -------
% corrected : [nRow x nCol x nFrame] single
% flow      : [nRow x nCol x nFrame x 2] single, (:,:,:,1) = x, (:,:,:,2) = y,
%             in pixels, in the convention MATLAB's estimateFlow returns.
% info      : struct with .warpSign, .signTest, .residualBefore, .residualAfter,
%             .improvement, .elapsedSeconds, .options
%
% Requires Computer Vision Toolbox (opticalFlowFarneback), Image Processing
% Toolbox (imwarp) and, for the default smoothing, Signal Processing Toolbox
% (sgolayfilt). MATLAB R2019b compatible.

if nargin < 3, options = struct(); end
options = setDefault(options, 'intensityScale',    'auto');
options = setDefault(options, 'temporalSmoothing', 'sgolay');
options = setDefault(options, 'sgolayOrder',       3);
options = setDefault(options, 'sgolayFrames',      9);
options = setDefault(options, 'movmeanFrames',     10);
options = setDefault(options, 'warpSign',          'auto');
options = setDefault(options, 'interpolation',     'cubic');
options = setDefault(options, 'zeroBorder',        true);
options = setDefault(options, 'verbose',           true);

validateattributes(movie,    {'numeric'}, {'3d',   'nonempty'}, mfilename, 'movie');
validateattributes(template, {'numeric'}, {'2d',   'nonempty'}, mfilename, 'template');
[nRow, nCol, nFrame] = size(movie);
assert(isequal(size(template), [nRow nCol]), 'MotionCorrection:TemplateSize', ...
    'Template is %dx%d but frames are %dx%d.', size(template,1), size(template,2), nRow, nCol);
assert(mod(options.sgolayFrames,2) == 1, 'MotionCorrection:SgolayEven', ...
    'sgolayFrames must be odd; got %d.', options.sgolayFrames);

tStart = tic;

% ---------------------------------------------------------------- preparation
% Flow is estimated on a smoothed, rescaled copy and applied to the raw frames.
% Smoothing along time matters more than it looks: Farneback assumes brightness
% constancy, and at 787 Hz the shot noise between consecutive frames violates
% that badly enough to put spurious structure in the flow field.
if ischar(options.intensityScale) || isstring(options.intensityScale)
    assert(strcmpi(options.intensityScale, 'auto'), 'MotionCorrection:BadScale', ...
        'intensityScale must be ''auto'' or a positive number.');
    if isinteger(movie)
        scale = double(intmax(class(movie)));
    else
        scale = max(max(double(movie(:))), max(double(template(:))));
    end
else
    scale = double(options.intensityScale);
    assert(scale > 0, 'MotionCorrection:BadScale', 'intensityScale must be positive.');
end
forFlow  = single(movie) / scale;
templateScaled = single(template) / scale;

% Guard the silent-failure mode. opticalFlowFarneback returns exactly zero flow
% for input outside [0, 1], without erroring, so the whole correction becomes a
% no-op that still produces plausible-looking output. Refuse rather than return
% a movie that looks corrected and is not.
peak = max(max(forFlow(:)), max(templateScaled(:)));
assert(peak <= 1 + 1e-6, 'MotionCorrection:ScaledOutOfRange', ...
    ['After dividing by %g the data peaks at %.3f, outside the [0, 1] range ' ...
     'opticalFlowFarneback requires. It would return zero flow and report no ' ...
     'error. Use intensityScale = ''auto'', or pass at least %g.'], ...
    scale, peak, scale*peak);
if options.verbose
    fprintf('  intensity scale %g, scaled peak %.3f\n', scale, peak);
end

switch lower(options.temporalSmoothing)
    case 'sgolay'
        % Matches the GPU path. sgolayfilt works down columns, so the movie is
        % flattened to [pixel x time] and filtered along dimension 2.
        n = options.sgolayFrames;
        if nFrame >= n
            flat = reshape(forFlow, nRow*nCol, nFrame);
            flat = sgolayfilt(double(flat), options.sgolayOrder, n, ones(n,1), 2);
            forFlow = single(reshape(flat, nRow, nCol, nFrame));
        elseif options.verbose
            fprintf(['  %d frames is fewer than the %d-frame Savitzky-Golay ' ...
                     'window; smoothing skipped\n'], nFrame, n);
        end
    case 'movmean'
        forFlow = movmean(forFlow, options.movmeanFrames, 3);
    case 'none'
        % leave as is
    otherwise
        error('MotionCorrection:BadSmoothing', ...
            'temporalSmoothing must be sgolay, movmean or none; got ''%s''.', ...
            options.temporalSmoothing);
end

% ------------------------------------------------------------ flow estimation
% estimateFlow is stateful: it returns the flow between the previous call's
% frame and this one. Re-priming with the template before every frame therefore
% gives flow(template -> frame) each time, rather than frame-to-frame flow. That
% re-priming is the source's idiom and is load-bearing, not redundant.
flow = zeros(nRow, nCol, nFrame, 2, 'single');
opticFlow = opticalFlowFarneback;
reportEvery = max(1, round(nFrame/10));
for ii = 1:nFrame
    estimateFlow(opticFlow, templateScaled);
    thisFlow = estimateFlow(opticFlow, forFlow(:,:,ii));
    flow(:,:,ii,1) = thisFlow.Vx;
    flow(:,:,ii,2) = thisFlow.Vy;
    if options.verbose && mod(ii, reportEvery) == 0
        fprintf('  flow %d/%d\n', ii, nFrame);
    end
end

% A second guard on the same failure, after the fact. The range check above
% should make this unreachable, but zero flow is catastrophic and invisible
% downstream, so it is worth confirming rather than assuming.
assert(any(flow(:) ~= 0), 'MotionCorrection:ZeroFlow', ...
    ['The estimated flow field is identically zero across all %d frames. The ' ...
     'usual cause is input outside [0, 1]; the intensity scale used was %g. ' ...
     'A genuinely motionless recording still produces small nonzero flow.'], ...
    nFrame, scale);

% --------------------------------------------------------------- warp sign
% Which sign undoes the motion depends on the flow convention, and the two
% source implementations disagree. Rather than trust either, warp a sample of
% frames both ways and keep whichever moves them closer to the template.
signTest = struct('tested', false, 'residualPlus', NaN, 'residualMinus', NaN, ...
    'nFramesTested', 0);
if ischar(options.warpSign) || isstring(options.warpSign)
    assert(strcmpi(options.warpSign, 'auto'), 'MotionCorrection:BadWarpSign', ...
        'warpSign must be ''auto'', +1 or -1.');
    nTest = min(nFrame, 25);
    testIdx = round(linspace(1, nFrame, nTest));
    residualPlus  = warpResidual(movie, template, flow, testIdx, +1, options);
    residualMinus = warpResidual(movie, template, flow, testIdx, -1, options);
    if residualMinus < residualPlus, warpSign = -1; else, warpSign = +1; end
    signTest.tested        = true;
    signTest.residualPlus  = residualPlus;
    signTest.residualMinus = residualMinus;
    signTest.nFramesTested = nTest;
    if options.verbose
        fprintf(['  warp sign: residual %+d = %.4g, %+d = %.4g over %d frames ' ...
                 '-> using %+d\n'], 1, residualPlus, -1, residualMinus, nTest, warpSign);
    end
else
    warpSign = sign(options.warpSign);
    assert(warpSign ~= 0, 'MotionCorrection:BadWarpSign', 'warpSign must be nonzero.');
end

% ------------------------------------------------------------------- warping
corrected = zeros(nRow, nCol, nFrame, 'single');
for ii = 1:nFrame
    corrected(:,:,ii) = imwarp(single(movie(:,:,ii)), ...
        warpSign * squeeze(flow(:,:,ii,:)), options.interpolation);
end
if options.zeroBorder
    corrected([1 end], :, :) = 0;
    corrected(:, [1 end], :) = 0;
end

% ------------------------------------------------------------- what it achieved
% Mean absolute deviation from the template, before and after. Border pixels are
% excluded because zeroBorder sets them to zero, which would otherwise register
% as a large "residual" that the correction created rather than removed.
inner = false(nRow, nCol);
inner(2:end-1, 2:end-1) = true;
info.residualBefore = frameResidual(single(movie),  template, inner);
info.residualAfter  = frameResidual(corrected,      template, inner);
info.improvement    = 1 - info.residualAfter / info.residualBefore;
info.warpSign       = warpSign;
info.signTest       = signTest;
info.elapsedSeconds = toc(tStart);
info.options        = options;

if options.verbose
    fprintf('  residual to template: %.4g -> %.4g (%.1f%% reduction), %.1f s\n', ...
        info.residualBefore, info.residualAfter, 100*info.improvement, ...
        info.elapsedSeconds);
end
end

%% ========================================================= local functions

function r = warpResidual(movie, template, flow, idx, warpSign, options)
% Mean absolute deviation from the template after warping a sample of frames.
inner = false(size(template));
inner(2:end-1, 2:end-1) = true;
total = 0;
for k = 1:numel(idx)
    ii = idx(k);
    warped = imwarp(single(movie(:,:,ii)), ...
        warpSign * squeeze(flow(:,:,ii,:)), options.interpolation);
    total = total + frameResidual(warped, template, inner);
end
r = total / numel(idx);
end

function r = frameResidual(movie, template, inner)
% Mean |frame - template| over inner pixels, averaged across frames. Both are
% converted to double first: these are uint16 counts, and the subtraction would
% otherwise saturate at zero wherever the frame is darker than the template.
nFrame = size(movie, 3);
t = double(template(inner));
total = 0;
for ii = 1:nFrame
    f = double(movie(:,:,ii));
    total = total + mean(abs(f(inner) - t));
end
r = total / nFrame;
end

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = value; end
end
