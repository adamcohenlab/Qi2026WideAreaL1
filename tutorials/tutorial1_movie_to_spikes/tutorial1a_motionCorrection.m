function results = tutorial1a_motionCorrection(nwbFile, startSample, nFrames, doPlot)
% TUTORIAL 1a  Optical-flow motion correction of a wide-field voltage movie.
%
% results = TUTORIAL1A_MOTIONCORRECTION(nwbFile, startSample, nFrames, doPlot)
%
% Registers raw frames onto a still-period template by dense optical flow, and
% checks the result against the motion traces the published pipeline produced.
%
% The problem
% -----------
% The animal is awake on a treadmill. When it runs, the brain moves under the
% cranial window. A cell's fluorescence then changes for two entirely different
% reasons -- it fired, or it moved -- and at 787 Hz over a 2.9 x 2.1 mm field
% those are not separable after the fact. Everything downstream depends on
% fixing this first.
%
% Two outputs, and the second matters as much as the first
% --------------------------------------------------------
% Registration produces a corrected movie, which is what the demixing consumes.
% It also produces a dense displacement field, and averaging that field over one
% cell's region gives that cell's own motion trace. Those per-cell traces are
% used by every later stage: to choose which frames are trustworthy, and as the
% regressors that remove residual motion artifact from the extracted traces. The
% flow field is not a by-product of the correction, it is half the point.
%
% Why dense flow rather than a rigid shift
% ----------------------------------------
% Brain motion under a window is not a rigid translation of the field of view.
% The tissue deforms, so a cell at one edge can move differently from a cell at
% the other, and one (dx, dy) per frame cannot represent that. It is also why
% per-cell traces are meaningful: if the motion were rigid, every cell would
% have the same trace.
%
% What this tutorial does
% -----------------------
%   1. Rebuilds the registration template from running speed and the raw movie.
%      The template is not deposited -- the source pipeline never saved it --
%      but it is exactly reproducible. See buildRegistrationTemplate.
%   2. Picks a window where the animal was running, because a motion-correction
%      demonstration on a still window shows nothing.
%   3. Estimates dense flow against the template and warps the frames back.
%   4. Determines the warp sign empirically rather than assuming it (see below).
%   5. Derives per-cell motion traces from the flow field and correlates them
%      against the deposited ones, which is the real check that this reproduces
%      the published registration.
%
% The warp-sign trap
% ------------------
% The source has two implementations and they use OPPOSITE signs:
% `optical_flow_motion_correction` does imwarp(..., +flow) and
% `optical_flow_motion_correction_gpu` does imwarp(..., -flow). Both are correct,
% because MATLAB's and OpenCV's flow conventions are opposite. But substituting
% one for the other without flipping DOUBLES the motion instead of removing it,
% and nothing errors -- you get a corrected-looking movie that is twice as bad.
% This tutorial measures the sign instead of trusting it, and shows the
% measurement.
%
% What this reproduces, and what it does not
% ------------------------------------------
% The published data came from the GPU path: a compiled MEX over OpenCV's CUDA
% Farneback, which is not distributable. This uses MATLAB's own
% `opticalFlowFarneback`. Same algorithm family, different implementation, and
% the two also differ in pre-smoothing and in frame tiling. Expect the per-cell
% traces to track the deposited ones closely and not to match them exactly. The
% deposited traces are the ones the published analysis used; recompute here to
% understand the method, then use the deposited ones for anything quantitative.
%
% A timebase warning
% ------------------
% Do not convert seconds to sample indices by dividing by the frame period. The
% clock is not uniform: each acquisition-chunk boundary advances 20 microseconds
% rather than a full 1.27 ms frame period, because the upstream concatenation
% discarded the physical pause between chunks. Dividing accumulates about three
% samples of error across the window. Go through the timestamp vector:
%     idx = interp1(timestamps, 1:numel(timestamps), seconds, 'nearest');
% Elapsed time across a boundary is likewise not physical elapsed time. Treat
% the chunks as separate epochs; see /intervals/recording_chunks.
%
% Inputs
% ------
% nwbFile : char
%     The raw-movie excerpt asset.
% startSample : numeric, optional
%     First sample of the window. Default [] selects the highest-running-speed
%     stretch, which is where there is motion to correct.
% nFrames : numeric, optional
%     Window length. Default 1500 (about 1.9 s). This is memory-bound rather
%     than time-bound: the flow field alone is nRow*nCol*nFrames*2 singles,
%     1.8 GB at the default, and the raw, smoothed and corrected copies roughly
%     double that.
% doPlot : logical, optional. Default true.
%
% Outputs
% -------
% results : struct with
%     .template, .templateInfo    the rebuilt reference and how it was found
%     .window                     .start, .stop, .nFrames, .meanSpeed
%     .raw, .corrected            [nRow x nCol x nFrames]
%     .flow                       [nRow x nCol x nFrames x 2]
%     .mc                         info from opticalFlowMotionCorrection
%     .cellTraces                 recomputed per-cell motion, both region rules
%     .comparison                 correlation against the deposited traces
%
% Requires Computer Vision, Image Processing and Signal Processing Toolboxes.
% MATLAB R2019b compatible. MatNWB is not required.

if nargin < 2, startSample = []; end
if nargin < 3 || isempty(nFrames), nFrames = 1500; end
if nargin < 4 || isempty(doPlot),  doPlot  = true; end
assert(isfile(nwbFile), 'Tutorial1a:MissingFile', 'NWB file not found: %s', nwbFile);
addpath(fileparts(mfilename('fullpath')));

MOVIE  = '/acquisition/raw_voltage_movie/data';
SPEED  = '/processing/behavior/behavioral_time_series/running_speed/data';
MOTX   = '/processing/ophys/motion_trace_x/data';
MOTY   = '/processing/ophys/motion_trace_y/data';
SEG    = '/processing/ophys/ImageSegmentation/plane_segmentation/initial_image_mask';

%% ===================================================================== setup
% MATLAB sees HDF5 dimensions reversed, so the movie reads as [column row time].
movieDims = h5info(nwbFile, MOVIE).Dataspace.Size;
nCol = movieDims(1); nRow = movieDims(2); nSampleTotal = movieDims(3);
timestamps = h5read(nwbFile, '/acquisition/raw_voltage_movie/timestamps');
framePeriod = median(diff(timestamps));
fprintf('movie %dx%d, %d samples, frame period %.4f ms\n', ...
    nRow, nCol, nSampleTotal, 1e3*framePeriod);

speed = double(h5read(nwbFile, SPEED));
speed = speed(:)';

%% ============================================== 1. rebuild the template
fprintf('\n--- 1. registration template ---\n');
[template, templateInfo] = buildRegistrationTemplate(nwbFile);

%% ============================================== 2. choose a moving window
fprintf('\n--- 2. analysis window ---\n');
if isempty(startSample)
    % A motion-correction tutorial needs motion. Running speed is the cause of
    % brain motion here, and it is cheap to read in full, so it is used to
    % locate the window rather than the deposited motion traces, which are
    % 320 x 1.15M and expensive to scan.
    startSample = highestSpeedWindow(speed, nFrames, templateInfo);
    fprintf('selected the highest-running-speed window\n');
end
startSample = max(1, min(startSample, nSampleTotal - nFrames + 1));
window = startSample:(startSample + nFrames - 1);
fprintf('frames %d-%d, mean speed %.4g, peak %.4g\n', ...
    window(1), window(end), mean(speed(window)), max(speed(window)));
warnIfSpansChunkBoundary(nwbFile, window);

%% ============================================== 3. register
fprintf('\n--- 3. optical flow registration ---\n');
raw = h5read(nwbFile, MOVIE, [1 1 startSample], [nCol nRow nFrames]);
raw = permute(raw, [2 1 3]);                       % -> [row column time]
[corrected, flow, mc] = opticalFlowMotionCorrection(raw, template);

%% ============================================== 4. per-cell motion traces
% The source averages the flow field over each cell's DMD stimulation polygon.
% Those polygons are not deposited, but Stage 0 defines each initial footprint
% only inside its cell's polygon, so the support of `initial_image_mask` IS the
% polygon. Both that rule and the simpler rDisk = 10 disk are computed here, and
% correlating each against the deposited trace shows which the source used.
fprintf('\n--- 4. per-cell motion traces ---\n');
tutorialCells = find(h5read(nwbFile, '/units/tutorial_cell') == 1)';
useCells = tutorialCells(1:min(6, numel(tutorialCells)));
cellColumn = double(h5read(nwbFile, '/units/x_pixel'));   % x_pixel is the column
cellRow    = double(h5read(nwbFile, '/units/y_pixel'));

cellTraces = struct('cellIndex', {}, 'footprintXY', {}, 'diskXY', {}, ...
    'nFootprintPixels', {}, 'nDiskPixels', {});
for k = 1:numel(useCells)
    ci = useCells(k);
    plane = readFootprint(nwbFile, SEG, ci, nCol, nRow);
    footprintMask = plane > 0;
    diskMask = diskAround(cellRow(ci), cellColumn(ci), 10, nRow, nCol);

    cellTraces(k).cellIndex        = ci;
    cellTraces(k).footprintXY      = averageFlowOver(flow, footprintMask);
    cellTraces(k).diskXY           = averageFlowOver(flow, diskMask);
    cellTraces(k).nFootprintPixels = sum(footprintMask(:));
    cellTraces(k).nDiskPixels      = sum(diskMask(:));
end

%% ============================================== 5. compare with the deposited
fprintf('\n--- 5. agreement with the deposited motion traces ---\n');
comparison = compareWithDeposited(nwbFile, MOTX, MOTY, window, cellTraces);
printComparison(comparison);

%% ============================================== 6. against the corrected movie
% The step above checks the flow FIELD against the deposited per-cell traces.
% This checks the WARPED MOVIE against the deposited corrected pixels, which is
% the other half of what Stage 1 produced and a more direct test of the warp
% itself, sign included. Both targets are same-stage: mcTrace and movMCMasked
% are written in the same loop of the source script, with nothing in between.
fprintf('\n--- 6. agreement with the deposited corrected pixels ---\n');
pixelCheck = compareCorrectedPixels(nwbFile, raw, corrected, window, 400);
printPixelCheck(pixelCheck);

%% ===================================================================== output
results.template     = template;
results.templateInfo = templateInfo;
results.window       = struct('start', window(1), 'stop', window(end), ...
    'nFrames', nFrames, 'meanSpeed', mean(speed(window)), ...
    'peakSpeed', max(speed(window)), 'speed', speed(window), ...
    'fullSpeed', speed, 'timestamps', timestamps(window)');
results.raw          = raw;
results.corrected    = corrected;
results.flow         = flow;
results.mc           = mc;
results.cellTraces   = cellTraces;
results.comparison   = comparison;
results.pixelCheck   = pixelCheck;

if doPlot
    plotMotionCorrectionDiagnostics(results);
end
end

%% ========================================================= local functions

function startSample = highestSpeedWindow(speed, nFrames, templateInfo)
% Window whose mean running speed is highest, excluding the template's own
% frames: those were chosen for being still, so they are the worst possible
% demonstration and would also make the correction look artificially perfect.
n = numel(speed);
if nFrames >= n, startSample = 1; return; end
% Mean over a sliding window, as a difference of cumulative sums.
c = [0 cumsum(speed)];
windowMean = (c(nFrames+1:end) - c(1:end-nFrames)) / nFrames;
blocked = false(size(windowMean));
lo = max(1, templateInfo.startSample - nFrames + 1);
hi = min(numel(windowMean), templateInfo.stopSample);
blocked(lo:hi) = true;
windowMean(blocked) = -Inf;
[~, startSample] = max(windowMean);
end

function warnIfSpansChunkBoundary(nwbFile, window)
% Chunk boundaries carry a recording artifact and a non-physical time step, so a
% window that straddles one is a poor demonstration. Warn rather than refuse.
try
    startTimes = h5read(nwbFile, '/intervals/recording_chunks/start_time');
    nChunk = numel(startTimes);
catch
    return
end
chunkLength = 288000;
boundaries = chunkLength * (1:nChunk-1);
crossed = boundaries(boundaries > window(1) & boundaries < window(end));
if ~isempty(crossed)
    warning('Tutorial1a:SpansChunkBoundary', ...
        ['Window spans acquisition-chunk boundary at sample %s. Those frames ' ...
         'carry a recording artifact and the timestamps there advance 20 us ' ...
         'rather than a frame period. Consider a window inside one chunk.'], ...
        mat2str(crossed));
end
end

function plane = readFootprint(nwbFile, segPath, cellIndex, nCol, nRow)
% One cell's footprint as [row column]. Stored MATLAB order is [column row cell].
slab = h5read(nwbFile, segPath, [1 1 cellIndex], [nCol nRow 1]);
plane = permute(slab, [2 1 3]);
end

function mask = diskAround(centreRow, centreCol, radius, nRow, nCol)
[cols, rows] = meshgrid(1:nCol, 1:nRow);
mask = (rows - centreRow).^2 + (cols - centreCol).^2 <= radius^2;
end

function xy = averageFlowOver(flow, mask)
% Mean flow inside a region, per frame. Returns [2 x nFrames]: row 1 x, row 2 y.
[nRow, nCol, nFrame, ~] = size(flow);
idx = find(mask);
if isempty(idx), xy = zeros(2, nFrame, 'single'); return; end
flat = reshape(flow, nRow*nCol, nFrame, 2);
xy = [mean(flat(idx,:,1), 1); mean(flat(idx,:,2), 1)];
end

function comparison = compareWithDeposited(nwbFile, motXPath, motYPath, window, cellTraces)
% Correlate each recomputed trace against the deposited one for the same cell
% and the same samples. Deposited arrays are [cell x time] in MATLAB order, so
% only this window and these cells are read.
comparison = struct('cellIndex', {}, 'footprintR', {}, 'diskR', {}, ...
    'depositedRangeX', {}, 'recomputedRangeX', {});
for k = 1:numel(cellTraces)
    ci = cellTraces(k).cellIndex;
    depX = double(h5read(nwbFile, motXPath, [ci window(1)], [1 numel(window)]));
    depY = double(h5read(nwbFile, motYPath, [ci window(1)], [1 numel(window)]));

    comparison(k).cellIndex  = ci;
    comparison(k).footprintR = traceCorrelation(cellTraces(k).footprintXY, depX, depY);
    comparison(k).diskR      = traceCorrelation(cellTraces(k).diskXY,      depX, depY);
    comparison(k).depositedRangeX  = [min(depX) max(depX)];
    comparison(k).recomputedRangeX = [min(cellTraces(k).footprintXY(1,:)) ...
                                      max(cellTraces(k).footprintXY(1,:))];
    comparison(k).depositedX = depX;
    comparison(k).depositedY = depY;
end
end

function r = traceCorrelation(recomputedXY, depX, depY)
% One number per cell: the mean of the x and y correlations. Reported as a pair
% too, because a large difference between the axes would itself be informative.
rx = safeCorr(double(recomputedXY(1,:)), depX);
ry = safeCorr(double(recomputedXY(2,:)), depY);
r = struct('x', rx, 'y', ry, 'mean', mean([rx ry]));
end

function r = safeCorr(a, b)
a = a(:); b = b(:);
if std(a) == 0 || std(b) == 0, r = NaN; return; end
c = corrcoef(a, b);
r = c(1,2);
end

function check = compareCorrectedPixels(nwbFile, raw, corrected, window, nPixelsToCheck)
% Correlate this correction against the published one, pixel by pixel.
%
% `mc_pixel_subset` holds the source's own motion-corrected movie for the 16
% tutorial cells' pixel disks. Sampling the same pixels out of the corrected
% movie here gives a direct comparison of the two registrations. The raw movie
% is scored the same way as a baseline: if the correction is doing anything, it
% should agree with the published corrected pixels better than the uncorrected
% frames do.
%
% Only a sample of pixels is read. The subset is now chunked [1 nTime] for fast
% per-pixel access, so a time slab across all 5,712 pixels touches 5,712 chunks;
% a few hundred pixels is plenty for this statistic and costs a fraction of that.
SUBSET = '/processing/ophys/mc_pixel_subset/data';
MAP    = '/processing/ophys/mc_pixel_map';

check = struct('nPixels', 0, 'medianCorrectedR', NaN, 'medianRawR', NaN, ...
    'fractionImproved', NaN, 'correctedR', [], 'rawR', [], 'available', false);
if ~datasetExists(nwbFile, SUBSET)
    fprintf('mc_pixel_subset not present; skipping the pixel comparison\n');
    return
end

pixelRow    = double(h5read(nwbFile, [MAP '/pixel_row']));
pixelColumn = double(h5read(nwbFile, [MAP '/pixel_column']));
pixelUnitId = double(h5read(nwbFile, [MAP '/unit_id']));
nSubset = numel(pixelRow);

% Stratify by owning cell. The columns are 16 contiguous per-cell disks, so a
% sample spread evenly over the column index would land unevenly across cells
% and could end up testing one or two neighbourhoods of a 2.9 x 2.1 mm field.
% Sampling within each cell covers the whole field for the same number of reads.
units = unique(pixelUnitId);
perUnit = max(1, floor(min(nPixelsToCheck, nSubset) / numel(units)));
pick = [];
for u = 1:numel(units)
    columnsHere = find(pixelUnitId == units(u));
    take = unique(round(linspace(1, numel(columnsHere), ...
        min(perUnit, numel(columnsHere)))));
    pick = [pick; columnsHere(take)]; %#ok<AGROW>
end
pick = unique(pick);

nFrames = numel(window);
correctedR = nan(numel(pick), 1);
rawR       = nan(numel(pick), 1);
for k = 1:numel(pick)
    p = pick(k);
    deposited = double(h5read(nwbFile, SUBSET, [p window(1)], [1 nFrames]));
    mine = double(squeeze(corrected(pixelRow(p), pixelColumn(p), :)))';
    theirsRaw = double(squeeze(raw(pixelRow(p), pixelColumn(p), :)))';
    correctedR(k) = safeCorr(mine,      deposited);
    rawR(k)       = safeCorr(theirsRaw, deposited);
end

check.nPixels          = numel(pick);
check.correctedR       = correctedR;
check.rawR             = rawR;
check.medianCorrectedR = median(correctedR, 'omitnan');
check.medianRawR       = median(rawR, 'omitnan');
check.fractionImproved = mean(correctedR > rawR, 'omitnan');
check.available        = true;
end

function printPixelCheck(check)
if ~check.available, return; end
fprintf('%d pixels sampled from the deposited corrected subset\n', check.nPixels);
fprintf('  median r, this correction vs deposited : %.3f\n', check.medianCorrectedR);
fprintf('  median r, raw frames    vs deposited   : %.3f\n', check.medianRawR);
fprintf('  pixels where correcting helped         : %.1f%%\n', ...
    100*check.fractionImproved);
if check.medianCorrectedR <= check.medianRawR
    fprintf(['WARNING: correcting did not improve agreement with the published\n' ...
             'corrected pixels. Suspect the warp sign, the template, or the\n' ...
             'row/column convention before trusting anything downstream.\n']);
end
end

function tf = datasetExists(nwbFile, path)
try
    h5info(nwbFile, path);
    tf = true;
catch
    tf = false;
end
end

function printComparison(comparison)
fprintf('%-6s %-22s %-22s\n', 'cell', 'footprint support', 'rDisk = 10 disk');
fprintf('%-6s %-22s %-22s\n', '', 'r(x)   r(y)   mean', 'r(x)   r(y)   mean');
for k = 1:numel(comparison)
    f = comparison(k).footprintR; d = comparison(k).diskR;
    fprintf('%-6d %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f\n', ...
        comparison(k).cellIndex, f.x, f.y, f.mean, d.x, d.y, d.mean);
end
% Rank on magnitude, not signed correlation. The deposited traces were averaged
% from OpenCV's flow field, whose sign convention is opposite to MATLAB's -- the
% same difference that makes the source warp with -flow on the GPU path and
% +flow on the CPU path. So a near-perfect match appears here as r close to -1,
% and comparing signed values would rank the better region rule as the worse.
signedMean    = mean(arrayfun(@(c) c.footprintR.mean, comparison), 'omitnan');
footprintMean = mean(abs(arrayfun(@(c) c.footprintR.mean, comparison)), 'omitnan');
diskMean      = mean(abs(arrayfun(@(c) c.diskR.mean,      comparison)), 'omitnan');
fprintf('\nmean |r| over cells: footprint %.3f, disk %.3f\n', footprintMean, diskMean);

if signedMean < 0
    fprintf(['Correlations are systematically NEGATIVE. That is expected, not a ' ...
             'failure:\nthe deposited traces come from OpenCV''s flow field, ' ...
             'whose sign convention is\nopposite to MATLAB''s. It is the same ' ...
             'difference that makes the source warp\nwith -flow on the GPU path ' ...
             'and +flow on the CPU path. Magnitude is the result.\n']);
end

if footprintMean > diskMean
    fprintf(['The footprint support tracks the deposited traces better than the ' ...
             'disk does,\nwhich is what we expect if it recovers the DMD polygon ' ...
             'the source averaged over.\n']);
else
    fprintf(['The disk tracks better, so the footprint support is NOT recovering ' ...
             'the DMD\naveraging region. Worth investigating before relying on it.\n']);
end
end
