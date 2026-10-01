function results = tutorial1b_demixing(nwbFile, cellIndices, doPlot, nSamplesToUse)
% TUTORIAL 1b  Demixing overlapping cells in wide-field voltage imaging.
%
% results = TUTORIAL1B_DEMIXING(nwbFile, cellIndices, doPlot)
%
% Separates one cell's voltage trace from its overlapping neighbours, using one NWB
% file and nothing else.
%
% The problem
% -----------
% This field of view is 2.9 x 2.1 mm holding 320 layer 1 interneurons, so many
% somata sit within a point spread function of each other. The pixels around any
% one cell also carry its neighbours' light. Averaging them with a fixed spatial
% mask -- the obvious thing to do -- gives a trace containing the neighbours'
% spikes, and in a study about coupling between these cells that is fatal: leaked
% spikes look exactly like synaptic or electrical coupling.
%
% This tutorial shows the contaminated trace and the demixed trace side by side,
% and measures how much neighbour signal each one contains.
%
% Reading the numbers honestly
% ----------------------------
% Demixing does NOT drive neighbour correlation to zero, and it should not. Layer 1
% interneurons are synaptically and electrically coupled, which is the subject of
% the manuscript, so neighbouring cells are genuinely correlated. Nothing in this
% comparison separates residual leakage from real coupling. The defensible
% statement is the CHANGE in shared signal, not the absolute residual. A reader who
% tunes the method until neighbour correlation approaches zero is deleting biology.
%
% Inputs
% ------
% nwbFile : char
%     The raw-movie excerpt asset. Must have been through
%     patch_initial_footprints_to_nwb, because the demixing is seeded with the
%     initial segmentation footprints and the unpatched asset carries only the
%     refined ones. Seeding with the refined footprint would initialise the method
%     with its own answer.
% cellIndices : numeric, optional
%     One-based cell indices. Must be cells whose pixels are in the deposited
%     subset, i.e. flagged `tutorial_cell`. Default: the first two.
% doPlot : logical, optional. Default true.
%
% Outputs
% -------
% results : struct array with .cellIndex, .neighbours, .demix, .leakage
%
% Requires demixCellNMF and its dependencies, including CaImAn-MATLAB on the path.
% MATLAB R2019b compatible.

if nargin < 2, cellIndices = []; end
if nargin < 3 || isempty(doPlot), doPlot = true; end
% Default to the whole deposited window, which is what the published run consumed
% per cell. This was briefly restricted to one chunk when a cell took 4,860 s, but
% that was not the algorithm: the asset had been chunked across the pixel axis,
% amplifying every single-pixel read by 5717x. With the asset chunked one pixel per
% chunk the full window takes about 60 s per cell, so there is nothing to trade
% away. Pass a sample count to shorten it anyway.
if nargin < 4 || isempty(nSamplesToUse), nSamplesToUse = Inf; end
assert(isfile(nwbFile), 'Tutorial1b:MissingFile', 'NWB file not found: %s', nwbFile);
addpath(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

SEG = '/processing/ophys/ImageSegmentation/plane_segmentation';
assert(hasDataset(nwbFile, SEG, 'initial_image_mask'), 'Tutorial1b:NotPatched', ...
    ['This asset has no initial_image_mask. Run ' ...
     'patch_initial_footprints_to_nwb first: the demixing must be seeded with the ' ...
     'initial segmentation footprints, not the refined ones it produces.']);

%% ===================================================================== setup
session = readDemixSession(nwbFile);
fprintf('%d cells, %d samples, frame %dx%d, %d subset pixels\n', ...
    session.nCells, session.nSamples, session.frameSize(1), session.frameSize(2), ...
    numel(session.pixelRow));

% Restrict the analysis window. Everything derived from the timebase is trimmed
% together so the frame sets stay consistent with the shortened patch.
nUse = min(nSamplesToUse, session.nSamples);
if nUse < session.nSamples
    fprintf(['using the first %d of %d samples (%.1f s of %.1f s); pass Inf for ' ...
             'the full window\n'], nUse, session.nSamples, ...
        nUse*session.framePeriod, session.nSamples*session.framePeriod);
    session.nSamples        = nUse;
    session.timestamps      = session.timestamps(1:nUse);
    session.motionFreeFrames = session.motionFreeFrames(session.motionFreeFrames <= nUse);
    session.badFrames        = session.badFrames(session.badFrames <= nUse);
end

if isempty(cellIndices)
    cellIndices = session.tutorialCells(1:min(2, numel(session.tutorialCells)));
end
cellIndices = cellIndices(:)';

% The overlap rule: a cell counts as a neighbour of the target when its footprint
% exceeds its own significance level anywhere in the target's pixel disk. The level
% is set by area rather than by amplitude, so it means the same thing for a bright
% cell and a dim one.
ROI_HALF = 12;
AREA_LIMIT = 81;                      % pixels, the footprint area treated as "the cell"
threshPercentile = 100 * (1 - AREA_LIMIT*pi/4/(2*ROI_HALF+1)^2);

results = struct('cellIndex', {}, 'neighbours', {}, 'demix', {}, 'leakage', {});
for k = 1:numel(cellIndices)
    ci = cellIndices(k);
    fprintf('\n=== cell %d ===\n', ci);

    isMine = session.pixelUnitId == (ci - 1);
    assert(any(isMine), 'Tutorial1b:CellNotDeposited', ...
        ['Cell %d has no pixels in the deposited subset. Choose a cell flagged ' ...
         'tutorial_cell; see /units/tutorial_cell.'], ci);
    pixelRows = find(isMine);
    linearIdx = sub2ind(session.frameSize, session.pixelRow(pixelRows), ...
        session.pixelColumn(pixelRows));
    fprintf('patch: %d pixels\n', numel(pixelRows));

    % Footprints, smoothed as the source pipeline smooths them before use.
    neighbours = findNeighbours(nwbFile, SEG, session, ci, linearIdx, ...
        ROI_HALF, threshPercentile);
    fprintf('overlapping neighbours: %s\n', mat2str(neighbours));
    components = [ci neighbours];
    footprints = zeros(numel(linearIdx), numel(components));
    for j = 1:numel(components)
        plane = readFootprint(nwbFile, SEG, 'initial_image_mask', components(j), session);
        smoothed = imgaussfilt(plane, 1);
        footprints(:, j) = smoothed(linearIdx);
    end

    % The movie patch: motion-corrected, not raw.
    patch = readPatch(nwbFile, pixelRows, session.nSamples);

    opt = struct('pixelIndices', linearIdx, 'frameSize', session.frameSize, ...
        'motionFreeFrames', session.motionFreeFrames, 'badFrames', session.badFrames, ...
        'chunkLength', session.chunkLength);
    demix = demixCellNMF(patch, footprints, opt);
    fprintf('rejected pixels: %d\n', numel(demix.rejectedPixels));

    leakage = measureLeakage(nwbFile, session, demix, neighbours);
    printLeakage(leakage, neighbours);

    results(k).cellIndex  = ci; %#ok<AGROW>
    results(k).neighbours = neighbours;
    results(k).demix      = demix;
    results(k).leakage    = leakage;
end

if doPlot
    for k = 1:numel(results)
        plotDemixingDiagnostics(results(k), session);
    end
end
end

%% ========================================================= local functions

function session = readDemixSession(nwbFile)
dims = h5info(nwbFile, '/processing/ophys/mc_pixel_subset/data').Dataspace.Size;
session.nSamples = dims(2);
maskDims = h5info(nwbFile, ...
    '/processing/ophys/ImageSegmentation/plane_segmentation/image_mask').Dataspace.Size;
% MATLAB order is [column row cell]; the frame is [row column].
session.frameSize = [maskDims(2) maskDims(1)];
session.nCells    = maskDims(3);

session.pixelRow    = double(h5read(nwbFile, '/processing/ophys/mc_pixel_map/pixel_row'));
session.pixelColumn = double(h5read(nwbFile, '/processing/ophys/mc_pixel_map/pixel_column'));
session.pixelUnitId = double(h5read(nwbFile, '/processing/ophys/mc_pixel_map/unit_id'));
session.tutorialCells = find(h5read(nwbFile, '/units/tutorial_cell') == 1)';
% x_pixel is the column index and y_pixel the row index, as defined where the
% centroids are computed in the source segmentation script.
session.cellColumn = double(h5read(nwbFile, '/units/x_pixel'));
session.cellRow    = double(h5read(nwbFile, '/units/y_pixel'));

t = double(h5read(nwbFile, '/processing/ophys/mc_pixel_subset/timestamps'));
session.framePeriod = median(diff(t));
session.timestamps  = t(:)';
nChunk = numel(h5read(nwbFile, '/intervals/recording_chunks/id'));
session.chunkLength = session.nSamples / nChunk;

% Motion-free frames are recomputed from the deposited motion traces rather than
% read, so the tutorial shows where the frame set comes from. The multiplier is
% the demixing-stage value, which differs from the spike-finding one.
mx = h5read(nwbFile, '/processing/ophys/motion_trace_x/data');
my = h5read(nwbFile, '/processing/ophys/motion_trace_y/data');
selection = selectMotionFreeFrames(mx, my, struct('multiplier', 4));
session.motionFreeFrames = selection.frames;
session.motionSelection  = selection;

guard = -400:100;
bad = (0:session.chunkLength:session.nSamples)' + guard;
bad = unique(bad(:));
session.badFrames = bad(bad >= 1 & bad <= session.nSamples)';
end

function tf = hasDataset(nwbFile, groupPath, name)
info = h5info(nwbFile, groupPath);
tf = ~isempty(info.Datasets) && any(strcmp({info.Datasets.Name}, name));
end

function plane = readFootprint(nwbFile, segPath, column, cellIndex, session)
raw = h5read(nwbFile, [segPath '/' column], ...
    [1 1 cellIndex], [session.frameSize(2) session.frameSize(1) 1]);
plane = double(squeeze(raw))';          % back to [row column]
end

function neighbours = findNeighbours(nwbFile, segPath, session, ci, linearIdx, ...
    roiHalf, threshPercentile)
% A cell is a neighbour when its footprint clears its own area-based significance
% level somewhere in the target's patch.
%
% Candidates are prefiltered by centroid distance first. A cell's footprint is only
% significant within roiHalf of its own centroid, and the target's patch extends
% rDisk from the target's centroid, so anything further than their sum plus margin
% cannot possibly overlap. That turns 320 footprint reads into about ten, which is
% the difference between a tutorial that runs in seconds and one that does not.
CANDIDATE_RADIUS = 40;    % pixels; comfortably exceeds rDisk + roiHalf*sqrt(2)

distance = hypot(session.cellColumn - session.cellColumn(ci), ...
                 session.cellRow    - session.cellRow(ci));
candidates = find(distance <= CANDIDATE_RADIUS);
candidates = candidates(candidates ~= ci);

neighbours = [];
for other = candidates(:)'
    plane = readFootprint(nwbFile, segPath, 'initial_image_mask', other, session);
    % The significance level is computed over the neighbour's OWN box, centred on
    % its own centroid, because it describes that cell's footprint, not the target's.
    rows = min(max(round(session.cellRow(other))    + (-roiHalf:roiHalf), 1), session.frameSize(1));
    cols = min(max(round(session.cellColumn(other)) + (-roiHalf:roiHalf), 1), session.frameSize(2));
    level = prctile(reshape(plane(rows, cols), [], 1), threshPercentile);
    smoothed = imgaussfilt(plane, 1);
    % No positivity guard on `level`, deliberately. A level of exactly zero is the
    % common case, not a degenerate one: Stage 0 defines each footprint only inside
    % its DMD polygon, which occupies under 10% of the 25x25 box, so the 89.8th
    % percentile of the box is often 0. The test then reads "any positive footprint
    % weight on the target's disk", which is what the source pipeline does -- it
    % applies no such guard either. An earlier `level > 0` condition here rejected
    % every one of those cells and reported zero overlapping neighbours for a cell
    % that has five, which would have made this tutorial demonstrate the opposite
    % of its point.
    if any(smoothed(linearIdx) > level)
        neighbours(end+1) = other; %#ok<AGROW>
    end
end
end

function patch = readPatch(nwbFile, pixelRows, nSamples)
% Rows of mc_pixel_subset are pixels; the tutorial wants samples x pixels.
patch = zeros(nSamples, numel(pixelRows), 'single');
for k = 1:numel(pixelRows)
    patch(:, k) = h5read(nwbFile, '/processing/ophys/mc_pixel_subset/data', ...
        [pixelRows(k) 1], [1 nSamples]);
end
end

function leakage = measureLeakage(nwbFile, session, demix, neighbours)
% How much of each neighbour's own trace appears in the naive average and in the
% demixed trace, measured in two bands.
%
% Band choice is not cosmetic here, it decides what the number means.
%
%   SPIKE BAND, a 20-sample (25 ms) high-pass. Spikes are about five samples wide,
%       so this isolates them. Shared signal here is mostly optical leakage,
%       because a fixed mask sums light from both somata.
%
%   SLOW BAND, a 1000-sample (1.27 s) high-pass. Dominated by subthreshold
%       voltage, which these cells genuinely share.
%
% NEITHER BAND IS PURELY LEAKAGE. Layer 1 interneurons are synaptically and
% electrically coupled, which is the subject of the manuscript. A gap-junction
% spikelet is fast, a few milliseconds, and arrives at zero lag from the
% presynaptic spike -- the same band and the same latency as residual optical
% leakage. No correlation or spike-triggered average can separate the two, and the
% spikelet analysis exists precisely because that separation needs its own method.
%
% So the defensible quantity from this comparison is the CHANGE in shared signal,
% in both bands. The residual is not a leakage rate and must not be reported as
% one. The spike band is the more sensitive measure of what demixing does, because
% that is where most of the contamination sits; the slow band is reported alongside
% it because the two together show how much of the naive trace's shared signal was
% contamination rather than physiology.
fast = @(x) x(:) - movmedian(x(:), 20);
slow = @(x) x(:) - movmedian(x(:), 1000);
naiveFast = fast(demix.naiveTrace); naiveSlow = slow(demix.naiveTrace);
mixedFast = fast(demix.trace);      mixedSlow = slow(demix.trace);

leakage = struct('neighbour', {}, 'naive', {}, 'demixed', {}, ...
    'naiveSlow', {}, 'demixedSlow', {}, 'trace', {});
for j = 1:numel(neighbours)
    nb = neighbours(j);
    raw = double(h5read(nwbFile, ...
        '/processing/ophys/Fluorescence/roi_response_series/data', ...
        [nb 1], [1 session.nSamples]));
    leakage(j).neighbour   = nb; %#ok<AGROW>
    leakage(j).naive       = corr(naiveFast, fast(raw));
    leakage(j).demixed     = corr(mixedFast, fast(raw));
    leakage(j).naiveSlow   = corr(naiveSlow, slow(raw));
    leakage(j).demixedSlow = corr(mixedSlow, slow(raw));
    % The spike-band trace is retained so the diagnostic figure can show the
    % neighbour's own spikes beside the naive and demixed traces, which is where
    % the leakage is actually visible.
    leakage(j).trace       = fast(raw)';
end
end

function printLeakage(leakage, neighbours)
if isempty(leakage), fprintf('no overlapping neighbours to measure\n'); return; end
fprintf('%-11s %s %s\n', '', '   spike band (leakage)', '  slow band (coupling too)');
fprintf('%-11s %8s %8s %8s %8s %8s %8s\n', 'neighbour', ...
    'naive', 'demixed', 'reduced', 'naive', 'demixed', 'reduced');
for j = 1:numel(leakage)
    a = abs(leakage(j).naive);     b = abs(leakage(j).demixed);
    c = abs(leakage(j).naiveSlow); d = abs(leakage(j).demixedSlow);
    fprintf('cell %-6d %8.3f %8.3f %7.0f%% %8.3f %8.3f %7.0f%%\n', neighbours(j), ...
        leakage(j).naive, leakage(j).demixed, 100*(1 - b/max(a, eps)), ...
        leakage(j).naiveSlow, leakage(j).demixedSlow, 100*(1 - d/max(c, eps)));
end
fprintf(['\nRead the CHANGE, not the residual. Shared spike-band signal is mostly\n' ...
         'optical leakage, so its reduction is the clearest measure of what ' ...
         'demixing\ndoes. But these cells are electrically coupled: a gap-junction ' ...
         'spikelet is\nfast and arrives at zero lag, in the same band and at the ' ...
         'same latency as\nresidual leakage. Neither correlation nor a ' ...
         'spike-triggered average can\nseparate them, so the residual in either ' ...
         'band is not a leakage rate.\n']);
end
