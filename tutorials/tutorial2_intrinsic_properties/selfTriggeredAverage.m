function [sta, nUsed] = selfTriggeredAverage(data, triggerFrames, timeWindow, varargin)
%SELFTRIGGEREDAVERAGE Per-cell average of a cell's own trace at its own events.
%
%   [sta, nUsed] = SELFTRIGGEREDAVERAGE(data, triggerFrames, timeWindow)
%   returns sta, nCells x (nBack+nFront+1), where row ii is the mean of that
%   cell's normalised voltage in the window timeWindow = [nBack nFront] around
%   each frame in triggerFrames{ii}. nUsed counts the triggers that contributed.
%
%   SELFTRIGGEREDAVERAGE(..., 'Name', value) accepts:
%
%     'Cells'    subset of cell indices, default all.
%     'Verbose'  true (default).
%
%   This is the streaming equivalent of the source's
%
%       get_sta_mat_self(triggers, traces_all', [nBack nFront])
%
%   restricted to the self-triggered case, which is all the intrinsic
%   properties need. The source builds a nSpike x nWindow index matrix against
%   a trace held whole in memory; that would be a per-cell read here, which
%   this asset's chunking punishes by 320x. Instead this makes one pass over
%   the acquisition chunks and accumulates.
%
%   WINDOW CONTAINMENT
%   Every trigger this tutorial uses is stimulus-locked, and the stimulus
%   occupies frames 394 to 7480 of each 8000-frame chunk, so both the +/-140
%   spike window and the -140/+390 blue-off window stay inside the chunk that
%   holds the trigger. The function asserts containment rather than stitching
%   across a chunk boundary, because a window that crossed one would span the
%   recording artifact the session notes flag there.
%
%   NORMALISATION
%   data.readChunk('voltage', ...) is already spike-height normalised, and a
%   triggered average is linear, so
%
%       STA(traces_all_n) == STA(traces_all) / spkHgtNoBlue == selfBlueStaFN
%
%   exactly. No spike-height variable appears anywhere in this tutorial.
%
%   MATCHING get_sta_mat_self
%   Two behaviours of the source helper are reproduced: a cell with fewer than
%   two triggers gets a row of zeros rather than NaN, and the mean ignores NaN
%   samples (nanmean). The doAlign branch is not ported; no intrinsic property
%   uses it.
%
%   See also STACKSTIMTRIALS, FITMEMBRANETIMECONSTANT.

opt = struct('Cells', [], 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
if isempty(opt.Cells)
    cells = (1:data.nCells)';
else
    cells = reshape(opt.Cells, [], 1);
end

nBack = timeWindow(1);
nFront = timeWindow(2);
nWindow = nBack + nFront + 1;
nFrame = data.nFrame;
timer = tic;

sta = zeros(data.nCells, nWindow);
acc = zeros(data.nCells, nWindow);
nUsed = zeros(data.nCells, 1);

% Bucket triggers by chunk once.
chunkOf = cell(data.nCells, 1);
for ii = cells'
    chunkOf{ii} = floor((triggerFrames{ii} - 1) / nFrame) + 1;
end

% get_sta_mat_self tests the cell's TOTAL trigger count, not the count in any
% one chunk, so the test has to be made before the streaming loop. A cell with
% one trigger in each of two chunks is included by the source.
nTrigTotal = zeros(data.nCells, 1);
for ii = cells'
    nTrigTotal(ii) = numel(triggerFrames{ii});
end
eligible = nTrigTotal >= 2;

touched = false(1, data.nChunk);
for ii = cells'
    if ~eligible(ii); continue; end
    touched(unique(chunkOf{ii})) = true;
end

for k = find(touched)
    block = [];
    for ii = cells'
        if ~eligible(ii); continue; end
        trig = triggerFrames{ii}(chunkOf{ii} == k);
        if isempty(trig); continue; end
        local = trig(:)' - (k-1) * nFrame;

        inside = local > nBack & local <= nFrame - nFront;
        if any(~inside)
            assert(false, 'Tutorial2:WindowCrossesChunk', ...
                ['Cell %d has a trigger at local frame %d, whose [%d %d] ' ...
                 'window leaves chunk %d. The stimulus should keep every ' ...
                 'window inside its chunk; check the protocol.'], ...
                ii, local(find(~inside, 1)), nBack, nFront, k);
        end
        if isempty(block)
            block = double(data.readChunk('voltage', k));
        end
        idx = local' + (-nBack:nFront);
        row = block(ii, :);
        window = row(idx);
        % get_sta_mat_self averages with nanmean, so a NaN sample would be
        % dropped per lag rather than poisoning the row. Neither deposited
        % session has a single non-finite sample, so this accumulates with a
        % plain sum and asserts the precondition instead of paying for
        % nanmean on every window.
        assert(all(isfinite(window(:))), 'Tutorial2:NonFiniteTrace', ...
            ['Cell %d chunk %d has non-finite samples in a trigger window. ' ...
             'get_sta_mat_self would nanmean past these; this function does ' ...
             'not, so the row would be wrong. Switch to nanmean before ' ...
             'trusting any STA from this asset.'], ii, k);
        acc(ii, :) = acc(ii, :) + sum(window, 1);
        nUsed(ii) = nUsed(ii) + numel(local);
    end
end

use = nUsed > 0;
sta(use, :) = acc(use, :) ./ nUsed(use);

if opt.Verbose
    fprintf('Self-triggered average over %d cells, window [-%d %d]\n', ...
        numel(cells), nBack, nFront);
    fprintf('  triggers/cell: median %g [%d, %d]; %d cells left at zero\n', ...
        median(nUsed(cells)), min(nUsed(cells)), max(nUsed(cells)), ...
        sum(nUsed(cells) == 0));
    fprintf('  %.0f s\n', toc(timer));
end
end
