function stacks = stackStimTrials(data, varargin)
%STACKSTIMTRIALS Trial-average every cell over its own stimulation repeats.
%
%   stacks = STACKSTIMTRIALS(data) makes one streaming pass over the asset and
%   returns the per-cell stimulus-locked averages the intrinsic properties are
%   read from:
%
%     fr_mat_stim_avg  nCells x nFrame, probability of a spike in each frame of
%                      the protocol, averaged over that cell's repeats.
%     traces_stim_avg  nCells x nFrame, trial-averaged SUBTHRESHOLD voltage in
%                      normalised (spike-height) units.
%     nTrial           nCells x 1, repeats contributing to each average.
%     spk_t_blue       nCells x 1 cell array of session-global frame indices of
%                      spikes that fell while that cell's own blue was on.
%
%   STACKSTIMTRIALS(data, 'Name', value, ...) accepts:
%
%     'Subthreshold'  'interpolate' (default) rebuilds the subthreshold with
%                     subthresholdByInterpolation, which is what the property
%                     definitions use. 'deposited' reads traces_all_sub_n
%                     instead, for the cross-check in tutorial2c. 'none' skips
%                     the voltage entirely and returns traces_stim_avg empty:
%                     the firing-rate properties need no trace, so tutorial2b
%                     runs in seconds off the spike lists already in memory.
%     'Cells'         subset of cell indices, default all. The I/O cost is
%                     unchanged: a time block is read whole either way.
%     'Verbose'       true (default).
%
%   WHAT THIS REPLACES
%   The source builds four nCells x nFrame x nTrial arrays (spk_t_mat_stim,
%   traces_all_stim, blue_all_stim, plus the two averages) by reshaping whole
%   session-length rows. Only the averages are ever used for a published
%   property, so this keeps the two averages and drops the three full stacks,
%   which is the difference between 500 MB and 20 MB. tutorial2a can still
%   materialise a single cell's stack for plotting.
%
%   I/O
%   One pass, reading each acquisition chunk once as a whole-population block.
%   Every series is chunked [nCells 10000] in MATLAB order, so this is the
%   cheap access pattern and a per-cell loop is a 320x amplification. See the
%   I/O note in README.md.
%
%   nTrial
%   The source pre-allocates nTrial = 10 and fills 1:sum(maskStim). Ten is
%   exactly the observed maximum in both sessions, but the minimum is 7, so
%   cells differ in how many repeats their average rests on. That is returned
%   rather than assumed.
%
%   See also SUBTHRESHOLDBYINTERPOLATION, TUTORIAL2A_STIMULUSPROTOCOL.

%% ---------------------------------------------------------------- options
opt = struct('Subthreshold', 'interpolate', 'Cells', [], 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
subMode = validatestring(opt.Subthreshold, {'interpolate', 'deposited', 'none'});
if isempty(opt.Cells)
    cells = (1:data.nCells)';
else
    cells = reshape(opt.Cells, [], 1);
end
assert(~isempty(data.stimChunk), 'Tutorial2:NoStimMap', ...
    'data.stimChunk is empty; load with ScanBlue true.');

nCells = data.nCells;
nFrame = data.nFrame;
timer = tic;

%% ------------------------------------------------------------ accumulate
spikeSum = zeros(nCells, nFrame);
subSum = zeros(nCells, nFrame);
nTrial = zeros(nCells, 1);
spk_t_blue = cell(nCells, 1);
blueOn = data.blue(:)' > 0;

% Spikes arrive as session-global frame indices; bucket them by chunk once so
% the pass does not rescan every cell's spike list 300 times.
chunkOfSpike = cell(nCells, 1);
for ii = cells'
    chunkOfSpike{ii} = floor((data.spk_t{ii} - 1) / nFrame) + 1;
    spk_t_blue{ii} = zeros(1, 0);
end

wantSub = ~strcmp(subMode, 'none');
for k = 1:data.nChunk
    active = cells(data.stimChunk(cells, k));
    if isempty(active); continue; end

    base = (k-1) * nFrame;

    % --- spikes -------------------------------------------------------
    for ii = active'
        local = data.spk_t{ii}(chunkOfSpike{ii} == k) - base;
        % Blue-evoked spikes: the source takes these BEFORE stripping t_noise,
        % so no chunk-boundary mask is applied here. It cannot matter either
        % way; blue is on over frames 394:7480 and t_noise is 1:140 and
        % 7861:8000. loadSparseOptoForEphys documents the asymmetry.
        isBlue = local >= 1 & local <= nFrame;
        isBlue(isBlue) = blueOn(local(isBlue));
        spk_t_blue{ii} = [spk_t_blue{ii}, local(isBlue) + base];
        spikeSum(ii, local(local >= 1 & local <= nFrame)) = ...
            spikeSum(ii, local(local >= 1 & local <= nFrame)) + 1;
    end

    % --- subthreshold -------------------------------------------------
    if wantSub
        switch subMode
            case 'deposited'
                block = data.readChunk('subthreshold', k);
                subSum(active, :) = subSum(active, :) + double(block(active, :));
            case 'interpolate'
                block = data.readChunk('voltage', k);
                for ii = active'
                    local = data.spk_t{ii}(chunkOfSpike{ii} == k) - base;
                    local = local(local >= 1 & local <= nFrame);
                    subSum(ii, :) = subSum(ii, :) + ...
                        double(subthresholdByInterpolation(block(ii, :), local));
                end
        end
    end

    nTrial(active) = nTrial(active) + 1;

    if opt.Verbose && mod(k, 50) == 0
        fprintf('  chunk %3d/%d  (%.0f s)\n', k, data.nChunk, toc(timer));
    end
end

%% -------------------------------------------------------------- finalise
% The source divides by sum(maskStim), the cell's own repeat count, not by a
% common nTrial. A cell with no stimulation would divide by zero; there are
% none in either session, but the guard keeps a subset run honest.
safe = max(nTrial, 1);
stacks = struct();
stacks.fr_mat_stim_avg = spikeSum ./ safe;
stacks.fr_mat_stim_avg(nTrial == 0, :) = NaN;
if wantSub
    stacks.traces_stim_avg = subSum ./ safe;
    stacks.traces_stim_avg(nTrial == 0, :) = NaN;
else
    stacks.traces_stim_avg = [];
end
stacks.nTrial = nTrial;
stacks.spk_t_blue = spk_t_blue;
% Carried through so later stages can enumerate a cell's stimulation repeats
% without re-reading blue. A repeat in which the cell happened not to fire is
% still a repeat, and rheobase counts it.
stacks.stimChunk = data.stimChunk;
stacks.nFrame = nFrame;
stacks.subthresholdMode = subMode;
stacks.cells = cells;

if opt.Verbose
    nBlue = cellfun(@numel, spk_t_blue(cells));
    fprintf('Stacked %d cells in %.0f s (%s subthreshold)\n', ...
        numel(cells), toc(timer), subMode);
    fprintf('  repeats/cell: median %g [%d, %d]\n', ...
        median(nTrial(cells)), min(nTrial(cells)), max(nTrial(cells)));
    fprintf('  blue-evoked spikes/cell: median %g [%d, %d]\n', ...
        median(nBlue), min(nBlue), max(nBlue));
end
end
