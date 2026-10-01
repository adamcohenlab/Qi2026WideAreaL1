function trig = buildOrthogonalTriggers(data, varargin)
%BUILDORTHOGONALTRIGGERS Hadamard-orthogonal spike triggers, per cell and per pair.
%
%   trig = BUILDORTHOGONALTRIGGERS(data) reproduces the trigger construction
%   at the top of `sparsePulseHad_spikelet_control.m` from an NWB load. It
%   returns the per-cell quantities; per-pair trigger sets come from
%   trig.pairTriggers(pre, post), which is a few microseconds once this has
%   run.
%
%   Options:
%     'GuardFrames'   140 (default). Half-width of the guard around the post
%                     cell's own Hadamard step edges and around chunk
%                     boundaries, and the amount each Hadamard step is shrunk
%                     at both ends.
%     'Verbose'       true (default).
%
%   FIELDS
%     tHad            nCells x 1 cell of frames. Cell i's spikes that fell
%                     inside a shrunk Hadamard step of i's own stimulation.
%                     Per-cell only; the pair guards are not applied yet.
%     blueHadOn       nCells x T logical. Hadamard step on, per cell.
%     blockedMask     nCells x T logical. Frames guarded out for that cell
%                     as a POST cell: within GuardFrames of its own step
%                     edges, or of a chunk boundary.
%     t_blueHadEdge   nCells x 1 cell of step-edge frames.
%     hadStepLevel    the blue level that marks a Hadamard step.
%     pairTriggers    @(pre,post) -> row vector of trigger frames.
%     poolTriggers    @(cells,post) -> pooled, guarded trigger frames of a
%                     set of cells under pair (cells, post) orthogonality.
%
%   WHAT MAKES A TRIGGER, AND WHY
%
%   The Hadamard block is every fourth recording chunk, and within it each
%   cell's blue stimulation follows one row of a Hadamard pattern, so at any
%   moment about half the cells are illuminated. That is what the method
%   exploits. For a pair (i, j) the usable spikes of i are those that fell
%   while **i was illuminated and j was not**:
%
%       orthoMask = blueHadOn(i,:) & ~blueHadOn(j,:)
%
%   Cell i is driven to spike, cell j sits at rest, and any deflection in j
%   locked to i's spike cannot be a shared response to shared light. That is
%   the whole design: it removes the common-illumination confound that the
%   parallel-stimulation case (`crossHadPara.mat`) carries, at the cost of
%   using only about half of i's evoked spikes.
%
%   Three exclusions apply on top:
%
%   1. Each Hadamard step is shrunk by GuardFrames at both ends before any
%      spike is counted. `cumsum` over rising and falling edges shifted by
%      +/-GuardFrames is an efficient way to write that. It keeps out the
%      onset transient of i's own stimulation.
%   2. Frames within GuardFrames of **j's** step edges are dropped. j's own
%      light switching on or off produces a large artifact in j's trace, and
%      it is not independent of i's spiking, because the Hadamard patterns
%      switch together.
%   3. Frames within GuardFrames of a chunk boundary are dropped, for every
%      pair, because the movie restarts there.
%
%   Exclusion 2 is what makes the trigger set a property of the *pair*
%   rather than of the presynaptic cell, and hence why the deposit stores it
%   per pair rather than folding it into a /units column.
%
%   A useful side effect of exclusion 3: no trigger is left within
%   GuardFrames of the start or end of the recording, so a spike-triggered
%   average is over exactly the same events whether the window is +/-20 or
%   +/-140 frames. The +/-20 tensor computed here and the +/-140 tensor in the
%   deposit are therefore directly comparable on their common samples.
%
%   IMPLEMENTATION NOTE
%   The source writes the per-pair set as
%   `setdiff(find(spk_t_mat_had(i,:) & orthoMask), [edges; breaks])`, which
%   scans all T frames per pair. Over the 9,124 near pairs that is 3e10
%   operations. Here the per-cell event list is found once and the pair
%   guards are applied by indexing into the mask rows at those events, which
%   is the same set by construction and about four orders of magnitude
%   faster. TUTORIAL4_VALIDATION.md checks the two against each other.
%
%   See also LOADSPARSEPULSEHADFORSPIKELET, TUTORIAL4A_ORTHOGONALTRIGGERS.

opt = struct('GuardFrames', 140, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
assert(~isempty(data.spk_t), 'Tutorial4:NoSpikes', ...
    'Load with ''LoadSpikes'', true before building triggers.');

nCells = data.nCells;
T      = data.T;
nFrame = data.nFrame;
g      = opt.GuardFrames;
timer  = tic;

%% ------------------------------------------- the Hadamard epoch, and levels
% Every fourth chunk, 1-based, exactly as stimTypeIdx{2} in the source.
hadChunkFrames = (1:nFrame) + (3*nFrame:nFrame*4:T)';
stimTypeMask = false(1, T);
stimTypeMask(hadChunkFrames(:)) = true;

% The Hadamard steps sit at the lowest nonzero blue level; the sparse pulses
% are brighter. Read the level off the data rather than hard-coding it.
blueVs = sort(unique(data.blue_all(1,:)));
blueVs = blueVs(2:end);
assert(~isempty(blueVs), 'Tutorial4:NoBlueLevels', 'Cell 1 has no nonzero blue level.');
hadStepLevel = blueVs(1);
blueHadOn = abs(data.blue_all - hadStepLevel) <= eps;

%% ------------------------------- per-cell Hadamard spikes, and the guards
% Chunk boundaries, shared by every cell.
tBreak = max(min((0:nFrame:T)' + (-g:g), T), 1);
breakMask = false(1, T);
breakMask(tBreak(:)) = true;

tHad          = cell(nCells,1);
t_blueHadEdge = cell(nCells,1);
blockedMask   = false(nCells, T);

for ii = 1:nCells
    % --- spikes inside a shrunk step of this cell's own stimulation
    blue_i = data.blue_all(ii,:);
    dRE = diff(blue_i); dRE = [0 dRE .* (dRE > 0)]; dRE = circshift(dRE,  g);
    dFE = diff(blue_i); dFE = [0 dFE .* (dFE < 0)]; dFE = circshift(dFE, -g);
    blue_i = cumsum(dRE + dFE);

    s = data.spk_t{ii};
    tHad{ii} = s(blue_i(s) > 0 & stimTypeMask(s));

    % --- this cell's own step edges, and the guard band around them
    e = find([false, abs(diff(blueHadOn(ii,:))) > 0]);
    t_blueHadEdge{ii} = e;
    blockedMask(ii,:) = intervalMask(e, g, T) | breakMask;
end

%% ------------------------------------------------------------- assemble
trig = struct();
trig.tHad          = tHad;
trig.blueHadOn     = blueHadOn;
trig.blockedMask   = blockedMask;
trig.t_blueHadEdge = t_blueHadEdge;
trig.breakMask     = breakMask;
trig.stimTypeMask  = stimTypeMask;
trig.hadStepLevel  = hadStepLevel;
trig.guardFrames   = g;
trig.nCells        = nCells;
trig.T             = T;
trig.pairTriggers  = @(pre, post) pairTriggers(trig, pre, post);
trig.poolTriggers  = @(cells, pre, post) poolTriggers(trig, cells, pre, post);

if opt.Verbose
    n = cellfun(@numel, tHad);
    fprintf('Hadamard triggers: step level %.3g, %d of %d frames in the Hadamard epoch\n', ...
        hadStepLevel, sum(stimTypeMask), T);
    fprintf('  per-cell Hadamard spikes before pair guards: median %d [%d, %d]\n', ...
        median(n), min(n), max(n));
    fprintf('  %.1f s\n', toc(timer));
end
end

%% ------------------------------------------------------------- helpers
function m = intervalMask(centres, g, T)
%INTERVALMASK Logical mask of every frame within g of one of `centres`.
%   Built with a difference array, which is O(T) rather than O(T * numel(g)).
m = false(1, T);
if isempty(centres)
    return
end
lo = max(centres - g, 1);
hi = min(centres + g, T);
d = zeros(1, T+1);
d(lo) = d(lo) + 1;
d(hi+1) = d(hi+1) - 1;
m = cumsum(d(1:T)) > 0;
end

function t = pairTriggers(trig, pre, post)
%PAIRTRIGGERS Cell `pre` spiking while `pre` was lit and `post` was not.
%   blueHadOn(pre,:) is tested even though tHad already lives inside a
%   shrunk step of pre's own stimulation, because the source writes it that
%   way and the two differ if a Hadamard chunk ever carries a blue level
%   other than the step level.
t = trig.tHad{pre};
t = t(trig.blueHadOn(pre, t) & ~trig.blueHadOn(post, t) & ~trig.blockedMask(post, t));
t = reshape(t, 1, []);
end

function t = poolTriggers(trig, cells, pre, post)
%POOLTRIGGERS Pooled guarded spikes of `cells` under pair (pre, post) optics.
%   The orthogonality mask belongs to the pair being tested, not to the cells
%   supplying the spikes: that is what keeps the post cell in the same
%   optical state in the null as in the observation.
t = unique([trig.tHad{cells}]);
t = t(trig.blueHadOn(pre, t) & ~trig.blueHadOn(post, t) & ~trig.blockedMask(post, t));
t = reshape(t, 1, []);
end
