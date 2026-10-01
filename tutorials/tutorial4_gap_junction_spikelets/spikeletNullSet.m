function out = spikeletNullSet(pre, post, trig, tracesSub, distMat, varargin)
%SPIKELETNULLSET Far-cell resampled null amplitudes for one directed pair.
%
%   out = SPIKELETNULLSET(pre, post, trig, tracesSub, distMat) reproduces the
%   control block of `sparsePulseHad_spikelet_control.m` for a single pair:
%   it draws nRand surrogate trigger sets from cells too far from the post
%   cell to be coupled to it, averages the post cell's trace on each, and
%   returns the resulting null distribution of spikelet amplitudes.
%
%   Options (defaults are the source literals for the published run):
%     'nRand'      5000
%     'rLim'       400    um; "far" means further than this FROM THE POST CELL
%     'nBack_c'    10     null window half-width, both sides
%     'dx'         6.5
%     'Seed'       []     rng(Seed,'twister') before drawing, [] leaves rng alone
%
%   FIELDS
%     nullAmp      nRand x 1 null spikelet amplitudes
%     triggers     the observed trigger frames for this pair
%     nSpk         numel(triggers)
%     nCtrl        size of the pool the surrogates were drawn from
%
%   WHAT THE NULL HOLDS FIXED, AND WHAT IT BREAKS
%
%   The question is whether the post cell's average deflection at offset
%   zero is bigger than it would be if the trigger spikes had come from a
%   cell that cannot be coupled to it. So the surrogate keeps:
%
%     - the post cell. Every null average is of the SAME trace, so the post
%       cell's own noise, spike rate and spike-subtraction residue are
%       matched exactly rather than modelled.
%     - the number of triggers. nSpk is read off the observed pair, so the
%       null and the observation average the same number of snippets and
%       have the same standard error. Without this, pairs with many spikes
%       would look significant for free.
%     - the stimulation condition. Surrogate spikes are drawn only from
%       frames satisfying the same orthogonality mask, blueHadOn(pre) and
%       not blueHadOn(post), and the same edge and chunk-boundary guards.
%       The post cell is therefore in the same optical state in the null as
%       in the observation.
%
%   and breaks exactly one thing: which cell fired. That is the whole
%   inference. Anything shared between a near pre cell and a far one —
%   network synchrony, illumination, the post cell's own statistics —
%   appears in both and cancels.
%
%   "Far" is measured from the POST cell, not the pre cell. The claim being
%   tested is that this post cell is coupled to this pre cell, so the
%   surrogate has to be a cell the post cell could not be coupled to.
%
%   The pooled draw is uniform over the pooled event list, so a far cell
%   that fires more contributes more surrogate triggers. That is deliberate:
%   the null should reflect the marginal timing statistics of Hadamard-epoch
%   spiking, not an equal vote per cell.
%
%   FAITHFUL ODDITY
%   The source draws indices with `round(rand*(nCtrl-1)+1)`, which gives the
%   first and last element of the pool half the weight of the others. The
%   effect is O(1/nCtrl) on a pool of tens of thousands. It is reproduced
%   here rather than silently fixed, so that the tutorial's p-values match
%   the published ones.
%
%   See also SPIKELETAMPLITUDE, TUTORIAL4D_SPIKELETDETECTION.

opt = struct('nRand', 5000, 'rLim', 400, 'nBack_c', 10, 'dx', 6.5, 'Seed', []);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
if ~isempty(opt.Seed)
    rng(opt.Seed, 'twister');
end

%% ------------------------------------------------- the observed triggers
triggers = trig.pairTriggers(pre, post);
nSpk = numel(triggers);

out = struct('nullAmp', zeros(opt.nRand,1), 'triggers', reshape(triggers,1,[]), ...
    'nSpk', nSpk, 'nCtrl', 0);
if nSpk == 0
    out.nullAmp = nan(opt.nRand,1);
    return
end

%% ------------------------------------------------------ the surrogate pool
idxFar = find(distMat(post,:) > opt.rLim / opt.dx);
assert(~isempty(idxFar), 'Tutorial4:NoFarCells', ...
    'No cell is further than %g um from post cell %d.', opt.rLim, post);
spk_ctrl = trig.poolTriggers(idxFar, pre, post);
nCtrl = numel(spk_ctrl);
assert(nCtrl > 1, 'Tutorial4:EmptyNullPool', ...
    'Surrogate pool for pair (%d,%d) holds %d events.', pre, post, nCtrl);
out.nCtrl = nCtrl;

%% ---------------------------------------------------------- draw and average
idxSampl = round(rand(nSpk, opt.nRand) * (nCtrl - 1) + 1);
spk_sampl = mat2cell(spk_ctrl(idxSampl)', ones(1, opt.nRand), nSpk);

nullSta = get_sta_mat_single(spk_sampl, double(tracesSub(post,:))', ...
    [opt.nBack_c opt.nBack_c], 1);          % nRand x 1 x nTau
out.nullAmp = spikeletAmplitude(squeeze(nullSta), opt.nBack_c);
out.nullAmp = out.nullAmp(:);
end
