function tracesSub = subtractSpikes(traces, spk_t, varargin)
%SUBTRACTSPIKES Replace each cell's own spikes with a linear interpolation.
%
%   tracesSub = SUBTRACTSPIKES(traces, spk_t) reproduces the spike-removal
%   block of `sparsePulseHad_spikelet_control.m`: for each cell, drop the
%   samples in a window around every one of that cell's own spikes and
%   linearly interpolate across the gap.
%
%   Options:
%     'SpikeWindow'  -3:5 (default), the source literal, in frames.
%     'Verbose'      true (default), prints progress every 25 cells.
%
%   WHY
%   The spikelet is a small deflection in the POST cell locked to the PRE
%   cell's spike. The post cell also fires on its own, and at 1.27 ms per
%   frame a post spike is tens of times larger than the spikelet. Post
%   spikes that happen to fall near a pre spike would dominate the average.
%   They are not removed because they are noise but because they are the
%   wrong signal: what survives here is the subthreshold trace, which is
%   where an electrical coupling shows up.
%
%   The window is asymmetric, -3:5, because the after-spike repolarisation
%   is slower than the rise.
%
%   NOTE ON THE DIAGONAL
%   Only the post cell's trace is spike-subtracted. The deposited tensor's
%   diagonal is a *self* average taken on the unsubtracted trace, so it is
%   about one spike height tall and is not comparable with anything computed
%   here. Stage 4c excludes it explicitly.
%
%   This is a full pass over an nCells-by-nFrames matrix and is the single
%   most expensive step in the tutorial: about 4 GB out and a few minutes for
%   320 cells at 3.17 M frames. Nothing downstream needs the unsubtracted
%   trace except the diagonal, so clearing `traces` afterwards is safe.
%
%   See also BUILDORTHOGONALTRIGGERS, TUTORIAL4B_SPIKELETWAVEFORM.

opt = struct('SpikeWindow', -3:5, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end

[nCells, T] = size(traces);
assert(numel(spk_t) == nCells, 'Tutorial4:SpikeRowMismatch', ...
    'traces has %d rows but spk_t has %d.', nCells, numel(spk_t));

tracesSub = zeros(size(traces), 'single');
w = opt.SpikeWindow(:);
timer = tic;

for ii = 1:nCells
    % setdiff on the spike-window union, exactly as the source writes it.
    t_i = setdiff(1:T, spk_t{ii} + w);
    tracesSub(ii,:) = interp1(t_i, double(traces(ii,t_i)), 1:T, 'linear');

    if opt.Verbose && (mod(ii,25) == 0 || ii == nCells)
        fprintf('  spike subtraction %3d/%d  (%.0f s)\n', ii, nCells, toc(timer));
    end
end

% interp1 leaves NaN outside the interpolation range, which happens only if a
% spike sits within the window of frame 1 or frame T. Carry the nearest good
% sample into those, so the array has no NaN for downstream averaging.
bad = isnan(tracesSub);
if any(bad(:))
    nBadCell = sum(any(bad,2));
    for ii = find(any(bad,2))'
        good = ~bad(ii,:);
        tracesSub(ii,~good) = interp1(find(good), double(tracesSub(ii,good)), ...
            find(~good), 'nearest', 'extrap');
    end
    if opt.Verbose
        fprintf('  filled %d edge NaN samples in %d cells\n', sum(bad(:)), nBadCell);
    end
end
end
