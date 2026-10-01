function staRow = spikeletStaRow(pre, trig, traceMat, varargin)
%SPIKELETSTAROW One presynaptic cell's spikelet waveform against every cell.
%
%   staRow = SPIKELETSTAROW(pre, trig, traceMat) returns nCells-by-nTau: for
%   each post cell, the average of that cell's spike-subtracted trace on the
%   orthogonal Hadamard triggers of the pair.
%
%   traceMat is T-by-nCells, the transpose of the usual cell-by-time layout,
%   because that is what get_sta_mat_single wants. Keep it `single`: as
%   `double` it is 8 GB rather than 4 GB, and each column is promoted here
%   anyway, which costs 25 MB per post cell instead of doubling the
%   resident set.
%
%   Options:
%     'nBack'   20 (default)
%     'nFront'  20 (default)
%     'Verbose' false
%
%   WHY A WHOLE ROW AND NOT A PAIR
%   The correction in stage 4c is estimated from the presynaptic cell's
%   relationship to *every* other cell: the distance trend needs cells at
%   all distances, and the common mode is defined as the average over cells
%   beyond 800 um. Neither can be computed pair by pair. The row is
%   therefore the natural unit of work on this path, and it is also why the
%   deposited tensor, which stops at 400 um, cannot be used as the input:
%   the far field it would need is stored as zeros.
%
%   The self entry is left at zero. The orthogonality mask for a cell
%   against itself is empty by construction, which is exactly what the
%   source produces, and stage 4c excludes it from both fits anyway.
%
%   About 3 to 4 s for 320 post cells at 3.17 M frames.
%
%   See also BUILDORTHOGONALTRIGGERS, CORRECTCOMMONMODE.

opt = struct('nBack', 20, 'nFront', 20, 'Verbose', false);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end

nCells = size(traceMat, 2);
staRow = zeros(nCells, opt.nBack + opt.nFront + 1);
timer = tic;

for jj = 1:nCells
    if jj == pre
        continue                    % empty trigger set, stays zero
    end
    staRow(jj,:) = get_sta_mat_single({trig.pairTriggers(pre, jj)}, ...
        double(traceMat(:,jj)), [opt.nBack opt.nFront], 1);
end

if opt.Verbose
    fprintf('  STA row for pre cell %d: %d post cells in %.1f s\n', pre, nCells, toc(timer));
end
end
