function isolated = isolatedBlueSpikes(spk_t_blue, nFrLim)
%ISOLATEDBLUESPIKES Keep blue-evoked spikes with no close neighbour either side.
%
%   isolated = ISOLATEDBLUESPIKES(spk_t_blue, nFrLim) returns, for each cell, the
%   subset of blue-evoked spikes separated from the previous AND the next spike
%   by more than nFrLim frames. nFrLim defaults to 20 (25.4 ms at 1.27 ms/frame).
%
%   Reproduces spk_t_blue_sta in sparseOpto_ephysProp_meas_cleanup.m:
%
%       ISI     = diff(spkT);
%       ISI_rev = spkT(1:end-1) - spkT(2:end);
%       idxKeep = find([ISI 0] > nFrLim & [0 ISI_rev] < -nFrLim);
%
%   WHY ISOLATION MATTERS HERE
%   The afterdepolarisation is read from the spike-triggered average 7 frames
%   (8.9 ms) after the peak. At the firing rates the ramp drives, a following
%   spike lands inside that window often enough to dominate the average, and
%   what would be measured is the mean inter-spike interval rather than the
%   ADP. Requiring isolation on both sides removes that contamination; it also
%   removes the preceding spike's own ADP from the baseline window at -11
%   frames.
%
%   THE PADDING CONVENTION
%   The source pads ISI with a trailing 0 and ISI_rev with a leading 0, so the
%   first and last spike of each cell can never be kept: [ISI 0] > nFrLim fails
%   on the last, and [0 ISI_rev] < -nFrLim fails on the first. That is a
%   deliberate edge convention, not an off-by-one, and it is preserved. Two
%   spikes per cell are lost, against the tens of thousands that remain.
%
%   See also SELFTRIGGEREDAVERAGE, TUTORIAL2D_ADPANDMEMBRANECONSTANT.

narginchk(1, 2);
if nargin < 2 || isempty(nFrLim)
    nFrLim = 20;
end

nCells = numel(spk_t_blue);
isolated = cell(nCells, 1);
for ii = 1:nCells
    spkT = sort(reshape(spk_t_blue{ii}, 1, []));
    if numel(spkT) < 2
        isolated{ii} = zeros(1, 0);
        continue;
    end
    ISI = diff(spkT);
    ISI_rev = spkT(1:end-1) - spkT(2:end);
    % The source writes idxKeep = find(...) and then spkT(idxKeep); the mask
    % is used directly here, which selects the same spikes.
    keep = [ISI 0] > nFrLim & [0 ISI_rev] < -nFrLim;
    isolated{ii} = spkT(keep);
end
end
