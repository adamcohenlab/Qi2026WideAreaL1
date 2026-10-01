function [baseline, nEvents] = farPulseBaseline(postCells, data, tBack, tFront, rCutNear, rCutFar, verbose)
%FARPULSEBASELINE Average response of a cell to distant cells' blue pulses.
%
%   [baseline, nEvents] = FARPULSEBASELINE(postCells, data, tBack, tFront, ...
%                                          rCutNear, rCutFar, verbose)
%
%   For each cell in postCells, averages that cell's own trace around blue
%   pulses delivered to cells FARTHER than rCutFar away, excluding any time at
%   which a cell within rCutNear was also pulsed, and excluding chunk-start
%   frames.
%
%   What this is for. A pulse aimed at cell i does not only stimulate cell i.
%   Scattered blue light and direct optical crosstalk reach the whole field, so
%   every cell's trace carries a stimulus-locked artifact on every pulse,
%   whether or not it is connected to the stimulated cell. Averaging over pulses
%   delivered far enough away that no synaptic effect is plausible isolates that
%   artifact. Subtracting it is the first of the two baseline corrections in
%   tutorial3b.
%
%   The result is indexed by the RECORDED (postsynaptic) cell, not the
%   stimulated one, which is why the subtraction in tutorial3b broadcasts along
%   the presynaptic axis.
%
%   Source: the Monte Carlo baseline loop in synapticConn_baselineGen.m.
%
%   NOTE ON AN EQUIVALENT SIMPLIFICATION. The source excludes near-cell times
%   from each far cell's list separately and then takes the union. Set
%   difference distributes over union, so taking the union first and
%   subtracting once is identical and much faster. generatePairNullSet.m
%   already does it this way.

if nargin < 7, verbose = false; end

nCells = data.nCells;
nTau = tBack + tFront + 1;
postCells = postCells(:)';
baseline = zeros(numel(postCells), nTau, 'single');
nEvents = zeros(numel(postCells), 1);

timer = tic;
for k = 1:numel(postCells)
    ii = postCells(k);

    idxCellFar  = find(data.distMat(ii,:) > rCutFar);
    idxCellNear = find(data.distMat(ii,:) <= rCutNear);

    farTimes = unique([data.t_blue{idxCellFar}]);
    farTimes = setdiff(farTimes, unique([data.t_blue{idxCellNear}]));
    farTimes = setdiff(farTimes, data.removedTimes(:)');

    if isempty(farTimes)
        error('Tutorial3:NoFarEvents', 'Cell %d has no usable far-pulse events.', ii);
    end

    baseline(k,:) = squeeze(get_sta_mat_single({farTimes}, ...
        data.traces_all(ii,:)', [tBack tFront], 1))';
    nEvents(k) = numel(farTimes);

    if verbose && (mod(k, 25) == 0 || k == numel(postCells))
        fprintf('    [%6.1f s] far-pulse baseline %d/%d\n', ...
            toc(timer), k, numel(postCells));
    end
end
end
