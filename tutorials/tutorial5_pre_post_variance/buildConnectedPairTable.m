function [tbl, dRef] = buildConnectedPairTable(sessions, labels, responseField)
%BUILDCONNECTEDPAIRTABLE One row per connected pair, pooled over sessions.
%
%   [tbl, dRef] = BUILDCONNECTEDPAIRTABLE(sessions, labels, responseField)
%   turns each session's cell-by-cell matrices (from loadSparsePulseHadForLme)
%   into a long table with one row per connected, directed pair whose
%   response is positive. responseField is 'ipspAmpTestN2' or 'ipspDecay'.
%
%   COLUMNS
%     animal            categorical session label (fixed effect)
%     pre, post         categorical cell identity, e.g. 'M27_unit007'. The
%                       session label is part of the name, so cell 7 of one
%                       mouse and cell 7 of the other are different levels.
%     preUnitId, postUnitId   zero-based NWB unit ids
%     distanceUm        pair distance in micrometres
%     logDistanceRatio  log(distanceUm / dRef)
%     distanceCentered  distanceUm - dRef
%     response          the raw response (amplitude or decay constant)
%     logResponse       log(response), the modelled variable
%     seTau             decay standard error for the pair (NaN if absent)
%
%   dRef is the median positive distance over ALL measured pairs of all
%   sessions (connected or not, diagonal excluded), in micrometres. It is
%   only a centring point for distance: changing it moves the intercept and
%   nothing else. It is about 1.15 mm, well beyond the 60-400 um span of the
%   connected pairs, because most measured pairs are far apart; it is kept
%   so the numbers match the published analysis exactly.
%
%   ROW, LEVEL AND COLUMN ORDER
%   Rows follow MATLAB column-major order within each session, sessions in
%   the order given, and categorical levels are listed in that same order.
%   Neither affects the fitted model, but both fix which random number goes
%   to which pair and which cell in the parametric bootstrap. animal is the
%   last column, as in the source table, because fitlme builds its
%   fixed-effect design in table-column order. Keeping all three identical
%   to the source analysis makes its bootstrap reproducible to the last
%   digit.

nSessions = numel(sessions);
assert(numel(labels) == nSessions, 'Tutorial5:Labels', ...
    'Give one label per session.');

% Reference distance, pooled over sessions.
allDistances = cell(nSessions, 1);
for s = 1:nSessions
    [measured, distanceUm] = measuredPairs(sessions{s}, responseField);
    allDistances{s} = distanceUm(measured & distanceUm > 0);
end
dRef = median(vertcat(allDistances{:}));

parts = cell(nSessions, 1);
preLevels = cell(nSessions, 1);
postLevels = cell(nSessions, 1);
for s = 1:nSessions
    session = sessions{s};
    [measured, distanceUm] = measuredPairs(session, responseField);
    response = session.(responseField);
    connected = measured & session.pMat & response > 0;

    [preRow, postColumn] = find(connected);          % column-major order
    pairIndex = sub2ind(size(response), preRow, postColumn);
    preUnitId = session.unitId(preRow);
    postUnitId = session.unitId(postColumn);
    preName = arrayfun(@(u) sprintf('%s_unit%03d', labels{s}, u), preUnitId, ...
        'UniformOutput', false);
    postName = arrayfun(@(u) sprintf('%s_unit%03d', labels{s}, u), postUnitId, ...
        'UniformOutput', false);
    preLevels{s} = unique(preName);
    postLevels{s} = unique(postName);

    d = distanceUm(pairIndex);
    parts{s} = table(repmat(labels(s), numel(pairIndex), 1), preName, postName, ...
        preUnitId, postUnitId, d, log(d ./ dRef), d - dRef, ...
        response(pairIndex), log(response(pairIndex)), session.seTau(pairIndex), ...
        'VariableNames', {'animal','pre','post','preUnitId','postUnitId', ...
        'distanceUm','logDistanceRatio','distanceCentered','response', ...
        'logResponse','seTau'});
end

tbl = vertcat(parts{:});
% fitlme orders fixed-effect columns by their position in the table. The
% source table carries animal after the distance terms, so animal moves to
% the last column here too. The model is the same either way, but a different
% column order changes floating-point rounding at the 1e-10 level.
tbl = movevars(tbl, 'animal', 'After', width(tbl));
tbl.animal = categorical(tbl.animal, labels(:));
tbl.pre = categorical(tbl.pre, vertcat(preLevels{:}));
tbl.post = categorical(tbl.post, vertcat(postLevels{:}));
end

function [measured, distanceUm] = measuredPairs(session, responseField)
% A pair is "measured" when its response and distance are finite. Self-pairs
% are never measured.
distanceUm = session.distMat * session.dx;
measured = isfinite(session.(responseField)) & isfinite(distanceUm);
measured(logical(eye(size(measured)))) = false;
end
