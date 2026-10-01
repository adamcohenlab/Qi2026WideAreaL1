% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: the IPSP detection and fitting code written for this study
% (folder 'psp fitting' in the author's analysis tree).
% --------------------------------------------------------------------------
function Ynull = generatePairNullSet(preCell, postCell, tBlue, tracesAll, distMat, tBack, tFront, rCutNear, rCutFar, nNull, removedTimes, stream)
%GENERATEPAIRNULLSET Construct far-stimulation null averages for one pair.
%   This is the retaining-free version of the loop in
%   crosstalkTemplate_nullSet.m: it creates the random far-stimulation
%   averages only for the current pair, so no giant all-pairs null tensor is
%   written to disk.  It requires get_sta_mat_single on the MATLAB path.
%   stream is an optional RandStream used only for selecting null events.

if nargin < 12
    stream = [];
end

validateattributes(preCell, {'numeric'}, {'scalar','integer','positive'});
validateattributes(postCell, {'numeric'}, {'scalar','integer','positive'});
validateattributes(tBlue, {'cell'}, {'vector','nonempty'});
validateattributes(tracesAll, {'numeric'}, {'2d','real','finite'});
validateattributes(distMat, {'numeric'}, {'2d','real','finite','square'});
validateattributes(tBack, {'numeric'}, {'scalar','integer','nonnegative'});
validateattributes(tFront, {'numeric'}, {'scalar','integer','nonnegative'});
validateattributes(nNull, {'numeric'}, {'scalar','integer','positive'});
if exist('get_sta_mat_single', 'file') ~= 2
    error('generatePairNullSet:MissingDependency', 'get_sta_mat_single must be on the MATLAB path.');
end
if preCell > numel(tBlue) || postCell > size(tracesAll,1)
    error('generatePairNullSet:CellIndex', 'Cell indices exceed tBlue or tracesAll.');
end

farCells = find(distMat(postCell,:) > rCutFar);
nearCells = find(distMat(postCell,:) <= rCutNear);
excludedTimes = unique([tBlue{nearCells}]);
farTimes = cell(numel(farCells), 1);
for k = 1:numel(farCells)
    farTimes{k} = setdiff(tBlue{farCells(k)}, excludedTimes);
end
farTimes = farTimes(~cellfun(@isempty, farTimes));
if isempty(farTimes)
    error('generatePairNullSet:NoFarEvents', 'No usable far-stimulation events are available for post cell %d.', postCell);
end
farTimes = setdiff(unique([farTimes{:}]), removedTimes(:).');
nStim = numel(tBlue{preCell});
if nStim == 0 || isempty(farTimes)
    error('generatePairNullSet:InsufficientEvents', 'The pair has no stimulation events or no remaining far events.');
end

if isempty(stream)
    sampleIndex = randi(numel(farTimes), nStim, nNull);
else
    sampleIndex = randi(stream, numel(farTimes), nStim, nNull);
end
eventCells = mat2cell(farTimes(sampleIndex).', ones(1, nNull), nStim);
raw = squeeze(get_sta_mat_single(eventCells, tracesAll(postCell,:).', [tBack tFront], 1));
nTime = tBack + tFront + 1;
if isvector(raw)
    raw = raw(:);
end
if size(raw, 1) == nTime
    Ynull = raw;
elseif size(raw, 2) == nTime
    Ynull = raw.';
else
    error('generatePairNullSet:UnexpectedStaShape', 'get_sta_mat_single returned [%d x %d], expected a dimension of %d.', size(raw,1), size(raw,2), nTime);
end
if size(Ynull,2) ~= nNull
    error('generatePairNullSet:UnexpectedNullCount', 'Expected %d null waveforms, received %d.', nNull, size(Ynull,2));
end
Ynull = Ynull - mean(Ynull, 1);
end
