% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% --------------------------------------------------------------------------
function staMat = get_sta_mat_raw(spikeTimes,tracesAll,timeWindow,include_sync_spk)
% function staMat = get_sta_mat_2(spikeTimes,tracesAll,timeWindow,include_sync_spk)
% makes spike triggered average matrix: nCells x nCells x timeWindow
% spikeTimes: nCells x 1 cell array
% tracesAll: Times x nCells
% timeWindow: [tBefore tAfter] (both positive numbers)
% staMat: column triggered by row

% [T,nCells] = size(tracesAll);
[T,nRCells] = size(tracesAll);
nTCells = length(spikeTimes);

if isempty(include_sync_spk)
    include_sync_spk = 1;
end;

nBack = timeWindow(1);
nFront = timeWindow(2);
% staMat = zeros(nTCells, nRCells, nBack+nFront+1);

if include_sync_spk == 0;  % remove synchronous spikes
    spikeTMat = zeros(nRCells, T);
    for j = 1:nCells;
        spikeTMat(j,spikeTimes{j}) = 1;

%         spikeT_j = max(min(spikeTimes{j}+(-timeWindow(1):timeWindow(2))',T),1); % test for wider window
%         spikeTMat(j,spikeT_j) = 1; % test for wider window
    end;
end;


for j = 1:nTCells
%     if isempty(spikeTimes{j});continue;end
    spikeT = spikeTimes{j};
    spikeT(spikeT < nBack + 1) = [];
    spikeT(spikeT > T - nFront) = [];
    if isempty(spikeTimes{j});continue;end

    nSpike = length(spikeT);
    tMat = spikeT' + [-nBack:nFront];  % rectangular matrix with time window around each spike
    staIdx = (1:T:T*nRCells)' - 1 + shiftdim(tMat,-1);  % nRCells x nSpikes x length of time window
    traceST = tracesAll(staIdx);  % same shape as staIdx
    if include_sync_spk == 0;  % remove synchronous spikes
        synchMask = spikeTMat(:,spikeT);
        synchMask(synchMask == 1) = NaN;
        synchMask(synchMask == 0) = 1;
        traceST = traceST.*synchMask;
    end;
    staMat= traceST;
%     j
end