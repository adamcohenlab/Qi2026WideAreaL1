% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% Modified from the HALS routines of CaImAn-MATLAB
% (https://github.com/flatironinstitute/CaImAn-MATLAB), and so
% distributed under the GNU General Public License v2; see license.txt
% in this folder.
% --------------------------------------------------------------------------
function [A,C,b,f] = ALS_nonneg(Ymc, A, b, nFrames,cellMovIdx,rowCol,maxIter,resLim,resDiffLim,doPlot)
%% Alternating least square to factorize voltage imaging movie
if ~exist('resLim','var')
    resLim = 0;
end
if ~exist('resDiffLim','var')
    resDiffLim = 0;
end
% doPlot draws the per-iteration diagnostic figures. Default off: they rebuild a
% full-frame image and call drawnow on every iteration, which costs real time and
% leaves stray figures behind when this is called in a loop over cells. Pass true
% to restore the original interactive behaviour. Numerics are unaffected.
if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = false;
end

%% YQ parameters

nCells = size(A,2);
nBg = size(b,2);
nR = rowCol(1); 
nC = rowCol(2);

res = zeros(maxIter,nCells+nBg);
%% find A and b on "motion-free" movie, Ymc
for miter=1:maxIter

    tmp = zeros(nR*nC,size(A,2) + size(b,2));
    tmp(cellMovIdx,:) = [A b];
    im_lo = imgaussfilt(toimg(tmp,nR,nC),1);
%     im_lo = toimg(tmp,nR,nC);
    A_lo = im_lo(:,:,1);
    A_lo = single(A_lo(cellMovIdx));
    b_lo = tovec(im_lo(:,:,2:end));
    b_lo = [single(b_lo(cellMovIdx,:))];
    temporal = HALS_temporal(Ymc, [A_lo b_lo], zeros(nCells+nBg,nFrames), 5,[],true);
    
%     temporal = HALS_temporal(Ymc, [A b], zeros(nCells+nBg,nFrames), 5,[],true);
    C = zscore(temporal(1:nCells,:),[],2);
    C = C - min(C,[],2);
    f = zscore(temporal(nCells+1:end,:),[],2);
    f = f - min(f,[],2);
    
    
    spatial = HALS_spatial(Ymc, [A b], [C; f], [], 5);
%     spatial = HALS_spatial(Ymc, [A b], fft_clean([C; f]')', [], 5);
    res(miter,:) = sum(abs(spatial - [A b]),1);
    A = spatial(:,1:nCells);
    b = spatial(:,nCells+1:end);

    if doPlot
        tmp = zeros(nR*nC,size(A,2) + size(b,2));
        tmp(cellMovIdx,:) = [A b];
        figure(1);clf;moviefixsc(zscore(toimg(tmp,nR,nC),[],[1 2]))
        figure(2);clf;plot([C;f]')
        figure(3);clf;
        subplot(2,1,1)
        plot(res)
        xlim([2 maxIter])
        subplot(2,1,2)
        plot(abs(diff(res)))
        xlim([2 maxIter])
        drawnow
    end

    if miter >= 3
        if res(miter,1) < resLim && abs(res(miter,1) - res(miter-1,1)) < resDiffLim
            break
        end
    end
end

temporal = HALS_temporal(Ymc, [A b], zeros(nCells+nBg,nFrames), 5,[],true);
C = temporal(1:nCells,:);
f = temporal(nCells+1:end,:);
