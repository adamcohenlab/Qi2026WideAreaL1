% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% Modified from the HALS routines of CaImAn-MATLAB
% (https://github.com/flatironinstitute/CaImAn-MATLAB), and so
% distributed under the GNU General Public License v2; see license.txt
% in this folder.
% --------------------------------------------------------------------------
function [A,C,b,f] = ALS_nonneg_BackgroundUpdateOnly(Y, A, b, maxIter,cellMovIdx,rowCol,doPlot)
%% Alternating least square to factorize voltage imaging movie

% doPlot draws the per-iteration diagnostic figures. Default off, for the same
% reason as in ALS_nonneg: a full-frame rebuild plus drawnow on every iteration.
% Numerics are unaffected.
if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = false;
end

%% YQ parameters

nCells = size(A,2);
nBg = size(b,2);

nR = rowCol(1); nC = rowCol(2);

C = zeros(nCells,size(Y,2));
f = zeros(nBg,size(Y,2));
%% find A and b on "motion-free" movie, Ymc

 
for miter=1:maxIter
%     spatial = SeeResiduals_vec(Ymc0,[C_hi;f_hi;m],0,1);
%     d = spatial(:,nCells+nBg+1:end);
%     Ymc = Ymc0 - d*m;
    
    temporal = HALS_temporal(max(Y,0), [A b], zeros(nCells+nBg,size(Y,2)), 5,[],true);
    C_pos = temporal(1:nCells,:);
    f_pos = temporal(nCells+1:end,:);
    
    temporal = HALS_temporal(max(-Y,0), [A b], zeros(nCells+nBg,size(Y,2)), 5,[],true);
    C_neg = -temporal(1:nCells,:);
    f_neg = -temporal(nCells+1:end,:);
    
    C = C_pos+C_neg;
    f = f_pos+f_neg;
    
    spatial = HALS_spatial(Y, [A b], [C; f], [], 5);
%     A = spatial(:,1);
%     b = spatial(:,2:end);
    b(:,end) = spatial(:,end);
    
    if doPlot
        tmp = zeros(nR*nC,size(A,2) + size(b,2));
        tmp(cellMovIdx,:) = [A b];
        figure(1);clf;moviefixsc(toimg(tmp./std(tmp,[],1),nR,nC))
        figure(2);clf;plot(zscore([C;f]'))
        drawnow
    end
end
%     temporal = HALS_temporal(Ymc, [A b], zeros(nCells+nBg,size(Y,2)), 5,[],false);
%     C = temporal(1:nCells,:);
%     f = temporal(nCells+1:end,:);
end