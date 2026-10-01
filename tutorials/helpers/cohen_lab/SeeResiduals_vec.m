% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB library (imaging/motion_correction).
% Modified from SeeResiduals, Cohen Lab, Harvard University.
% --------------------------------------------------------------------------
function [res,C] = SeeResiduals_vec(vec, refSig, remOffset,outChoice)
% function out = SeeResiduals_vec(vec, refSig, remOffset,outChoice)
% vec: N X T (number of pixels x time)
% refSign: matrix orientation doesn't matter
% outChoice: 0, out = vec - C*refSig; 1, out = C;
    [~,nFrame] = size(vec);
    if ~exist('remOffset','var'), remOffset = 1;end
    if ~exist('outChoice','var'), outChoice = 0;end

    [refY, refX] = size(refSig);
    if refY == nFrame 
        refSig = refSig';
    elseif refX ~= nFrame
        'Size of refSig must be [N, # of frames]'
        return
    end
    
    if remOffset
%         dRef = bsxfun(@minus, refSig, mean(refSig, 2));  % remove the DC offset from refSig because we will account for offset separately.
        dRef = refSig - mean(refSig,2);
        I = [ones(1, nFrame); dRef];  % Add a row of ones to measure the DC offset at each pixel.
    else
        I = refSig;
    end

%     C = vec*I'*inv(I*I');  % Linear algebra to find the pseudo-inverse
    C = vec*I' / (I*I');  % Linear algebra to find the pseudo-inverse
%     C = (vec' \ I')'; % least square regression
    switch outChoice
        case 0
            residVec = vec - C*I;  % Look at the residuals
            res = residVec;
        case 1
            res = C;
        case 2
%             residVec = vec - C*I;  % Look at the residuals
%             res = vecnorm(residVec,2,2);
            res = vec - C*I;  % Look at the residuals
    end