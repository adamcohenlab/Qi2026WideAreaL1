% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Authors: Yitong Qi and AEC (as credited in the original header).
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% --------------------------------------------------------------------------
function template = fitSelfBlueSta(staMat, tPrev, tOn, tOff);
% function template = fitSelfBlueSta(staMat, tPrev, tOn, tOff);
% Fit a template to the self blue-stim-triggered STA for a cell;
% tPrev = time of the previous DMD switch;
% ton = time of blue light on;
% toff = time of blue light off;
% YQ and AEC 12/20/2022

trace = median(staMat);
template = zeros(size(trace));
if tPrev > 1;
    template(1:tPrev-1) = mean(trace(1:tPrev-1));
end;
template(tPrev:tOn-1) = mean(trace((tOn - 22):(tOn-2)));
tau = tOn:(tOff-1);
tau = tau - mean(tau);
P = polyfit(tau, trace(tOn:(tOff-1)), 3);
template(tOn:(tOff-1)) = polyval(P, tau);
template(tOff:end) = smooth(trace(tOff:end), 5);

