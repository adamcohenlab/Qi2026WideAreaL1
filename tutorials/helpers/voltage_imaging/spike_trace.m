% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% --------------------------------------------------------------------------
function [t_rec,template] = spike_trace(trace,spikes,halfWindow)
T = length(trace);
spikes = spikes(spikes - halfWindow >= 1 & spikes + halfWindow <= T);
spkMat = reshape(spikes,[],1) + (-halfWindow:halfWindow);
spkMat = trace(spkMat);
template = mean(spkMat,1);
spkTMat = zeros(size(trace));
spkTMat(spikes) = 1;
t_rec = conv(spkTMat,template,'same');
