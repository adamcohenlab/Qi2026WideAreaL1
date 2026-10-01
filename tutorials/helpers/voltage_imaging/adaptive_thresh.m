% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% An independent MATLAB implementation of the method used in VolPy:
%   Cai C, Friedrich J, Singh A, et al. (2021). VolPy: Automated and
%   scalable analysis pipelines for voltage imaging datasets.
%   PLoS Comput Biol 17(4): e1008806.
% --------------------------------------------------------------------------
function [thresh, falsePosRate, detectionRate, low_spikes] = adaptive_thresh(pks, clip, pnorm, min_spikes)

%     """ Adaptive threshold method for deciding threshold given heights of all peaks.
%     Args:
%         pks: 1-d array
%             height of all peaks
%         clip: int
%             maximum number of spikes for producing templates
%         pnorm: float, between 0 and 1, default is 0.5
%             a variable deciding the amount of spikes chosen for adaptive threshold method
%             
%         min_spikes: int
%             minimal number of spikes to be detected
%     Returns:
%         thresh: float
%             threshold for choosing spikes
%         falsePosRate: float
%             possibility of misclassify noise as real spikes
%         detectionRate: float
%             possibility of real spikes being detected
%         low_spikes: boolean
%             true if number of spikes is smaller than minimal value
%     """
% Method from VolPy (Cai et al., 2021); see the note at the top of this file.
% Yitong Qi, 2022

if ~exist('pnorm','var')
    pnorm = .5;
end
if ~exist('min_spikes', 'var')
    min_spikes = 10;
end

% find median of the kernel density estimation of peak heights
spread = [min(pks) max(pks)];
spread = spread + diff(spread) .* [-.05 .05];
low_spikes = false;
pts = linspace(spread(1), spread(2), 2001);
[f,xi] = ksdensity(pks, pts);
center = find(xi > median(pks),1,'first');

fmodel = [f(1:center+1) flip(f(1:center))]; % noise distribution
if length(fmodel) < length(f)
    fmodel = [fmodel ones(1, length(f) - length(fmodel)) * min(fmodel)];
else
    fmodel = fmodel(1:length(f));
end

% adjust the model so it doesn't exceed the data (YQ: this will force there
% to be a threshold)
csf = cumsum(f)/sum(f);
csmodel = cumsum(fmodel)/max(sum(f),sum(fmodel));
lastpt = find( (csf(1:end-1) > csmodel(1:end-1)+eps) & (csf(2:end) < csmodel(2:end)));
if isempty(lastpt)
    lastpt = center;
else
    lastpt = lastpt(1);
end
fmodel(1:lastpt+1) = f(1:lastpt+1);
fmodel(lastpt:end) = min(fmodel(lastpt:end),f(lastpt:end));

% find threshold
csf = cumsum(f);
csmodel = cumsum(fmodel);
csf2 = csf(end) - csf;
csmodel2 = csmodel(end) - csmodel;
obj = csf2.^pnorm - csmodel2.^pnorm;
[~,maxind] = max(obj);
thresh = xi(maxind);

if sum(pks > thresh) < min_spikes
    low_spikes = true;
    warning('Few spikes were detected. Adjusting threshold to take %i largest spikes',min_spikes);
    thresh = prctile(pks,100 * (1-min_spikes/length(pks)));
elseif (sum(pks > thresh) > clip) && (clip > 0)
    warning('Selecting top %i spikes for template',clip);
    thresh = prctile(pks, 100 * (1-clip/length(pks)));
end
[~,ix] = min(abs(xi - thresh));
falsePosRate = csmodel2(ix)/csf2(ix);
detectionRate = (csf2(ix) - csmodel2(ix))/max(csf2 - csmodel2);

