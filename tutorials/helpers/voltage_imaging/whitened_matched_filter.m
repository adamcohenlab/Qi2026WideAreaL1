% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Author: Yitong Qi.
% Origin: the author's MATLAB voltage-imaging library (imaging/voltage).
% An independent MATLAB implementation of the method used in VolPy:
%   Cai C, Friedrich J, Singh A, et al. (2021). VolPy: Automated and
%   scalable analysis pipelines for voltage imaging datasets.
%   PLoS Comput Biol 17(4): e1008806.
% --------------------------------------------------------------------------
function datafilt = whitened_matched_filter(data, locs, window)
% """
%     Function for using whitened matched filter to the original signal for better
%     SNR. Use welch method to approximate the spectral density of the signal.
%     Rescale the signal in frequency domain. After scaling, convolve the signal with
%     peak-triggered-average to make spikes more prominent.
%     
%     Args:
%         data: 1-d array
%             input signal
%         locs: 1-d array
%             spike times
%         window: 1-d array
%             window with size of temporal filter
%     Returns:
%         datafilt: 1-d array
%             signal processed after whitened matched filter
%     
%     """
% Method from VolPy (Cai et al., 2021); see the note at the top of this file.
% Yitong Qi, 2022

N = ceil(log2(length(data)));
censor = zeros(1,length(data));
censor(locs) = 1;
censor = conv(censor,ones(1,length(window)),'same');
censor = censor < .5;
noise = data(censor);

pxx = pwelch(noise,2^N/16,[],2^N,2*pi);

Nf2 = [reshape(pxx,1,[]) flip(reshape(pxx(2:end-1),1,[]))];
scaling_vector = 1./sqrt(Nf2);

cc = padarray(reshape(data,1,[]) ,[0 2^N - length(data)], 0, 'post');
dd = fft(cc) .* scaling_vector;
dataScaled = real(ifft(dd));
PTDscaled = dataScaled(reshape(locs,[],1) + reshape(window,1,[]));
PTAscaled = mean(PTDscaled,1);
datafilt = conv(dataScaled,flip(PTAscaled),'same');
% datafilt = conv(dataScaled,flip(PTAscaled-mean(PTAscaled)),'same'); % for comparison purpose only
datafilt = datafilt(1:length(data));

%% 
