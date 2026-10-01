function out = estimateSpikesPerPulse(traceRow, tBlue, prm)
%ESTIMATESPIKESPERPULSE Spikes evoked per blue pulse, from spectral width.
%
%   out = ESTIMATESPIKESPERPULSE(traceRow, tBlue, prm) estimates the mean number
%   of spikes a single blue pulse evokes in one cell, WITHOUT detecting spikes.
%
%   Why not just count spikes. During the 20 ms pulse the cell sits on a large
%   optogenetic depolarisation plus optical crosstalk. Spike detection there is
%   unreliable, which is the same fact that makes the crosstalk template's rise
%   constant untrustworthy (see tutorial3c). So the rate is inferred from the
%   residual's temporal structure instead of from events.
%
%   The idea. Subtract the stereotyped, trial-invariant pulse response and what
%   remains, within the pulse window, is mostly spikes. Spikes are brief, so the
%   more of them there are, the broader the residual's power spectrum. The
%   spectral centroid of the residual autocovariance is therefore a monotone
%   proxy for firing rate; multiplying by the pulse duration converts rate to a
%   count.
%
%   INPUTS
%     traceRow  1 x nFrames voltage for one cell. Scale-free: see SCALE below.
%     tBlue     vector of pulse trigger frames for that cell, convention A.
%     prm       struct with fields tPrev, tOn, tOff, nBack, nFront, dt, tPulse.
%
%   OUTPUT
%     out.nSpikesPulse  the estimate
%     out.P             one-sided spectrum after median subtraction
%     out.autoMat       per-trial autocovariance, for plotting
%     out.stMat         per-trial pulse-triggered matrix
%     out.template      the fitted stereotyped response
%     out.nTrials
%
%   SCALE. Every step -- median across trials, polyfit, smooth, regression
%   residual, xcov, median subtraction, and the centroid ratio -- is positively
%   homogeneous, and the final quantity is a RATIO of spectra. The per-cell
%   spkHgtNoBlue factor therefore cancels exactly, so this returns the same
%   answer from deposited normalized_voltage as from the source traces_all_bc.
%
%   This implements the block headed "%% Correct version" in
%   bluePulseSpikeCount.m, which is the one whose output reached results.mat.
%   The earlier block in the same file uses a two-sided fft2 across trials and a
%   different frequency scale, and gives different numbers; see
%   TUTORIAL3_VALIDATION.md.
%
%   Requires get_sta_mat_raw, fitSelfBlueSta and SeeResiduals_vec on the path.

nLag = prm.tOff - prm.tOn + 1;

% Per-trial pulse-triggered matrix. get_sta_mat_raw returns the raw trials, not
% an average, unlike its get_sta_mat_* siblings.
stMat = squeeze(get_sta_mat_raw({tBlue}, traceRow(:), [prm.nBack prm.nFront], 1));

% The stereotyped response: flat before the pulse, cubic during it, smoothed
% after. Fitted to the across-trial median so spikes do not drag it.
template = fitSelfBlueSta(stMat, prm.tPrev, prm.tOn, prm.tOff);

% Regress that template out of every trial. What is left in the pulse window is
% spiking plus noise.
stMatCor = SeeResiduals_vec(stMat, template);

nTrials = size(stMat,1);
autoMat = zeros(nTrials, 2*(prm.tOff - prm.tOn) + 3);
for t = 1:nTrials
    s1 = stMatCor(t, prm.tOn-1:prm.tOff);
    autoMat(t,:) = xcov(s1, s1);
end

% Average the per-trial autocovariance spectra, drop the broadband floor, and
% keep the one-sided half.
P = abs(mean(fft(autoMat, [], 2), 1));
Pcor = max(P - median(P), 0);
Pcor = Pcor(1:nLag+1);

% Spectral centroid in Hz, times pulse duration, gives spikes per pulse.
freq = (0:nLag) * (1/prm.dt/(nLag+1)/2);
out.nSpikesPulse = sum(freq .* Pcor) / sum(Pcor) * prm.tPulse;

out.P = Pcor;
out.freq = freq;
out.autoMat = autoMat;
out.stMat = stMat;
out.stMatCor = stMatCor;
out.template = template;
out.nTrials = nTrials;
end
