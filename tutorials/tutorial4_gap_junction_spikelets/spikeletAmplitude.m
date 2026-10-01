function amp = spikeletAmplitude(sta, nBack)
%SPIKELETAMPLITUDE Peak-minus-flank amplitude of a spikelet waveform.
%
%   amp = SPIKELETAMPLITUDE(sta, nBack) takes an array whose LAST dimension
%   is the time offset, with offset zero at index nBack+1, and returns the
%   source definition:
%
%       mean over offsets [0 1]  -  mean over offsets [-4 -3 -2 4 5 6]
%
%   The output has the size of sta with the last dimension removed.
%
%   WHY THIS DEFINITION
%   A spikelet through a gap junction is a low-pass filtered copy of the
%   presynaptic spike: it peaks within a frame or two of the pre-spike peak
%   and is a few frames wide at 1.27 ms sampling. Two samples at the peak
%   capture it without being at the mercy of a one-frame alignment error.
%
%   The flank average is the local baseline. It straddles the peak, three
%   samples each side, with a two-sample gap left on either side so the
%   shoulders of the spikelet itself do not enter the baseline. Taking it
%   symmetrically is what makes the measure insensitive to a linear drift
%   through the window, which matters because the drift is exactly what the
%   common-mode correction in stage 4c is removing and neither correction is
%   perfect.
%
%   The definition is applied identically to the observed waveform and to
%   every null waveform. That is the only thing that has to be true for the
%   p-values to mean something, and it is why the tutorial keeps it in one
%   function rather than writing the index arithmetic out four times.
%
%   NOTE
%   The flank offsets are fixed at -4:-2 and 4:6 regardless of the window
%   half-width, so nBack must be at least 4 and nFront at least 6. The source
%   uses a +/-20 window for the observed waveform and +/-10 for the nulls;
%   both satisfy that, and the amplitude is unaffected by the surplus.
%
%   See also CORRECTCOMMONMODE, SPIKELETNULLSET.

nTau = size(sta, ndims(sta));
peakIdx  = nBack + 1 + (0:1);
flankIdx = nBack + 1 + [-4:-2, 4:6];
assert(min(flankIdx) >= 1 && max(flankIdx) <= nTau, 'Tutorial4:WindowTooShort', ...
    'Window of %d samples with nBack = %d cannot hold the -4:6 flanks.', nTau, nBack);

sz = size(sta);
flat = reshape(sta, [], nTau);
amp = mean(flat(:, peakIdx), 2) - mean(flat(:, flankIdx), 2);

if numel(sz) > 2
    amp = reshape(amp, sz(1:end-1));
else
    amp = reshape(amp, sz(1), 1);
end
end
