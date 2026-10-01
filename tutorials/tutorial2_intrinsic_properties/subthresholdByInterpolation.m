function sub = subthresholdByInterpolation(trace, spikeFrames, spkWidth)
%SUBTHRESHOLDBYINTERPOLATION Remove spikes by linear interpolation across them.
%
%   sub = SUBTHRESHOLDBYINTERPOLATION(trace, spikeFrames) cuts a window around
%   every spike out of trace and fills the gap by linear interpolation from the
%   surviving samples, reproducing
%
%       t_i = setdiff(1:nFramesTotal, spk_t{ii} + spkWidth');
%       traces_all_sub(ii,:) = interp1(t_i, traces_all(ii,t_i), 1:nFramesTotal, 'linear');
%
%   from sparseOpto_ephysProp_meas_cleanup.m. spkWidth defaults to -3:5.
%
%   trace is 1 x nT, spikeFrames is a vector of 1-based indices into it, and
%   sub is 1 x nT of the same class as trace.
%
%   WHY INTERPOLATION AND NOT THE DEPOSITED SUBTHRESHOLD
%   The deposit's traces_all_sub_n is a different estimate, built in
%   dataPrepsForMAT.m by replacing the spike window with a 16-frame running
%   median and then applying sgolayfilt(...,2,17). Rheobase does not want that
%   one. It averages ten stimulus repeats before reading the trace, and that
%   average already removes the noise the Savitzky-Golay pass targets; applying
%   a 17-frame polynomial filter as well only distorts the slow ramp the cubic
%   fit is there to measure. Plain interpolation is the estimate the property
%   was defined on.
%
%   EDGE BEHAVIOUR
%   interp1 returns NaN outside the span of the retained samples, so a spike
%   inside the first or last spkWidth samples leaves NaN at the edge. The
%   source has the same behaviour over the whole session; here the function is
%   applied per acquisition chunk, so it can occur at each chunk edge instead
%   of only twice per session. Those frames are 8 at most and lie inside the
%   chunk-boundary artifact region the stimulus never reaches, but the count is
%   returned through a warning rather than being silently filled.
%
%   CHUNK-LOCAL VERSUS SESSION-GLOBAL
%   Interpolation is local: a sample's value depends only on the nearest
%   retained sample either side. Applying it per chunk therefore gives exactly
%   the session-global result except within spkWidth of a chunk edge. The
%   stimulus occupies frames 394 to 7480 of 8000, so the difference cannot
%   reach any window this tutorial measures. tutorial2c asserts this on a
%   sample of cells rather than asserting it here on every call.
%
%   See also TUTORIAL2C_RHEOBASE.

narginchk(2, 3);
if nargin < 3 || isempty(spkWidth)
    spkWidth = -3:5;
end
trace = reshape(trace, 1, []);
nT = numel(trace);

% spk_t{ii} + spkWidth' is a numel(spkWidth) x nSpk matrix in the source; the
% orientation does not matter once it is flattened for setdiff.
cut = reshape(spikeFrames(:)' + spkWidth(:), 1, []);
cut = cut(cut >= 1 & cut <= nT);
keep = setdiff(1:nT, cut);

assert(numel(keep) >= 2, 'Tutorial2:NoSubthresholdSamples', ...
    'Spike windows cover all but %d of %d samples; nothing to interpolate from.', ...
    numel(keep), nT);

sub = interp1(keep, double(trace(keep)), 1:nT, 'linear');

nEdgeNaN = sum(isnan(sub));
if nEdgeNaN > 0
    warning('Tutorial2:SubthresholdEdgeNaN', ...
        ['%d frame(s) at the trace edge are outside the interpolation span ' ...
         'and are NaN. Expected only when a spike lands within %d frames of ' ...
         'an edge.'], nEdgeNaN, max(abs(spkWidth)));
end

sub = cast(sub, 'like', trace);
end
