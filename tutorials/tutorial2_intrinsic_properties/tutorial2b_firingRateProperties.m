%% Tutorial 2b - maximum firing rate and pulse adaptation
%
% Two properties, both read off the trial-averaged spike train and neither
% needing the voltage trace at all. This is the cheapest stage: the spike lists
% are already in memory after the load, so it costs no trace I/O.
%
%   frMax        the firing rate over the last 635 ms of the ramp, when the
%                drive is at its maximum. An upper bound on sustained output.
%
%   frPulseAdp   within a single 500 ms pulse, how much the rate falls from the
%                first 50 ms to the last 250 ms, as a fraction of the late
%                rate. Spike-frequency adaptation on a timescale the ramp is
%                too slow to expose.
%
% Both are rates, so neither involves the spike-height normalisation at all;
% they would be identical on raw traces.
%
% Runtime: about 30 s, all of it the load.

close all;

% Set nwbFile before running to point at your copy of a sparseOpto session,
% or edit the placeholder below. Both deposited sparseOpto sessions work.
if ~exist('nwbFile', 'var') || isempty(nwbFile)
    nwbFile = 'path/to/M_YQ0201_27_sparseOpto.nwb';
end

data = loadSparseOptoForEphys(nwbFile);
protocol = blueProtocol(data.blue, data.dt);
stacks = stackStimTrials(data, 'Subthreshold', 'none');

dt = data.dt;

%% -------------------------------------------------------- the two windows
% 50 ms at pulse onset against 250 ms at pulse offset. The asymmetry is
% deliberate: the onset window has to be short enough to catch the transient
% before adaptation sets in, while the offset window can be long because the
% rate there is by then stationary, and a longer window buys precision.
nFrameAvg1 = round(50e-3/dt);     % 39 frames
nFrameAvg2 = round(250e-3/dt);    % 197 frames
fprintf('\nAdaptation windows: first %d frames (%.1f ms), last %d frames (%.1f ms)\n', ...
    nFrameAvg1, nFrameAvg1*dt*1e3, nFrameAvg2, nFrameAvg2*dt*1e3);

nPulse = numel(protocol.pulseAmp);
frPulseStart = zeros(data.nCells, nPulse);
frPulseEnd = zeros(data.nCells, nPulse);
for jj = 1:nPulse
    % The source locates these with find(diff(blue) == blueVs(jj)), an exact
    % float equality against the amplitude. That works, but it means the
    % window indices depend on the amplitudes being distinct. Taking them from
    % the segmented protocol instead is equivalent here and does not.
    onsetWindow  = protocol.pulseOnset(jj)  + (0:nFrameAvg1-1);
    offsetWindow = protocol.pulseOffset(jj) - (nFrameAvg2-1:-1:0);
    frPulseStart(:, jj) = mean(stacks.fr_mat_stim_avg(:, onsetWindow), 2);
    frPulseEnd(:, jj)   = mean(stacks.fr_mat_stim_avg(:, offsetWindow), 2);
end

frPulseAdpAll = (frPulseStart - frPulseEnd) ./ frPulseEnd;

% The source keeps column 3, the third pulse. That is the 0.793 amplitude
% step: high enough that most cells fire through the whole pulse, so the late
% rate in the denominator is not near zero, but below saturation.
adpPulse = 3;
frPulseAdp = frPulseAdpAll(:, adpPulse);
fprintf('Adaptation taken from pulse %d, amplitude %.4f\n', ...
    adpPulse, protocol.pulseAmp(adpPulse));

% Two degenerate cases are worth naming, because both land on a finite-looking
% number. A cell that fires nothing in the early window gives exactly -1: that
% is the floor of the measure, reached whenever the cell is silent at pulse
% onset, and it does not distinguish "slightly negative adaptation" from
% "no response at all". A cell silent in the LATE window puts a zero in the
% denominator and gives Inf or NaN.
floored = frPulseStart(:, adpPulse) == 0;
degenerate = frPulseEnd(:, adpPulse) == 0;
fprintf('  %d cells are silent in the early window, so frPulseAdp is exactly -1 (%d pass QC)\n', ...
    sum(floored), sum(floored & data.idxUse));
fprintf('  %d cells are silent in the late window, giving a zero denominator (%d pass QC)\n', ...
    sum(degenerate), sum(degenerate & data.idxUse));

%% --------------------------------------------------------- maximum rate
% The source finds the ramp offset with
%   tFrMax = find(diff(blue.*(1:nFrame>4000)) < -blueVs(end-1));
% which is the ramp's own down-step, masked so the pulses cannot match. That
% is protocol.rampOff.
maxWindow = protocol.rampOff + (-500:0);
frMax = mean(stacks.fr_mat_stim_avg(:, maxWindow), 2) / dt;
fprintf('\nMaximum rate over the last %d frames (%.0f ms) of the ramp\n', ...
    numel(maxWindow), numel(maxWindow)*dt*1e3);

%% ------------------------------------------------------- against the deposit
% Neither property depends on the voltage trace or on the spike set outside
% the stimulus window, which is why both reproduce the deposit exactly despite
% the missing motion correction.
report = @(name, mine, dep) fprintf( ...
    '  %-12s max|diff| %.3e   r %.6f   %5.1f%% within 1e-6\n', ...
    name, max(abs(mine - dep)), corr(mine, dep), 100*mean(abs(mine - dep) < 1e-6));

fprintf('\nAgainst the deposited unit columns (all %d cells)\n', data.nCells);
report('frMax', frMax, data.deposited.frMax);
report('frPulseAdp', frPulseAdp, data.deposited.frPulseAdp);

%% ----------------------------------------------------------------- plots
plotIntrinsicPropertyDiagnostics('rates', data, protocol, stacks, ...
    struct('frMax', frMax, 'frPulseAdp', frPulseAdp, ...
           'frPulseStart', frPulseStart, 'frPulseEnd', frPulseEnd, ...
           'adpPulse', adpPulse));

fprintf('\nNext: tutorial2c_rheobase.m\n');
