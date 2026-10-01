%% Tutorial 2a - the stimulus, and what a "trial" is here
%
% Intrinsic properties are what a cell does to an input you control, so before
% measuring any of them you have to know exactly what the input was and which
% repeats belong to which cell. That is the whole of this stage. It produces
% the scaffold the other three stages read, and it establishes four facts that
% the source script assumes without stating:
%
%   1. One stimulus waveform is shared by every cell and every repeat. The
%      "sparse" in sparseOpto is which cells are lit in a given chunk, not what
%      they are lit with.
%   2. The protocol is four 500 ms square pulses of rising amplitude, then a
%      5 s ramp. Pulses measure adaptation and maximum rate; the ramp measures
%      rheobase; the down-steps measure the membrane time constant.
%   3. Every analysis window fits inside one 8000-frame acquisition chunk, so
%      nothing has to be stitched across a chunk boundary.
%   4. The chunk-boundary artifact mask t_noise cannot touch any of them.
%
% Runtime: about 60 s, almost all of it the one streaming pass over blue.
%
% See README.md for the I/O note before you rewrite any of this as a per-cell
% loop.

close all;

% Set nwbFile before running to point at your copy of a sparseOpto session,
% or edit the placeholder below. Both deposited sparseOpto sessions work.
if ~exist('nwbFile', 'var') || isempty(nwbFile)
    nwbFile = 'path/to/M_YQ0201_27_sparseOpto.nwb';
end

%% ------------------------------------------------------------------ load
% ScanBlue makes one pass over the 3 GB blue series. It returns which cells
% were stimulated in which chunk, and the stimulus waveform itself, asserting
% on the way that there is only one waveform to return. After this, nothing
% downstream reads blue again.
data = loadSparseOptoForEphys(nwbFile);

%% -------------------------------------------------- recover the protocol
protocol = blueProtocol(data.blue, data.dt);

fprintf('\nProtocol\n');
fprintf('  %d square pulses\n', numel(protocol.pulseAmp));
for p = 1:numel(protocol.pulseAmp)
    fprintf('    %d: frames %4d-%4d (%5.1f ms) at amplitude %.4f\n', p, ...
        protocol.pulseOnset(p), protocol.pulseOffset(p), ...
        (protocol.pulseOffset(p) - protocol.pulseOnset(p) + 1) * data.dt * 1e3, ...
        protocol.pulseAmp(p));
end
fprintf('  ramp frames %d-%d (%.0f ms); findpeaks onset %d\n', ...
    protocol.rampFirst, protocol.rampOff, ...
    (protocol.rampOff - protocol.rampFirst + 1) * data.dt * 1e3, protocol.rampOn);
fprintf('  down-steps at %s\n', mat2str(protocol.offFrames));

% The source's nFrameRamp = round(5/dt) is the ramp length in frames. It is
% computed in both rheobase blocks and never used; the value is confirmed here
% so a reader does not go looking for what it does.
fprintf('  source nFrameRamp = %d frames; measured ramp = %d frames (unused in source)\n', ...
    round(5/data.dt), protocol.rampOff - protocol.rampFirst + 1);

%% --------------------------------------- fact 3: windows stay in a chunk
% The two triggered averages this tutorial computes use windows of [-140 +140]
% around an evoked spike and [-140 +390] around a blue down-step. Both have to
% fit between frame 1 and frame nFrame of the chunk that holds the trigger.
spikeSpan = [min(protocol.pulseOnset) - 140, protocol.rampOff + 140];
offSpan   = [min(protocol.offFrames) - 140, max(protocol.offFrames) + 390];
fprintf('\nWindow containment (chunk is 1..%d)\n', data.nFrame);
fprintf('  spike STA  spans %d..%d\n', spikeSpan(1), spikeSpan(2));
fprintf('  blue-off STA spans %d..%d\n', offSpan(1), offSpan(2));
assert(spikeSpan(1) >= 1 && spikeSpan(2) <= data.nFrame && ...
       offSpan(1) >= 1 && offSpan(2) <= data.nFrame, ...
    'A triggered-average window leaves its acquisition chunk.');
fprintf('  both inside the chunk, so no boundary stitching is needed\n');

%% ------------------------------------------ fact 4: t_noise is irrelevant
% The source strips t_noise = (0:nFrame:T) + (-140:140) from the spike matrix
% before computing firing rates. In chunk-local terms that is frames 1:140 and
% 7861:8000. The stimulus never reaches either, so the mask is a no-op for
% every property measured here and is not implemented.
tNoiseLocal = unique(mod((0:data.nFrame:data.T)' + (-140:140) - 1, data.nFrame) + 1);
stimFrames = find(data.blue > 0);
overlap = intersect(tNoiseLocal, stimFrames);
fprintf('\nt_noise covers %d local frames; stimulus covers %d..%d; overlap %d frames\n', ...
    numel(tNoiseLocal), min(stimFrames), max(stimFrames), numel(overlap));
assert(isempty(overlap), 't_noise intersects the stimulus; it can no longer be skipped.');

%% ------------------------------------------------------- stack the trials
% One pass, all 320 cells. 'none' skips the voltage read entirely: this stage
% and 2b need only the spike-locked average, which is built from spike lists
% already in memory. 2c and 2d pay for the voltage.
stacks = stackStimTrials(data, 'Subthreshold', 'none');

fprintf('\nRepeat counts\n');
fprintf('  the source pre-allocates nTrial = 10; observed max is %d, min %d\n', ...
    max(stacks.nTrial), min(stacks.nTrial));
fprintf('  so %d of %d cells rest on fewer than 10 repeats\n', ...
    sum(stacks.nTrial < 10), data.nCells);

%% ----------------------------------------------------------------- plots
plotIntrinsicPropertyDiagnostics('protocol', data, protocol, stacks);

%% ------------------------------------------------------------ hand-off
save(fullfile(tempdir, 'tutorial2a_scaffold.mat'), 'protocol', '-v7.3');
fprintf('\nProtocol saved to %s\n', fullfile(tempdir, 'tutorial2a_scaffold.mat'));
fprintf('Next: tutorial2b_firingRateProperties.m\n');
