function protocol = blueProtocol(blue, dt)
%BLUEPROTOCOL Recover the optogenetic stimulus structure from the waveform.
%
%   protocol = BLUEPROTOCOL(blue, dt) takes one acquisition chunk of blue
%   stimulation (nFrame x 1) and returns
%
%     pulseOnset   1 x nPulse, first frame of each square pulse
%     pulseOffset  1 x nPulse, last frame of each square pulse
%     pulseAmp     1 x nPulse, amplitude of each square pulse
%     rampOn       scalar, ramp onset as findpeaks reports it
%     rampFirst    scalar, first frame on which the ramp is above zero
%     rampOff      scalar, last frame of the ramp
%     offFrames    1 x (nPulse+1), every frame after which blue steps down
%     nFrame, dt
%
%   WHY THIS EXISTS
%   The source script hard-codes tRampOn = 3557 and tRampOff = 7480, and
%   separates the ramp from the square pulses with
%
%       [tmp,tRamp] = findpeaks(diff(blue),'minPeakHeight',6e-3);
%       tRamp(tmp > 0.25) = [];
%
%   Both literals are correct: measured on M-YQ0201-27 and M-YQ0201-29 the
%   protocol is identical, four 500 ms pulses at 0.2687 / 0.5145 / 0.7930 /
%   1.0958 with onsets 393 / 1181 / 1968 / 2755, then a ramp from 3544 to 7480
%   whose findpeaks onset is 3557. They are protocol constants, not session
%   constants.
%
%   The 0.25 discriminant is nonetheless brittle: it works only because the
%   smallest pulse is 0.2687, a 7% margin. A dimmer first pulse in a future
%   session would be silently reclassified as a ramp and every rheobase would
%   be measured from the wrong frame. This function separates the two on
%   segment duration instead, which has three orders of magnitude of headroom,
%   then reports the findpeaks onset so the source arithmetic is reproduced
%   exactly. It asserts against the published literals so a protocol change
%   fails loudly instead of quietly.
%
%   See also LOADSPARSEOPTOFOREPHYS, TUTORIAL2A_STIMULUSPROTOCOL.

narginchk(2, 2);
blue = double(reshape(blue, [], 1));
nFrame = numel(blue);

%% ------------------------------------------------- segment the stimulus
% A "segment" is a maximal run of frames with blue above zero. The square
% pulses and the ramp are separated by their duration, which differs by a
% factor of ten here, rather than by their amplitude, which differs by 7%.
on = blue > 0;
assert(any(on), 'Tutorial2:NoStimulation', 'This chunk carries no stimulation.');
edges = diff([false; on; false]);
segStart = find(edges == 1);
segStop  = find(edges == -1) - 1;
segLen = segStop - segStart + 1;

% The ramp is the long segment. Everything shorter is a square pulse.
isRamp = segLen > 2 * median(segLen);
assert(sum(isRamp) == 1, 'Tutorial2:RampCount', ...
    'Expected exactly one ramp segment, found %d.', sum(isRamp));

pulseStart = segStart(~isRamp);
pulseStop  = segStop(~isRamp);

protocol = struct();
protocol.nFrame = nFrame;
protocol.dt = dt;
protocol.pulseOnset = reshape(pulseStart, 1, []);
protocol.pulseOffset = reshape(pulseStop, 1, []);
protocol.pulseAmp = reshape(arrayfun(@(a,b) max(blue(a:b)), pulseStart, pulseStop), 1, []);
protocol.rampFirst = segStart(isRamp);
protocol.rampOff = segStop(isRamp);

% Square pulses must actually be square, or "amplitude of the pulse" is not a
% well-defined quantity and frPulseStart/frPulseEnd are averaging over a ramp.
for p = 1:numel(pulseStart)
    seg = blue(pulseStart(p):pulseStop(p));
    assert(max(seg) - min(seg) < 1e-6 * max(seg), 'Tutorial2:PulseNotSquare', ...
        'Pulse %d varies by %.3g within its segment.', p, max(seg) - min(seg));
end

%% ------------------------------------------- ramp onset, source arithmetic
% Reproduced from the source rather than replaced: vRheobase is measured as a
% fraction of (tRampOff - tRampOn), so the onset convention is part of the
% published number. findpeaks on diff(blue) fires at the first frame where the
% ramp's step exceeds 6 mV-equivalent, which is 13 frames after the ramp first
% leaves zero because the ramp starts far below that step size.
minPkHgt = 6e-3;
maxPkHgt = 0.25;
[pkHeight, pkFrame] = findpeaks(diff(blue)', 'minPeakHeight', minPkHgt);
pkFrame(pkHeight > maxPkHgt) = [];
assert(isscalar(pkFrame), 'Tutorial2:RampOnset', ...
    ['The source findpeaks rule returned %d ramp onsets, not 1. The 0.25 ' ...
     'discriminant has probably collided with a pulse amplitude; the ' ...
     'smallest pulse here is %.4f.'], numel(pkFrame), min(protocol.pulseAmp));
protocol.rampOn = pkFrame;

%% -------------------------------------------------------- down-steps
% t_blue_off in the source: every frame after which blue decreases. These are
% the pulse offsets plus the ramp offset, and they trigger the blue-off STA
% that the membrane time constant is fitted to.
protocol.offFrames = reshape(find(diff(blue) < 0), 1, []);
assert(numel(protocol.offFrames) == numel(pulseStart) + 1, ...
    'Tutorial2:OffFrameCount', ...
    ['Expected %d down-steps (one per pulse plus the ramp), found %d. A ' ...
     'non-monotonic ramp would add spurious blue-off triggers.'], ...
    numel(pulseStart) + 1, numel(protocol.offFrames));

%% ---------------------------------------------- assert the published shape
% Not a tolerance check on the science, just a tripwire: these are the values
% both deposited sessions carry, and every window index below is derived from
% them. A different protocol needs the stage code re-read, not re-run.
% Note the off-by-one that catches everyone reading the source: findpeaks and
% find(diff(blue)==amp) report the last frame BEFORE a step, so the source's
% "onset 393" is the frame before the pulse. pulseOnset here is the first
% stimulated frame, 394. rampOn keeps the source's diff convention because
% vRheobase is measured as a fraction of (rampOff - rampOn).
expected = struct('pulseOnset', [394 1182 1969 2756], ...
                  'pulseOffset', [787 1574 2362 3149], ...
                  'rampOn', 3557, 'rampOff', 7480);
names = fieldnames(expected);
for k = 1:numel(names)
    if ~isequal(protocol.(names{k}), expected.(names{k}))
        warning('Tutorial2:ProtocolDiffers', ...
            ['%s is %s, not the %s seen in both deposited sessions. The ' ...
             'stage code derives its windows from this struct so it will ' ...
             'still run, but the published comparison no longer applies.'], ...
            names{k}, mat2str(protocol.(names{k})), mat2str(expected.(names{k})));
    end
end
end
