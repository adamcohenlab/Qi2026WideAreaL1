%% Tutorial 4b — the spikelet waveform, and why the raw amplitude is not usable
%
% Stage 4a produced trigger sets that are clean in the sense that matters:
% the post cell was dark when the pre cell fired. Average the post cell's
% trace on them and a coupled pair shows a small positive deflection at
% offset zero, a few frames wide. That is the spikelet.
%
% This stage computes those averages, defines the amplitude, and then shows
% the problem that stage 4c exists to solve: the average is NOT flat for
% uncoupled pairs. Two artifacts survive orthogonality, and one of them is
% strongly distance dependent, which is fatal because distance is exactly
% the variable a connectivity result is read against.
%
% Runtime: about 2.5 min, most of it the one-off spike subtraction.

clear; close all

%% ----------------------------------------------------------- configuration
% Helper functions bundled with this repository (tutorials/helpers).
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));

cfg = struct();
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
cfg.nBack = 20;
cfg.nFront = 20;
cfg.examplePre = 152;

%% ------------------------------------------- 1. load, triggers, subtraction
data = loadSparsePulseHadForSpikelet(cfg.nwbFile);
trig = buildOrthogonalTriggers(data);

% The post cell's own spikes are removed first. They are tens of times
% larger than a spikelet, and a post spike that happens to land near a pre
% spike would dominate the average. What is left is the subthreshold trace,
% which is where an electrical coupling shows up.
fprintf('\nspike subtraction (one-off, about 70 s):\n');
tracesSub = subtractSpikes(data.traces_all, data.spk_t);
data.traces_all = [];                       % free 4 GB; not needed below
traceMat = tracesSub';                      % T x nCells, kept single: 4 GB, not 8
clear tracesSub

%% ----------------------------------------- 2. one presynaptic cell's row
staRow = spikeletStaRow(cfg.examplePre, trig, traceMat, ...
    'nBack', cfg.nBack, 'nFront', cfg.nFront, 'Verbose', true);
tau = -cfg.nBack:cfg.nFront;
distRow = data.distMat(cfg.examplePre,:);

% Every post cell, coloured near to far, with the amplitude windows shaded:
% two samples at the peak in red, the six flank samples in blue.
plotSpikeletDiagnostics('waveformRow', tau, staRow, distRow, data.dx, ...
    cfg.examplePre, data.dt, cfg.nBack);

%% ------------------------------------------------- 3. the amplitude measure
%   amplitude = mean over offsets [0 1] - mean over offsets [-4 -3 -2 4 5 6]
%
% Two samples at the peak, because a gap-junction spikelet is low-pass
% filtered and a couple of frames wide at 1.27 ms, and because two samples
% survive a one-frame alignment error. A symmetric flank baseline, with a
% two-sample gap on each side so the spikelet's own shoulders stay out of
% it; symmetric so that a linear drift through the window cancels.
ampRaw = spikeletAmplitude(staRow - mean(staRow,2), cfg.nBack);

%% --------------------------------- 4. the problem: amplitude versus distance
% If the average were clean, an uncoupled pair would sit at zero and only
% genuinely coupled pairs would lift off it. Instead there is a smooth
% positive trend that grows towards short distances, and a pedestal that
% does not vanish even at the far edge of the field.
fprintf('\n--- raw amplitude versus distance, pre cell %d ---\n', cfg.examplePre);
edges = [0 100 200 400 800 1600 3000];
d = distRow(:) * data.dx;
notSelf = (1:data.nCells)' ~= cfg.examplePre;
for k = 1:numel(edges)-1
    m = d >= edges(k) & d < edges(k+1) & notSelf;
    if any(m)
        fprintf('  %4d - %4d um : n = %3d, median amplitude %+.3e\n', ...
            edges(k), edges(k+1), sum(m), median(ampRaw(m)));
    end
end

% Two artifacts, with different shapes, and it is worth separating them
% because stage 4c removes them in two different ways.
%
% 1. A smooth, slow, DISTANCE-DEPENDENT trend. Optical crosstalk and
%    scattered light from the pre cell's own stimulation reach the post
%    cell, and how much depends on how far away it is. Orthogonality does
%    not touch this: the pre cell really is lit, and its light really does
%    spread. This is the dangerous one, because it is largest exactly where
%    coupling is most likely.
%
% 2. A largely DISTANCE-INDEPENDENT common mode. The pre cell's spike is
%    not an isolated event; it is correlated with network activity, and the
%    whole field shares optics and illumination. Cells past 800 um are too
%    far to be coupled, so whatever they still show at offset zero is
%    mostly this.
far = d > 800 & notSelf;
veryFar = d > 1600 & notSelf;
fprintf('\n  far-field (> 800 um) median amplitude: %+.3e over %d cells\n', ...
    median(ampRaw(far)), sum(far));
fprintf('  beyond 1600 um                       : %+.3e over %d cells\n', ...
    median(ampRaw(veryFar)), sum(veryFar));
fprintf('  this pedestal is what the far-field average in stage 4c estimates.\n');
%
% Note the two numbers are not equal. The 800 um cut is the pipeline's
% operational choice, and at that radius the distance trend has flattened
% but has not vanished, so the "common mode" estimated there still carries
% a little of the trend's tail. A stricter cut would estimate a smaller
% pedestal on fewer cells. This matters less than it looks, because stage
% 4c fits a per-post-cell coefficient on that waveform rather than
% subtracting it outright, so the shape is what has to be right, not the
% scale.

% The far-field waveform itself, which is the thing stage 4c regresses out.
figure('Color','w','Name','the two artifacts');
plot(tau*data.dt*1e3, mean(staRow(far,:) - mean(staRow(far,:),2), 1), 'k', 'LineWidth', 2);
hold on
nearMask = d > 0 & d < 200 & notSelf;
plot(tau*data.dt*1e3, mean(staRow(nearMask,:) - mean(staRow(nearMask,:),2), 1), ...
    'r', 'LineWidth', 2);
xlabel('Peri-spike time (ms)'); ylabel('Normalized spike height')
legend({'> 800 \mum (common mode only)', '< 200 \mum (common mode + trend + coupling)'}, ...
    'Box','off','Location','northwest')
title(sprintf('pre cell %d', cfg.examplePre))

%% ----------------------------------------------- 5. why not just subtract it
% It is tempting to subtract the far-field average from everything and stop.
% That removes artifact 2 but not artifact 1, and artifact 1 is the one that
% masquerades as connectivity: it is monotonic in distance, so subtracting a
% distance-independent constant leaves a distance-dependent residue that any
% connectivity-versus-distance plot will read as real coupling.
%
% Stage 4c therefore does both, in the order that matters, and excludes the
% peak region from both fits so that neither correction can eat the signal.
%
% Keep this stage's outputs: 4c re-uses them.
save(fullfile(tempdir, 'tutorial4b_state.mat'), 'cfg', 'staRow', 'ampRaw', ...
    'tau', 'distRow', '-v7.3');
fprintf('\nStage 4b done. Next: tutorial4c_commonModeCorrection.\n');
