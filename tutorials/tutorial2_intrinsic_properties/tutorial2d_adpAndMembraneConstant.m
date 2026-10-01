%% Tutorial 2d - afterdepolarisation and membrane time constant
%
% The two waveform properties. Both are read off a spike-height-normalised
% triggered average, and both are exactly the kind of quantity that makes the
% normalisation worth having: they are voltages, so they need a scale, and the
% only scale available in an optical recording is the cell's own spike.
%
%   adp_blue    the afterdepolarisation. Average the cell's voltage around its
%               own optogenetically evoked spikes, then subtract a pre-spike
%               baseline from a post-spike window. What is left is the bump
%               that follows the spike.
%
%   membraneC   the membrane time constant. Average the voltage around the
%               moment the blue light steps down, and fit a single exponential
%               to the relaxation.
%
% Both averages are computed by streaming over acquisition chunks rather than
% by pulling whole cell traces, for the reason in README.md.
%
% Runtime: about 40 s.

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

%% ================================================== afterdepolarisation
%
% Which spikes. Only blue-evoked ones, so that the waveform is measured under
% a controlled drive, and only those with no other spike within 20 frames
% either side. The ADP is read 7 frames after the peak; at ramp firing rates a
% following spike lands inside that window often enough to dominate the
% average, and what would be measured is the mean interspike interval.
nFrLim = 20;
spk_t_blue_sta = isolatedBlueSpikes(stacks.spk_t_blue, nFrLim);

nBlue = cellfun(@numel, stacks.spk_t_blue);
nIso = cellfun(@numel, spk_t_blue_sta);
fprintf('\nIsolation (> %d frames, %.1f ms, on both sides)\n', nFrLim, nFrLim*dt*1e3);
fprintf('  blue-evoked spikes/cell: median %g\n', median(nBlue));
fprintf('  surviving isolation:     median %g (%.1f%% kept)\n', ...
    median(nIso), 100*median(nIso./nBlue));

%% ------------------------------------------------------------ the average
nBack = 140;
nFront = 140;
selfBlueStaFN_ADP = selfTriggeredAverage(data, spk_t_blue_sta, [nBack nFront]);

% The windows, in the source's frame offsets and in milliseconds.
adpWindow = nBack + 1 + 7 + (-5:5);
baseWindow = nBack + 1 - 11 + (-5:5);
fprintf('\nADP = mean(%+.1f to %+.1f ms) - mean(%+.1f to %+.1f ms) relative to peak\n', ...
    (adpWindow(1)-nBack-1)*dt*1e3, (adpWindow(end)-nBack-1)*dt*1e3, ...
    (baseWindow(1)-nBack-1)*dt*1e3, (baseWindow(end)-nBack-1)*dt*1e3);

adp_blue = mean(selfBlueStaFN_ADP(:, adpWindow), 2) - ...
           mean(selfBlueStaFN_ADP(:, baseWindow), 2);

%% -------------------------------------------- against the deposited STA
% The deposit carries the whole waveform, not just the scalar, so the STA can
% be checked sample by sample. Agreement is at float32 storage precision for
% every cell whose spike set the missing motion correction did not change.
staDiff = max(abs(selfBlueStaFN_ADP - data.deposited.selfBlueStaFN_ADP), [], 2);
fprintf('\nAgainst /processing/.../blue_evoked_spike_triggered_average\n');
fprintf('  %.1f%% of cells match to 1e-5 (float32 storage precision)\n', ...
    100*mean(staDiff < 1e-5));
fprintf('  the rest differ by up to %.2e, which is the motion-correction gap\n', ...
    max(staDiff));
fprintf('  cells affected: %s\n', mat2str(find(staDiff >= 1e-5)'));

fprintf('\nAgainst /units/adp_blue_norm\n');
fprintf('  max|diff| %.3e   r %.6f   %.1f%% within 1e-6\n', ...
    max(abs(adp_blue - data.deposited.adp_blue)), ...
    corr(adp_blue, data.deposited.adp_blue), ...
    100*mean(abs(adp_blue - data.deposited.adp_blue) < 1e-6));

%% ============================================ membrane time constant
%
% Triggers: every frame after which the blue steps down. That is the four
% pulse offsets plus the ramp offset, five per stimulated chunk. The source
% finds them with find(diff(blue_all(ii,:))<0) over the whole session row;
% here they come from the protocol and the stim map, which is the same set.
t_blue_off = cell(data.nCells, 1);
for ii = 1:data.nCells
    stimChunks = find(data.stimChunk(ii, :));
    t_blue_off{ii} = reshape((stimChunks(:)-1)*data.nFrame + protocol.offFrames, 1, []);
end

tBack = 140;
tFront = 390;
blueOffStaF_n = selfTriggeredAverage(data, t_blue_off, [tBack tFront]);

% Note what is being pooled: the four pulse offsets and the ramp offset are
% averaged together, so the estimate spans a range of preceding
% depolarisations rather than being a small-signal constant. That is the
% source's definition and it is kept.
fprintf('\nBlue-off triggers: %d per stimulated chunk, median %g per cell\n', ...
    numel(protocol.offFrames), median(cellfun(@numel, t_blue_off)));

mem = fitMembraneTimeConstant(blueOffStaF_n, tBack, dt);
membraneC = mem.membraneC;

fprintf('\nAgainst /units/membrane_time_constant_ms\n');
fprintf('  max|diff| %.3e   r %.6f\n', ...
    max(abs(membraneC - data.deposited.membraneC)), ...
    corr(membraneC, data.deposited.membraneC));

%% ------------------------------- the zero sentinels in M-YQ0201-27
% The deposited membraneC for M-YQ0201-27 contains seven exact zeros. Zero is
% not a measurable membrane time constant; it is what coefAll = zeros(nCells,3)
% leaves behind when the fit loop skips a cell. The current source code has
% only one skip path, isnan(std(blueOffStaF_n(ii,:))), and no cell in either
% deposited asset has a single non-finite sample, so that path cannot fire.
% The fits succeed here, several with adjusted R^2 above 0.98.
%
% Two of the seven pass QC, so they would enter a population histogram as
% 0 ms. See ISSUE_membraneC_zero_sentinels.md.
depZero = data.deposited.membraneC == 0;
if any(depZero)
    fprintf('\nDeposited zero sentinels: %d cells (%d pass QC)\n', ...
        sum(depZero), sum(depZero & data.idxUse));
    fprintf('  cell   thisRun   adjR^2   QC\n');
    for ii = find(depZero)'
        fprintf('  %4d  %8.2f   %.3f    %d\n', ii, membraneC(ii), mem.rSqMem(ii), data.idxUse(ii));
    end
    agree = ~depZero;
    fprintf('  excluding them: max|diff| %.3e   r %.6f\n', ...
        max(abs(membraneC(agree) - data.deposited.membraneC(agree))), ...
        corr(membraneC(agree), data.deposited.membraneC(agree)));
end

%% ----------------------------------------------------------------- plots
plotIntrinsicPropertyDiagnostics('waveforms', data, protocol, stacks, ...
    struct('staADP', selfBlueStaFN_ADP, 'staOff', blueOffStaF_n, ...
           'adp_blue', adp_blue, 'membraneC', membraneC, 'mem', mem, ...
           'nBack', nBack, 'tBack', tBack));

fprintf('\nAll five properties are now measured. See TUTORIAL2_VALIDATION.md.\n');
