%% Tutorial 3a -- how many spikes does a blue pulse evoke?
%
% Every published IPSP amplitude in this study is normalised "per presynaptic
% spike". That denominator is nSpikesPulse, and it is worth understanding before
% anything else, because it is not a spike count.
%
% The problem: during the 20 ms pulse the stimulated cell sits on a large
% optogenetic depolarisation and an optical crosstalk artifact. Individual
% spikes there cannot be detected reliably. The same fact shapes the crosstalk
% template in tutorial3c, where the rise constant is discarded for the same
% reason.
%
% The solution: remove the stereotyped pulse response, then read the firing rate
% off the *spectral width* of what remains. Spikes are brief events, so more of
% them broadens the residual's power spectrum. The centroid of that spectrum,
% times the pulse duration, is the expected count.
%
% This tutorial reproduces the deposited /units/estimated_spikes_per_blue_pulse
% from the NWB file alone, and checks it.
%
% Source: bluePulseSpikeCount.m, the block headed "%% Correct version".

clear; clc

%% 1. Configuration
cfg = struct;
cfg.nwbFile = 'path/to/M_YQ0201_27_sparsePulseHad.nwb';   % edit: your local copy
cfg.nCellsToRun = 24;        % Inf reproduces all 320; 24 is a quick check
cfg.plotCell = 3;            % a cell index to show diagnostics for, or [] for none

% Helper functions from the source pipeline, bundled in tutorials/helpers.
% These are used unchanged.
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helpers')));
for fn = {'get_sta_mat_raw','fitSelfBlueSta','SeeResiduals_vec'}
    assert(exist(fn{1},'file') == 2, 'Tutorial3a:MissingHelper', ...
        '%s is not on the MATLAB path.', fn{1});
end

% Window literals, exactly as in the source. With nBack = 140, tau = 0 sits at
% index 141, so tOn = 142 is the first frame of the pulse (+1 under convention
% A) and tOff = 158 is the last (+17, i.e. about 20 ms at 1.27 ms/frame).
prm = struct('tPrev', 1, 'tOn', 142, 'tOff', 158, ...
             'nBack', 140, 'nFront', 140, 'dt', 1.27e-3, 'tPulse', 20e-3);

%% 2. Load
% The full trace matrix is about 4 GB as single. Both per-cell datasets are
% chunked [nCells 10000], so reading one cell costs a whole-file pass anyway;
% loading once is both simpler and faster than reading cell by cell.
data = loadSparsePulseHadForIpsp(cfg.nwbFile);

nRun = min(cfg.nCellsToRun, data.nCells);
if isinf(cfg.nCellsToRun), nRun = data.nCells; end

%% 3. Estimate, cell by cell
fprintf('\nEstimating spikes per pulse for %d cells...\n', nRun);
estimate = nan(nRun,1);
nTrials = nan(nRun,1);
timer = tic;
for ii = 1:nRun
    out = estimateSpikesPerPulse(data.traces_all(ii,:), data.t_blue{ii}, prm);
    estimate(ii) = out.nSpikesPulse;
    nTrials(ii) = out.nTrials;
    if mod(ii, 10) == 0 || ii == nRun
        fprintf('  [%5.1f s] %d/%d\n', toc(timer), ii, nRun);
    end
end

%% 4. Check against the deposit
deposited = data.nSpikesPulse(1:nRun);
absErr = abs(estimate - deposited);
relErr = absErr ./ abs(deposited);

fprintf('\n--- reproduction of /units/estimated_spikes_per_blue_pulse ---\n');
fprintf('cells compared        : %d\n', nRun);
fprintf('pulses per cell       : median %d [%d, %d]\n', ...
    median(nTrials), min(nTrials), max(nTrials));
fprintf('deposited value       : median %.4f [%.4f, %.4f]\n', ...
    median(deposited), min(deposited), max(deposited));
fprintf('reproduced value      : median %.4f [%.4f, %.4f]\n', ...
    median(estimate), min(estimate), max(estimate));
fprintf('max absolute error    : %.3e\n', max(absErr));
fprintf('max relative error    : %.3e\n', max(relErr));
fprintf('correlation           : %.8f\n', corr(estimate, deposited));

tol = 1e-6;
if max(relErr) < tol
    fprintf('\nPASS: reproduced to better than %g relative.\n', tol);
    fprintf('Note this used normalized_voltage, while the source used raw\n');
    fprintf('traces_all_bc. The estimator is a ratio of spectra, so the\n');
    fprintf('per-cell spkHgtNoBlue factor cancels exactly.\n');
else
    fprintf('\nMISMATCH above %g. Check TUTORIAL3_VALIDATION.md before\n', tol);
    fprintf('proceeding; the two variants in bluePulseSpikeCount.m differ.\n');
end

%% 5. Optional per-cell diagnostics
if ~isempty(cfg.plotCell)
    ii = cfg.plotCell;
    out = estimateSpikesPerPulse(data.traces_all(ii,:), data.t_blue{ii}, prm);
    plotSpikeCountDiagnostics(out, ii, prm);
end

%% 6. What to take away
% nSpikesPulse is a rate estimate, not an event count, and it is the
% denominator of every normalised IPSP amplitude downstream. Two consequences
% worth carrying forward:
%
%  - It is a per-PREsynaptic-cell quantity. In tutorial3b and 3d the observed
%    pair waveform is divided by nSpikesPulse(pre), and so must every null
%    waveform be. Missing that is the single easiest way to invalidate the
%    p-values.
%
%  - It exists because spikes are unmeasurable during the pulse. Tutorial 3c
%    inherits the same constraint and handles it a different way.
