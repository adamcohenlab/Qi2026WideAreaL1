%% Tutorial 2c - optical rheobase
%
% Rheobase is the input a cell needs before it fires. In current clamp you
% ramp the injected current and read the current at the first spike. Here the
% ramp is optogenetic and the read-out is optical, which changes two things:
% the "input" is not in amperes, and the voltage is not in millivolts. Both
% problems are solved the same way, by expressing the answer as a voltage in
% spike-height units, read at the moment the cell reaches threshold.
%
% The method:
%
%   1. Take the time from ramp onset to the Nth blue-evoked spike, per repeat.
%   2. Take the median over that cell's repeats. The Nth spike rather than the
%      first, so a single early spike cannot set the answer.
%   3. Fit a cubic to the trial-averaged subthreshold voltage across the ramp.
%   4. Evaluate the cubic at the median time. That is vRheobase.
%
% This stage is where the subthreshold estimate matters, and where this
% tutorial deliberately does not use the subthreshold series in the deposit.
% See the comparison at the end.
%
% Runtime: about 60 s.

close all;

% Set nwbFile before running to point at your copy of a sparseOpto session,
% or edit the placeholder below. Both deposited sparseOpto sessions work.
if ~exist('nwbFile', 'var') || isempty(nwbFile)
    nwbFile = 'path/to/M_YQ0201_27_sparseOpto.nwb';
end

data = loadSparseOptoForEphys(nwbFile);
protocol = blueProtocol(data.blue, data.dt);

%% ------------------------------------------------- the subthreshold trace
% 'interpolate' cuts a [-3 +5] frame window around every spike and fills it by
% linear interpolation. That is the estimate the property was defined on.
%
% The deposit carries a different one, traces_all_sub_n, built with a
% median-filter substitution and then sgolayfilt(...,2,17). It is the better
% estimate for a single trace, but this stage averages ten repeats before
% reading anything, and that average already removes the noise the filter
% targets. Applying a 17-frame polynomial filter on top only distorts the slow
% ramp that the cubic is there to measure. The two are compared at the end.
stacks = stackStimTrials(data, 'Subthreshold', 'interpolate');

%% ---------------------------------------------------------- the estimate
rheo = estimateRheobase(stacks, protocol, 'nSpike', 10);

%% ----------------------------------------------- where the method strains
% Two failure modes are worth seeing rather than trusting.
%
% The leak: the source computes the time to the Nth spike over the whole
% session, not within the repeat, so a repeat in which the cell fired fewer
% than N times after ramp onset borrows its Nth spike from the next
% stimulation chunk thousands of frames later. The median over repeats
% normally absorbs it.
fprintf('\nLeaked repeats (Nth spike found past the end of the ramp)\n');
fprintf('  %d repeats over %d cells; %d of those cells pass QC\n', ...
    sum(rheo.leaked), sum(rheo.leaked > 0), sum(rheo.leaked > 0 & data.idxUse));

% The rail: if the median time still lands outside the ramp, the cubic is
% being extrapolated and the value is not a measurement.
fprintf('Cells whose median time falls outside the ramp: %d\n', sum(rheo.railed));

%% ---------------------------------------------- the spike criterion N
% The source contains two rheobase blocks. The first sweeps nSpkRh = 5:10 and
% plots the histograms; the second fixes nSpkRh = 10 and is the one whose
% output is saved. The sweep is reproduced here as a sensitivity check rather
% than as a second pipeline: if the answer moved a lot with N, the property
% would be measuring the firing rate instead of the threshold.
nSweep = 5:10;
vSweep = zeros(data.nCells, numel(nSweep));
for k = 1:numel(nSweep)
    r = estimateRheobase(stacks, protocol, 'nSpike', nSweep(k), 'Verbose', false);
    vSweep(:, k) = r.vRheobase;
end
qc = data.idxUse;
fprintf('\nSensitivity to the spike criterion N (QC cells only)\n');
for k = 1:numel(nSweep)
    fprintf('  N = %2d: median %.4f, r with N=10 = %.4f\n', nSweep(k), ...
        median(vSweep(qc, k)), corr(vSweep(qc, k), vSweep(qc, end)));
end

%% ------------------------------------------------------- against the deposit
fprintf('\nAgainst /units/optical_rheobase_norm (all %d cells)\n', data.nCells);
fprintf('  max|diff| %.3e   r %.6f   %.1f%% within 1e-6\n', ...
    max(abs(rheo.vRheobase - data.deposited.vRheobase)), ...
    corr(rheo.vRheobase, data.deposited.vRheobase), ...
    100*mean(abs(rheo.vRheobase - data.deposited.vRheobase) < 1e-6));

%% ------------------------------- what the deposited subthreshold would give
% Re-run the whole stage on traces_all_sub_n to measure what the filtering
% choice costs. This is the quantified version of the claim above.
stacksDep = stackStimTrials(data, 'Subthreshold', 'deposited', 'Verbose', false);
rheoDep = estimateRheobase(stacksDep, protocol, 'nSpike', 10, 'Verbose', false);
fprintf('\nSame property from the deposited (Savitzky-Golay) subthreshold\n');
fprintf('  max|diff| from published %.3e   r %.6f\n', ...
    max(abs(rheoDep.vRheobase - data.deposited.vRheobase)), ...
    corr(rheoDep.vRheobase, data.deposited.vRheobase));
fprintf('  interpolated vs filtered: r %.6f, median |diff| %.3e\n', ...
    corr(rheo.vRheobase, rheoDep.vRheobase), ...
    median(abs(rheo.vRheobase - rheoDep.vRheobase)));

%% --------------------------- chunk-local interpolation is session-global
% subthresholdByInterpolation runs per acquisition chunk, not over the whole
% session. Interpolation is local, so the two agree except within [-3 +5] of a
% chunk edge. Checked here on one cell rather than asserted.
% Interpolate the same chunk twice: once on its own, and once inside a
% three-chunk window so the middle chunk has real neighbours on both sides.
% If chunk-local interpolation were losing anything, the two would differ
% inside the stimulus window.
probe = find(data.idxUse, 1);
chunkIdx = find(data.stimChunk(probe, :), 1);
chunkIdx = min(max(chunkIdx, 2), data.nChunk - 1);
nFrame = data.nFrame;
spkAll = data.spk_t{probe};

% One cell's row from three consecutive chunks, in time order.
rows = cell(1, 3);
spkWide = [];
for j = 1:3
    k = chunkIdx - 2 + j;
    block = data.readChunk('voltage', k);
    rows{j} = block(probe, :);
    localSpk = spkAll(floor((spkAll-1)/nFrame)+1 == k) - (k-1)*nFrame;
    spkWide = [spkWide, localSpk + (j-1)*nFrame]; %#ok<AGROW>
end

subLocal = subthresholdByInterpolation(rows{2}, spkWide(spkWide > nFrame & ...
    spkWide <= 2*nFrame) - nFrame);
subWide = subthresholdByInterpolation([rows{1} rows{2} rows{3}], spkWide);
subWideMid = subWide(nFrame + (1:nFrame));

stimWindow = protocol.pulseOnset(1):protocol.rampOff;
fprintf('\nChunk-local versus three-chunk interpolation, cell %d, chunk %d\n', ...
    probe, chunkIdx);
fprintf('  max|diff| over the stimulus window %d..%d: %.3e\n', ...
    stimWindow(1), stimWindow(end), ...
    max(abs(double(subLocal(stimWindow)) - double(subWideMid(stimWindow)))));
fprintf('  max|diff| over the whole chunk (edges included): %.3e\n', ...
    max(abs(double(subLocal) - double(subWideMid))));
clear rows subWide;

%% ----------------------------------------------------------------- plots
plotIntrinsicPropertyDiagnostics('rheobase', data, protocol, stacks, rheo);

fprintf('\nNext: tutorial2d_adpAndMembraneConstant.m\n');
