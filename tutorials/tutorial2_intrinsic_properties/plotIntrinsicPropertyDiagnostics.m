function plotIntrinsicPropertyDiagnostics(stage, data, protocol, stacks, extra)
%PLOTINTRINSICPROPERTYDIAGNOSTICS Figures for tutorial 2.
%
%   PLOTINTRINSICPROPERTYDIAGNOSTICS(stage, data, protocol, stacks, extra)
%   draws the diagnostics for one stage. stage is 'protocol', 'rates',
%   'rheobase' or 'waveforms'.
%
%   This function computes nothing that affects a result. Everything it draws
%   is passed in already measured, so deleting the call cannot change a
%   number. That separation is a convention across these tutorials.

narginchk(4, 5);
if nargin < 5; extra = struct(); end
dt = data.dt;
tMs = @(f) (1:f) * dt * 1e3;

switch validatestring(stage, {'protocol', 'rates', 'rheobase', 'waveforms'})

%% ====================================================== protocol
case 'protocol'
    figure('Name', 'Tutorial 2a - stimulus', 'Color', 'w', ...
        'Position', [80 80 1000 720]);

    subplot(3, 2, [1 2]);
    plot(tMs(protocol.nFrame), data.blue, 'Color', [0.1 0.45 0.9], 'LineWidth', 1);
    hold on;
    for p = 1:numel(protocol.pulseAmp)
        xline(protocol.pulseOnset(p)*dt*1e3, ':', 'Color', [0.5 0.5 0.5]);
    end
    xline(protocol.rampOn*dt*1e3, 'r--', 'ramp on');
    xline(protocol.rampOff*dt*1e3, 'r--', 'ramp off');
    xlabel('time in chunk (ms)'); ylabel('blue (source units)');
    title(sprintf('one acquisition chunk: %d pulses then a %.1f s ramp', ...
        numel(protocol.pulseAmp), (protocol.rampOff-protocol.rampFirst+1)*dt));
    axis tight; box off;

    % Which cells are lit when. The sparseness is the point: each chunk lights
    % about ten cells, and each cell gets about ten chunks.
    subplot(3, 2, 3);
    imagesc(data.stimChunk);
    colormap(gca, [1 1 1; 0.1 0.45 0.9]);
    xlabel('acquisition chunk'); ylabel('cell');
    title('stimulation map');

    subplot(3, 2, 4);
    histogram(sum(data.stimChunk, 2), 'FaceColor', [0.1 0.45 0.9]);
    xlabel('repeats per cell'); ylabel('cells');
    title(sprintf('median %g repeats', median(sum(data.stimChunk, 2))));
    box off;

    % Population average response, the sanity check that the stack is aligned.
    subplot(3, 2, [5 6]);
    yyaxis left;
    plot(tMs(protocol.nFrame), mean(stacks.fr_mat_stim_avg, 1) / dt, 'LineWidth', 1);
    ylabel('population rate (Hz)');
    yyaxis right;
    plot(tMs(protocol.nFrame), data.blue, 'Color', [0.6 0.6 0.6]);
    ylabel('blue');
    xlabel('time in chunk (ms)');
    title('trial-averaged population rate, all cells');
    axis tight; box off;

%% ========================================================= rates
case 'rates'
    figure('Name', 'Tutorial 2b - firing-rate properties', 'Color', 'w', ...
        'Position', [100 100 1000 660]);

    jj = extra.adpPulse;
    on = protocol.pulseOnset(jj):protocol.pulseOffset(jj);

    subplot(2, 2, 1);
    plot((on - protocol.pulseOnset(jj)) * dt * 1e3, ...
        mean(stacks.fr_mat_stim_avg(data.idxUse, on), 1) / dt, 'k', 'LineWidth', 1);
    hold on;
    n1 = round(50e-3/dt); n2 = round(250e-3/dt);
    xregion = @(a, b, c) patch([a b b a], [0 0 1e4 1e4], c, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    xregion(0, n1*dt*1e3, [0.9 0.3 0.2]);
    xregion((numel(on)-n2)*dt*1e3, numel(on)*dt*1e3, [0.2 0.4 0.9]);
    ylim([0 max(mean(stacks.fr_mat_stim_avg(data.idxUse, on), 1)/dt)*1.2]);
    xlabel('time in pulse (ms)'); ylabel('rate (Hz)');
    title(sprintf('pulse %d (amp %.3f): early vs late window', jj, protocol.pulseAmp(jj)));
    box off;

    subplot(2, 2, 2);
    plot(extra.frPulseEnd(data.idxUse, jj)/dt, extra.frPulseStart(data.idxUse, jj)/dt, ...
        '.', 'Color', [0.2 0.2 0.2]);
    hold on; lim = [0 max(extra.frPulseStart(data.idxUse, jj)/dt)*1.05];
    plot(lim, lim, 'k--'); xlim(lim); ylim(lim);
    xlabel('late rate (Hz)'); ylabel('early rate (Hz)');
    title('adaptation: points above the line adapt'); box off;

    subplot(2, 2, 3);
    edges = linspace(min(extra.frMax), max(extra.frMax), 40);
    histogram(extra.frMax(data.idxUse & data.gfpIdx), edges, 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.1 0.5 0.2], 'LineWidth', 1.2); hold on;
    histogram(extra.frMax(data.idxUse & ~data.gfpIdx), edges, 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.6 0.2 0.6], 'LineWidth', 1.2);
    xlabel('frMax (Hz)'); ylabel('cells'); legend({'NPY+', 'NPY-'}, 'Box', 'off');
    title('maximum firing rate'); box off;

    subplot(2, 2, 4);
    plot(extra.frMax(data.idxUse), data.deposited.frMax(data.idxUse), '.', ...
        'Color', [0.2 0.2 0.2]);
    hold on; lim = [min(extra.frMax) max(extra.frMax)];
    plot(lim, lim, 'r--');
    xlabel('this run (Hz)'); ylabel('deposited (Hz)');
    title('frMax reproduces exactly'); axis square; box off;

%% ====================================================== rheobase
case 'rheobase'
    figure('Name', 'Tutorial 2c - rheobase', 'Color', 'w', ...
        'Position', [120 120 1000 660]);

    rampIdx = protocol.rampOn:protocol.rampOff;
    x = linspace(0, 1, numel(rampIdx));

    % One cell, worked through: the averaged subthreshold, the cubic, and the
    % point the median spike time picks out.
    probe = find(data.idxUse & isfinite(extra.vRheobase), 1);
    subplot(2, 2, 1);
    plot(x, stacks.traces_stim_avg(probe, rampIdx), 'Color', [0.6 0.6 0.6]); hold on;
    plot(x, polyval(extra.vFit{probe}, x), 'k', 'LineWidth', 1.5);
    plot(extra.fraction(probe), extra.vRheobase(probe), 'r.', 'MarkerSize', 28);
    xlabel('position along ramp'); ylabel('subthreshold V (spike heights)');
    title(sprintf('cell %d: cubic fit and rheobase point', probe)); box off;

    subplot(2, 2, 2);
    histogram(extra.fraction(data.idxUse), 40, 'FaceColor', [0.1 0.45 0.9]);
    xlabel('median time to Nth spike / ramp length'); ylabel('cells');
    title('where along the ramp threshold is reached'); box off;

    subplot(2, 2, 3);
    edges = linspace(min(extra.vRheobase), max(extra.vRheobase), 40);
    histogram(extra.vRheobase(data.idxUse & data.gfpIdx), edges, ...
        'DisplayStyle', 'stairs', 'EdgeColor', [0.1 0.5 0.2], 'LineWidth', 1.2); hold on;
    histogram(extra.vRheobase(data.idxUse & ~data.gfpIdx), edges, ...
        'DisplayStyle', 'stairs', 'EdgeColor', [0.6 0.2 0.6], 'LineWidth', 1.2);
    xlabel('vRheobase (spike heights)'); ylabel('cells');
    legend({'NPY+', 'NPY-'}, 'Box', 'off'); title('optical rheobase'); box off;

    subplot(2, 2, 4);
    plot(extra.vRheobase(data.idxUse), data.deposited.vRheobase(data.idxUse), '.', ...
        'Color', [0.2 0.2 0.2]); hold on;
    lim = [min(extra.vRheobase) max(extra.vRheobase)];
    plot(lim, lim, 'r--');
    xlabel('this run'); ylabel('deposited');
    title('vRheobase reproduces exactly'); axis square; box off;

%% ===================================================== waveforms
case 'waveforms'
    figure('Name', 'Tutorial 2d - ADP and membrane constant', 'Color', 'w', ...
        'Position', [140 140 1040 720]);

    nBack = extra.nBack;
    tau = (-nBack:nBack) * dt * 1e3;
    qc = data.idxUse;

    subplot(2, 3, 1);
    plot(tau, mean(extra.staADP(qc & data.gfpIdx, :), 1), 'Color', [0.1 0.5 0.2], ...
        'LineWidth', 1.4); hold on;
    plot(tau, mean(extra.staADP(qc & ~data.gfpIdx, :), 1), 'Color', [0.6 0.2 0.6], ...
        'LineWidth', 1.4);
    xline(7*dt*1e3, 'k:'); xline(-11*dt*1e3, 'k:');
    xlim([-30 40]); xlabel('time from spike (ms)'); ylabel('V (spike heights)');
    legend({'NPY+', 'NPY-'}, 'Box', 'off');
    title('blue-evoked spike average'); box off;

    subplot(2, 3, 2);
    edges = linspace(min(extra.adp_blue), max(extra.adp_blue), 40);
    histogram(extra.adp_blue(qc & data.gfpIdx), edges, 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.1 0.5 0.2], 'LineWidth', 1.2); hold on;
    histogram(extra.adp_blue(qc & ~data.gfpIdx), edges, 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.6 0.2 0.6], 'LineWidth', 1.2);
    xlabel('ADP (spike heights)'); ylabel('cells'); title('afterdepolarisation'); box off;

    subplot(2, 3, 3);
    staDiff = max(abs(extra.staADP - data.deposited.selfBlueStaFN_ADP), [], 2);
    semilogy(1:data.nCells, max(staDiff, 1e-9), '.', 'Color', [0.2 0.2 0.2]); hold on;
    yline(1e-5, 'r--', 'float32 precision');
    xlabel('cell'); ylabel('max |diff| from deposited STA');
    title(sprintf('%.1f%% bit-exact; rest is the motion-correction gap', ...
        100*mean(staDiff < 1e-5)));
    box off;

    tBack = extra.tBack;
    tOff = (-tBack:(size(extra.staOff, 2)-tBack-1)) * dt * 1e3;
    subplot(2, 3, 4);
    plot(tOff, mean(extra.staOff(qc & data.gfpIdx, :), 1), 'Color', [0.1 0.5 0.2], ...
        'LineWidth', 1.4); hold on;
    plot(tOff, mean(extra.staOff(qc & ~data.gfpIdx, :), 1), 'Color', [0.6 0.2 0.6], ...
        'LineWidth', 1.4);
    xline(0, 'k:'); xline(40*dt*1e3, 'r:', 'fit window end');
    xlabel('time from blue off (ms)'); ylabel('V (spike heights)');
    title('blue-off step response'); box off;

    subplot(2, 3, 5);
    good = ~extra.mem.skipped & ~extra.mem.atBound & qc;
    edges = linspace(0, max(extra.membraneC(good))*1.05, 40);
    histogram(extra.membraneC(good & data.gfpIdx), edges, 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.1 0.5 0.2], 'LineWidth', 1.2); hold on;
    histogram(extra.membraneC(good & ~data.gfpIdx), edges, 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.6 0.2 0.6], 'LineWidth', 1.2);
    xlabel('tau (ms)'); ylabel('cells');
    title(sprintf('membrane constant (%d railed/skipped excluded)', ...
        sum((extra.mem.atBound | extra.mem.skipped) & qc)));
    box off;

    subplot(2, 3, 6);
    depZero = data.deposited.membraneC == 0;
    plot(extra.membraneC(~depZero), data.deposited.membraneC(~depZero), '.', ...
        'Color', [0.2 0.2 0.2]); hold on;
    plot(extra.membraneC(depZero), data.deposited.membraneC(depZero), 'ro', ...
        'MarkerSize', 7, 'LineWidth', 1.2);
    lim = [0 max(extra.membraneC)*1.05];
    plot(lim, lim, 'r--'); xlim(lim); ylim(lim);
    xlabel('this run (ms)'); ylabel('deposited (ms)');
    if any(depZero)
        title(sprintf('%d deposited zero sentinels (red)', sum(depZero)));
    else
        title('membraneC reproduces');
    end
    axis square; box off;
end
end
