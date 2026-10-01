function out = estimateRheobase(stacks, protocol, varargin)
%ESTIMATERHEOBASE Optical rheobase from the ramp response.
%
%   out = ESTIMATERHEOBASE(stacks, protocol) returns
%
%     vRheobase    nCells x 1, the trial-averaged subthreshold voltage at the
%                  moment the cell reaches its spike criterion during the ramp,
%                  in normalised (spike-height) units.
%     tRheobase    nCells x 1 cell array, frames from ramp onset to the Nth
%                  blue-evoked spike, one entry per repeat.
%     vFit         nCells x 1 cell array of cubic coefficients.
%     fraction     nCells x 1, median(tRheobase) as a fraction of ramp length.
%                  This is the abscissa the cubic is evaluated at.
%     leaked       nCells x 1, repeats whose Nth spike fell outside the ramp.
%     railed       nCells x 1 logical, fraction outside [0 1], where the cubic
%                  is extrapolating.
%
%   ESTIMATERHEOBASE(..., 'Name', value) accepts:
%
%     'nSpike'   spike criterion, default 10 (the source's nSpkRh).
%     'Verbose'  true (default).
%
%   THE METHOD
%   A slow optogenetic ramp drives the cell from rest through threshold. The
%   frame at which it first fires reliably is a current-clamp rheobase read
%   optically: take the time from ramp onset to the Nth evoked spike, take the
%   median over repeats, and read the trial-averaged subthreshold voltage there.
%   A cubic in normalised ramp position smooths the read, because the ramp is
%   5 s long and the average still carries noise at single-frame scale.
%
%   Two details of the source are preserved because the published number
%   depends on them. The abscissa is
%
%       median(tRheobase{ii}) / (tRampOff - tRampOn)
%
%   against a cubic fitted on x = linspace(0,1,tRampOff-tRampOn+1), so the
%   divisor is the ramp length in frames while the fit spans length+1 samples;
%   the resulting half-frame offset is part of the definition. And tRheobase
%   counts the Nth spike, not the first, which makes it robust to a single
%   early spike but means a cell that never reaches N spikes in one repeat
%   borrows from the next.
%
%   THE LEAK
%   The source computes tDiff = spk_t_blue{ii} - tRamp(iStim) over the WHOLE
%   session and keeps every positive difference, so when a repeat contains
%   fewer than N evoked spikes after ramp onset the Nth spike is found in a
%   later stimulation chunk, thousands of frames away. The median over repeats
%   usually absorbs it. This function reproduces the behaviour exactly, and
%   counts the occurrences in `leaked` and `railed` rather than hiding them.
%
%   See also STACKSTIMTRIALS, TUTORIAL2C_RHEOBASE.

opt = struct('nSpike', 10, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
nSpkRh = opt.nSpike;

assert(~isempty(stacks.traces_stim_avg), 'Tutorial2:NoSubthreshold', ...
    'stacks.traces_stim_avg is empty; stack with Subthreshold interpolate.');

nCells = size(stacks.fr_mat_stim_avg, 1);
tRampOn = protocol.rampOn;
tRampOff = protocol.rampOff;
rampLen = tRampOff - tRampOn;
nFrame = protocol.nFrame;

vRheobase = NaN(nCells, 1);
fraction = NaN(nCells, 1);
leaked = zeros(nCells, 1);
tRheobase = cell(nCells, 1);
vFit = cell(nCells, 1);

% Ramp onsets are the same local frame in every stimulated chunk, so they can
% be enumerated from the stim map instead of running findpeaks over the whole
% 2.4 million-frame blue row as the source does.
x = linspace(0, 1, rampLen + 1);

for ii = 1:nCells
    spkBlue = sort(stacks.spk_t_blue{ii});
    if isempty(spkBlue); continue; end

    % Every stimulated chunk is a repeat, including one in which the cell did
    % not fire. Taking the repeats from the spike list instead would silently
    % drop those and shift the median.
    chunksStimulated = find(stacks.stimChunk(ii, :));
    tRamp = (chunksStimulated(:)' - 1) * nFrame + tRampOn;

    tRheobase_i = zeros(1, numel(tRamp));
    for iStim = 1:numel(tRamp)
        tDiff = spkBlue - tRamp(iStim);
        tDiff = sort(tDiff(tDiff > 0));
        if isempty(tDiff)
            % The source would error here on tDiff(0). Only reachable for the
            % final repeat of a cell that stops firing; recorded, not fudged.
            tRheobase_i(iStim) = NaN;
            leaked(ii) = leaked(ii) + 1;
            continue;
        end
        tRheobase_i(iStim) = tDiff(min(nSpkRh, numel(tDiff)));
        if tRheobase_i(iStim) > rampLen
            leaked(ii) = leaked(ii) + 1;
        end
    end
    tRheobase{ii} = tRheobase_i;

    trace = stacks.traces_stim_avg(ii, tRampOn:tRampOff);
    if any(~isfinite(trace)); continue; end
    p = polyfit(x, trace, 3);
    vFit{ii} = p;

    fraction(ii) = median(tRheobase_i, 'omitnan') / rampLen;
    vRheobase(ii) = polyval(p, fraction(ii));
end

out = struct();
out.vRheobase = vRheobase;
out.tRheobase = tRheobase;
out.vFit = vFit;
out.fraction = fraction;
out.leaked = leaked;
out.railed = fraction < 0 | fraction > 1;
out.nSpike = nSpkRh;

if opt.Verbose
    ok = isfinite(vRheobase);
    fprintf('Rheobase (N = %d spikes) for %d cells\n', nSpkRh, sum(ok));
    fprintf('  vRheobase range [%.3f %.3f]\n', min(vRheobase(ok)), max(vRheobase(ok)));
    fprintf('  repeats whose Nth spike fell past the ramp: %d over %d cells\n', ...
        sum(leaked), sum(leaked > 0));
    fprintf('  cells extrapolating outside the ramp: %d\n', sum(out.railed));
end
end
