function plotPairWaveformDiagnostics(preCellIndex, postIndices, tau, raw, afterFar, ...
    afterAnnulus, connectedMask)
%PLOTPAIRWAVEFORMDIAGNOSTICS Show what each baseline correction removes.
%
%   plotPairWaveformDiagnostics(preCellIndex, postIndices, tau, raw, ...
%       afterFar, afterAnnulus, connectedMask)
%
%   Three panels, all showing every nearby pair for one stimulated cell:
%
%     1. the raw pulse-triggered averages;
%     2. after subtracting the far-pulse baseline, which removes the
%        field-wide artifact;
%     3. after also subtracting the neighbour-annulus median, which removes
%        what is left of that cell's scattered light at mid range. Pairs the
%        deposit calls connected are drawn in red.
%
%   The useful thing to notice is that panel 3 still contains a large positive
%   transient. The corrections remove the field-wide and mid-range artifact,
%   not the local one, and nearby pairs are exactly where scattered light is
%   strongest -- so the residual artifact is largest precisely where the
%   synaptic signal lives. That is why tutorial3c models the remainder instead
%   of trying to subtract it.
%
%   raw, afterFar and afterAnnulus are each [nPost x nTau] for this cell.
%
%   Draws only.
%
%   See also TUTORIAL3B_BLUEPULSESTA.

% Drop the self pair. Cell i's response to its own pulse is an order of
% magnitude larger than any cross-pair -- it is essentially pure crosstalk plus
% the cell's own depolarisation -- so leaving it in sets the y-scale and makes
% the three panels look identical. It is never a detection candidate either:
% the detector excludes self pairs. Its shape is the subject of tutorial3c.
keep = postIndices(:)' ~= preCellIndex;
raw = raw(keep,:); afterFar = afterFar(keep,:); afterAnnulus = afterAnnulus(keep,:);
if ~isempty(connectedMask), connectedMask = connectedMask(keep); end
nPair = sum(keep);

figure('Name', sprintf('Tutorial 3b -- pulses on cell %d', preCellIndex), 'Color','w');
grey = [.6 .6 .6];
xl = [-.05 .15];
% one shared y-scale, so the panels are actually comparable
yl = [min(afterAnnulus(:)) max(raw(:))];
yl = yl + 0.08*diff(yl)*[-1 1];

subplot(1,3,1)
plot(tau, raw', 'Color', grey);
xlim(xl); ylim(yl); xlabel('time from pulse (s)'); ylabel('normalised \DeltaV')
title(sprintf('raw: %d nearby pairs (self excluded)', nPair))

subplot(1,3,2)
plot(tau, afterFar', 'Color', grey);
xlim(xl); ylim(yl); xlabel('time from pulse (s)')
title('minus far-pulse baseline')

subplot(1,3,3)
plot(tau, afterAnnulus', 'Color', grey); hold on
nConn = 0;
if ~isempty(connectedMask) && any(connectedMask)
    plot(tau, afterAnnulus(connectedMask,:)', 'r', 'LineWidth', 1.2);
    nConn = sum(connectedMask);
end
xlim(xl); ylim(yl); xlabel('time from pulse (s)')
title(sprintf('minus annulus median (%d connected, red)', nConn))
end
