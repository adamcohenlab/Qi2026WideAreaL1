function session = loadSparsePulseHadForLme(nwbFile, varargin)
%LOADSPARSEPULSEHADFORLME Load the pairwise IPSP metrics of one sparsePulseHad session.
%
%   session = LOADSPARSEPULSEHADFORLME(nwbFile) reads the static pairwise
%   table /analysis/pairwise_connectivity_metrics and the quality-control flag
%   from /units, and returns the cell-by-cell matrices the Tutorial 5 mixed
%   models need. Field names follow the source pipeline (pMat, ipspAmpTestN2,
%   ipspDecay, seTau, distMat, dx) so the tutorial reads like the original
%   analysis.
%
%   session = LOADSPARSEPULSEHADFORLME(nwbFile, 'Name', value, ...) accepts:
%
%     'MinDistanceUm'  connections between cells closer than this are dropped
%                      from pMat (default 60). See DISTANCE CUT below.
%     'Verbose'        true (default) prints a short load report.
%
%   Only small datasets are read (102,400-element columns), so loading takes
%   about a second. No voltage data is touched.
%
%   FIELDS
%     pMat           nCells x nCells logical, row = presynaptic, column =
%                    postsynaptic. Detected inhibitory connection, after the
%                    distance cut.
%     ipspAmpTestN2  IPSP amplitude, in postsynaptic spike heights per
%                    presynaptic spike (normalized units).
%     ipspDecay      IPSP decay time constant. The deposit does not record its
%                    time unit (the median is about 0.047, i.e. 47 ms if in
%                    seconds). The models work on log(decay), so the unit only
%                    shifts the intercept and has no effect on any variance.
%     seTau          standard error of ipspDecay, same unit. NaN where no
%                    decay fit was attempted.
%     distMat        cell-to-cell distance in pixels.
%     dx             micrometres per pixel.
%     unitId         nCells x 1, zero-based NWB unit id of each kept cell.
%     nCells         number of cells that passed quality control.
%     audit          counts: cells before and after QC, connections removed
%                    by the distance cut.
%
%   QUALITY CONTROL
%   Only cells with /units/passed_qc = 1 are kept (290 of 320 in M-YQ0201-27,
%   286 of 320 in M-YQ0201-29). This is the source idxUse, and it is the cell
%   set the published analysis used.
%
%   DISTANCE CUT
%   The published analysis file fig3_pre_post.mat was exported by
%   sup_fig3_prepareData_connMatStatModel.m, which applies
%   pMat(~distMask) = false at line 355. Measured on both sessions, that mask
%   removes exactly the connections between cells closer than 60 um: 26 in
%   M-YQ0201-27 and 12 in M-YQ0201-29, and every connection it keeps is 60 um
%   or more apart. Very close pairs are the ones most exposed to shared
%   blue-light crosstalk. The deposited synaptic_connection column is the
%   uncut detection result, so the cut is re-applied here.
%
%   PRE/POST ORIENTATION
%   Values are placed by the table's presynaptic_unit_id and
%   postsynaptic_unit_id columns. NWB files written before 2026-09-26 carry
%   those two columns swapped, which would silently exchange pre and post
%   and therefore exchange the two variance components this tutorial
%   estimates. The converter always serialised the metrics with the
%   presynaptic index varying fastest; the defective files label them with
%   the postsynaptic index varying fastest instead. That layout is refused.
%
%   See also TUTORIAL5A_AMPLITUDEVARIANCE, TUTORIAL5B_DECAYVARIANCE.

opt = struct('MinDistanceUm', 60, 'Verbose', true);
for k = 1:2:numel(varargin)
    name = validatestring(varargin{k}, fieldnames(opt));
    opt.(name) = varargin{k+1};
end
assert(exist(nwbFile, 'file') == 2, 'Tutorial5:MissingFile', ...
    'NWB file not found: %s', nwbFile);

tablePath = '/analysis/pairwise_connectivity_metrics/';
readColumn = @(name) double(h5read(nwbFile, [tablePath name]));

%% Units: quality control and pixel size
passedQc = logical(h5read(nwbFile, '/units/passed_qc'));
passedQc = passedQc(:);
nAll = numel(passedQc);

% x_um is written as x_pixel * dx, so their ratio recovers dx without parsing
% the free-text /general/notes.
xPixel = double(h5read(nwbFile, '/units/x_pixel'));
xUm = double(h5read(nwbFile, '/units/x_um'));
ratio = xUm(xPixel ~= 0) ./ xPixel(xPixel ~= 0);
dx = median(ratio);
assert(max(abs(ratio - dx)) < 1e-9 * dx, 'Tutorial5:PixelSize', ...
    'x_um / x_pixel is not constant across units.');

%% Pairwise table: place every row by its explicit pre/post unit ids
pre = readColumn('presynaptic_unit_id');
post = readColumn('postsynaptic_unit_id');
assert(numel(pre) == nAll^2 && numel(post) == nAll^2, 'Tutorial5:TableSize', ...
    'The pairwise table should have %d rows, one per directed pair.', nAll^2);
assert(all(pre >= 0 & pre < nAll & post >= 0 & post < nAll), ...
    'Tutorial5:UnitId', 'Pairwise unit ids fall outside the units table.');
linearIndex = sub2ind([nAll nAll], pre + 1, post + 1);
assert(numel(unique(linearIndex)) == nAll^2, 'Tutorial5:DuplicatePair', ...
    'The pairwise table does not hold exactly one row per directed pair.');

% Refuse the pre-2026-09-26 layout, in which the id columns are transposed
% relative to the stored values. In a correct file the presynaptic id varies
% fastest down the table, matching MATLAB's column-major metric(:).
if pre(1) == pre(2) && post(2) == post(1) + 1
    error('Tutorial5:SwappedPrePostIds', ...
        ['%s\nhas the presynaptic_unit_id / postsynaptic_unit_id columns swapped ' ...
        '(the layout written before 2026-09-26). Reading it would exchange pre ' ...
        'and post. Use a reconverted file.'], nwbFile);
end

toMatrix = @(values) placeValues(values, linearIndex, nAll);
pMatAll = toMatrix(readColumn('synaptic_connection')) > 0;
ampAll = toMatrix(readColumn('ipsp_amplitude_norm_per_presynaptic_spike'));
decayAll = toMatrix(readColumn('ipsp_decay_source_units'));
seAll = toMatrix(readColumn('ipsp_decay_standard_error_source_units'));
distAll = toMatrix(readColumn('distance_pixels'));

%% Keep QC cells, then apply the distance cut to the connection matrix
keep = find(passedQc);
session = struct();
session.nwbFile = nwbFile;
session.dx = dx;
session.unitId = keep - 1;
session.nCells = numel(keep);
session.distMat = distAll(keep, keep);
session.ipspAmpTestN2 = ampAll(keep, keep);
session.ipspDecay = decayAll(keep, keep);
session.seTau = seAll(keep, keep);

pMat = pMatAll(keep, keep);
tooClose = session.distMat * dx < opt.MinDistanceUm;
session.pMat = pMat & ~tooClose;

session.audit = struct( ...
    'nCellsInFile', nAll, ...
    'nCellsPassedQc', session.nCells, ...
    'nConnectionsBeforeCut', nnz(pMat), ...
    'nConnectionsRemovedByCut', nnz(pMat & tooClose), ...
    'nConnections', nnz(session.pMat), ...
    'minDistanceUm', opt.MinDistanceUm);

if opt.Verbose
    [~, name] = fileparts(nwbFile);
    fprintf(['%s: %d of %d cells passed QC; %d connections, of which %d closer ' ...
        'than %g um were removed, leaving %d.\n'], name, session.nCells, nAll, ...
        session.audit.nConnectionsBeforeCut, session.audit.nConnectionsRemovedByCut, ...
        opt.MinDistanceUm, session.audit.nConnections);
end
end

function matrix = placeValues(values, linearIndex, n)
matrix = nan(n);
matrix(linearIndex) = values;
end
