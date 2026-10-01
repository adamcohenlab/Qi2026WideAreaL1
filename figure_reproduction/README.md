# Figure reproduction

These scripts rebuild the data panels of manuscript Figures 2 to 6 from the
public NWB files. Each figure has one script that you run, and one loader
that reads the NWB file for it.

The published figures were assembled and styled outside MATLAB. These scripts
reproduce the **analysis and the plotted data**, not the final layout, fonts,
colours, scale bars or schematics. Expect raw MATLAB figures that carry the
same content as the published panels.

## Before you start

1. Add the repository to the MATLAB path:

   ```matlab
   addpath(genpath('path/to/this/repository'))
   ```

2. Download the NWB files for the figure you want (see the table below). The
   Data section of the [main README](../README.md) says where to get them.
3. Check that you have the MATLAB toolboxes listed for that figure.

## Which files each figure needs

| Figure | Script | NWB files | Toolboxes | Time |
|---|---|---|---|---|
| 2 | `reproduce_figure2_from_nwb` | `M_YQ0201_29_sparseOpto.nwb` | Statistics and Machine Learning; Curve Fitting; UMAP (optional) | not timed |
| 3 | `reproduce_figure3_from_nwb` | `M_YQ0201_27_sparsePulseHad.nwb`, `M_YQ0201_29_sparsePulseHad.nwb` | Statistics and Machine Learning | about 40 s |
| 4 | `reproduce_figure4_from_nwb` | `M_YQ0201_27_behavior.nwb`, `M_YQ0201_29_behavior.nwb` | Statistics and Machine Learning | about 3 min |
| 5 | `reproduce_figure5_from_nwb` | `M_YQ0201_27_sparseOpto.nwb`, `M_YQ0201_27_sparsePulseHad.nwb` | Statistics and Machine Learning | 1 to 2 min |
| 6 | `reproduce_figure6_from_nwb` | `M_YQ0201_29_behavior.nwb`, `M_YQ0201_27_behavior.nwb`, `M_YQ0201_29_sparseOpto.nwb`, `M_YQ0201_27_sparsePulseHad.nwb` | Signal Processing; Statistics and Machine Learning; Curve Fitting; UMAP | 3 to 4 min |

The file names are those used during analysis. The names on DANDI may
differ. `M_YQ0201_27` and `M_YQ0201_29` are the two mice. The three session
types are:

- **sparseOpto**: intrinsic-excitability mapping. A few cells at a time get
  four 500 ms light steps and a 5 s ramp.
- **sparsePulseHad**: connectivity mapping. Cells get short blue pulses, and
  every fourth recording chunk uses a Hadamard stimulation pattern for
  gap-junction mapping.
- **behavior**: spontaneous activity with no stimulation, with running and
  whisking recorded.

UMAP is the [UMAP package for MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)
from File Exchange. Install it and put it on the path. Figure 2 still runs
without it; it only leaves out the UMAP scatter. Figure 6 needs it for panels
E, F and G.

Times were measured on a workstation with the files on a local drive. Reading
from a network drive is much slower.

## Running a script

Every argument is optional, but the default file paths are placeholders
(`path/to/...`). So in practice you always pass your own file paths.

```matlab
% Figure 2: one sparseOpto file
results = reproduce_figure2_from_nwb('D:/data/M_YQ0201_29_sparseOpto.nwb');

% Figure 3: the two sparsePulseHad files, mouse 27 first
results = reproduce_figure3_from_nwb( ...
    'D:/data/M_YQ0201_27_sparsePulseHad.nwb', ...
    'D:/data/M_YQ0201_29_sparsePulseHad.nwb');

% Figure 4: the two behavior files, mouse 27 first
results = reproduce_figure4_from_nwb( ...
    'D:/data/M_YQ0201_27_behavior.nwb', ...
    'D:/data/M_YQ0201_29_behavior.nwb');

% Figure 5: the mouse 27 sparseOpto and sparsePulseHad files
results = reproduce_figure5_from_nwb( ...
    'D:/data/M_YQ0201_27_sparseOpto.nwb', ...
    'D:/data/M_YQ0201_27_sparsePulseHad.nwb');

% Figure 6: a list of behavior files, then sparseOpto, then sparsePulseHad
results = reproduce_figure6_from_nwb( ...
    {'D:/data/M_YQ0201_29_behavior.nwb'; 'D:/data/M_YQ0201_27_behavior.nwb'}, ...
    'D:/data/M_YQ0201_29_sparseOpto.nwb', ...
    'D:/data/M_YQ0201_27_sparsePulseHad.nwb');
```

The full argument lists are:

```
reproduce_figure2_from_nwb(nwbFile, outputDirectory)
reproduce_figure3_from_nwb(nwbFile27, nwbFile29, outputDirectory, panels)
reproduce_figure4_from_nwb(nwbFile27, nwbFile29, outputDirectory, panels)
reproduce_figure5_from_nwb(optoNwbFile, connNwbFile, outputDirectory, panels)
reproduce_figure6_from_nwb(behaviorNwbFiles, optoNwbFile, connNwbFile, outputDirectory, panels)
```

For Figure 6, the **first** behavior file is the session shown in panels E, F
and H to K. Panel G pools every behavior file in the list.

### Running only some panels

For Figures 3 to 6 you can ask for a subset of panels. Give one letter, a
comma-separated list, or a cell array. Pass `[]` for any argument you want to
leave at its default. Only the files those panels need are opened.

```matlab
reproduce_figure3_from_nwb(f27, f29, [], 'F')
reproduce_figure5_from_nwb([], connFile, [], 'D,E,F,G')   % never opens the sparseOpto file
reproduce_figure6_from_nwb([], [], connFile, [], 'D')     % needs only the sparsePulseHad file
reproduce_figure4_from_nwb(f27, f29, [], 'STATS')         % the numbers quoted in the text
```

Figure 4 has an extra selector, `'STATS'`. It recomputes the numbers the
manuscript text quotes for Figure 4: length constants, short-range
correlations, correlation widths and the cross-STA depolarization.

## What you get

Each script:

- draws each panel in its own figure window and **leaves the windows open**
  so you can inspect them;
- saves each panel as a PNG in `outputDirectory`;
- saves a reference file, `figureN_nwb_reference_results.mat`, holding the
  numbers behind every panel;
- returns the same numbers as a `results` structure.

`outputDirectory` defaults to a `figureN_nwb_reproduction` folder next to the
script. A run of only some panels writes a reference file named after those
panels, so it does not overwrite the full one.

If the output folder is inside a Dropbox or OneDrive folder, the sync client
can briefly lock new files and make a save fail. Point `outputDirectory`
somewhere outside the synced folder if that happens.

## Checking your run

These numbers should come out exactly. They are the sample sizes and
statistics printed on the published panels or in the captions.

| Figure | Quantity | Expected |
|---|---|---|
| 2 | cells per cluster (Ward clustering of 320 cells) | 94, 116, 110 |
| 3A | trials in the synaptic example | 283 |
| 3C | spikes in the gap-junction example | 7215 |
| 3J | presynaptic identity vs IPSP amplitude | R = 0.61, p = 1.052e-51, n = 561 |
| 3K | postsynaptic identity vs IPSP amplitude | R = 0.26, p = 4.341e-09, n = 561 |
| 3M | presynaptic identity vs IPSP decay | R = 0.07, p = 5.617e-02, n = 706 |
| 3N | postsynaptic identity vs IPSP decay | R = 0.33, p = 4.714e-18, n = 706 |
| 4G | cells, two mice | 505 (247 + 258) |
| 4H | neighbour distances labelled on the two examples | 88, 115, 130, 183, 198 µm and 134, 136, 136, 177, 177 µm |
| 5A-C | cells in the quintiles | 280 |
| 5D-E | connected pairs | 800 |
| 5F-G | gap-junction pairs | 3553 |
| 6D | pairs, (−)→(−), (−)↔(+), (+)→(+) | 1420, 1818, 2224 |
| 6G | NPY status vs UMAP group | 157, 100 / 46, 202 (505 cells) |
| 6H | NPY(−) and NPY(+) cells | 141, 117 |
| 6I | oscillation-peak events | 5888 |
| 6J | whisking-onset events | 547 |

For Figure 3, n is the number of pairs that pass the selection. The number of
points actually plotted is a little smaller, because cells with no other
qualifying pair drop out of the leave-one-out average.

## Panels not reproduced

| Panel | Why |
|---|---|
| 2A cartoon, 2B table | Drawn outside MATLAB. The Figure 2 script reproduces the data traces in 2A and all five properties behind 2B. |
| 3I | Hand-drawn circuit diagrams. |
| 3L, 3O | Mixed-model variance analysis. This is reproduced, with bootstrap intervals, by [tutorial 5](../tutorials/README.md). |
| 4J | Schematic. |
| 6A | Confocal image. The source images are not in the NWB release. |

Figures 4 and 6 stop with an error that says why if you ask for one of these
panels. Figure 3 treats any letter it does not implement as invalid.

## Known differences from the published figures

None of these changes a result. They are listed so that small differences do
not surprise you.

- **Motion-corrected spike times.** Some original scripts removed a few
  spikes during large movements. The motion traces that step needs are not in
  the NWB files, so the scripts use the stored spike times. The step was
  judged unimportant for these figures. In tutorial 2, where it was measured,
  it changes only the after-depolarization of about 3% of cells.
- **Two definitions of the quiet state.** Most panels use the stored
  `quiet_mask`, a two-Gaussian fit to whisking energy. Figure 5B and the state
  bar in Figure 6H keep the original scripts' simpler median rule. The two
  agree on 84% and 88% of frames. Both are saved in `results`.
- **UMAP is random.** UMAP and k-means are not seeded. Figure 2's embedding
  may come out rotated or mirrored; the three clusters come from separate
  hierarchical clustering and do not depend on it. In Figure 6, each run
  relabels the two UMAP groups to match the NPY subtypes, so the confusion
  matrix in panel G cannot swap its columns. The published matrix reproduces
  exactly.
- **Figure 3B and 3D** colour each pair by amplitude and show significant
  pairs as a binary call. The published panels shade by p-value, which is not
  stored.
- **Figure 3 distance bins.** Panel F uses ten 34-µm bins and panel H twelve
  30-µm bins. This is what the analysis code does. The caption says ten 30-µm
  bins.
- **Figure 4A** shows a 10 s window. The published panel is a crop of it.

## Files

| File | Role |
|---|---|
| `reproduce_figureN_from_nwb.m` | The script you run for figure N. Analysis and plotting, one local function per panel. |
| `loadFigureNSessionFromNWB.m` | Reads the NWB file and returns ordinary MATLAB arrays. Figure 5's loader reuses Figure 3's for connectivity files, and Figure 6's reuses Figure 5's. |
| `utils/colorcet.m` | Perceptually uniform colour maps, by Peter Kovesi. |
| `utils/raster_plot.m` | Spike raster drawing. |
| `utils/scatter_kde.m` | Scatter plot coloured by point density, by Nils Haentjens. |

Each script's header lists the original analysis script and line numbers
behind every panel, and every choice kept on purpose.

## How the loaders read the NWB files

The loaders use MATLAB's built-in HDF5 functions (`h5read`, `h5info`), not
MatNWB. They return arrays in the orientation the original analysis used,
for example presynaptic cell × postsynaptic cell × time for pairwise
averages.

Each NWB file carries a table, `/analysis/matlab_variable_lookup`, that maps
every original MATLAB variable name to its place in the file. Comments in the
code use those original names (`traces_all_n`, `pMat`, `spk_t`, and so on).
The table tells you where each one lives.

Two things to know if you write your own reader:

- **Voltage is chunked by time block.** The large voltage arrays are stored in
  blocks of 10,000 frames covering all cells. Reading one cell's whole trace
  still touches every block, so it costs about as much as reading all cells.
  The scripts stream across time instead.
- **Use the stored timestamps.** Recording chunks do not always join on a
  regular frame grid. In the behavior files, for example, they are joined with
  a 20 µs step rather than the real pause between them. Converting seconds to
  frame numbers with `round(t/dt)` can then drift by about a frame per chunk.
  Look times up in the timestamp vector instead.
