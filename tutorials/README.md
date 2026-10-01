# Method tutorials

Five tutorials that teach the analysis methods behind the study, one step at
a time. Each one runs on the public NWB files and checks its own answer
against the numbers stored in them, so you can see that the method you are
learning is the one that produced the published data.

The tutorials teach methods. They do not rebuild manuscript figures; for
that, see [`../figure_reproduction`](../figure_reproduction/README.md).

| | Tutorial | What you learn | NWB file |
|---|---|---|---|
| 1 | [From raw movie to spike times](#tutorial-1--from-raw-movie-to-spike-times) | motion correction, demixing overlapping cells, spike detection with measured error rates | `M_YQ0201_29_behavior_rawexcerpt.nwb` |
| 2 | [Intrinsic properties](#tutorial-2--intrinsic-properties-measured-with-light) | five per-cell properties measured with light instead of an electrode | either `sparseOpto` file |
| 3 | [Synaptic connectivity](#tutorial-3--synaptic-connectivity-and-ipsp-fitting) | detecting inhibitory connections despite a much larger light artifact, and fitting IPSPs | `M_YQ0201_27_sparsePulseHad.nwb` |
| 4 | [Gap junctions](#tutorial-4--spikelets-and-gap-junctions) | detecting electrical coupling from spikelets | `M_YQ0201_27_sparsePulseHad.nwb` |
| 5 | [Pre vs post variance](#tutorial-5--is-the-ipsp-set-by-the-presynaptic-or-the-postsynaptic-cell) | a mixed model that asks whether the sending or the receiving cell sets IPSP size and duration | both `sparsePulseHad` files |

The order is the order of the analysis, but each tutorial stands alone.
Tutorial 5 is the quickest to try: it reads only a small table and runs in
minutes. Stage 1c is the quickest part of tutorial 1.

## Setting up

1. Add the repository to the MATLAB path:

   ```matlab
   addpath(genpath('path/to/this/repository'))
   ```

   Each tutorial script also adds its own folder and `helpers/` when it runs.

2. Download the NWB file the tutorial needs. The file names are those used
   during analysis; the names on DANDI may differ.

3. Tell the tutorial where the file is. How depends on the tutorial:

   | Tutorial | How to give it the file |
   |---|---|
   | 1 | Functions. Pass the path as the first argument. |
   | 2 | Scripts. Set a variable `nwbFile` in the workspace, then run the script. |
   | 3, 4, 5 | Scripts that begin with `clear`. Open the script and edit the `cfg.nwbFile` (or `cfg.nwbFiles`) line near the top. |

The default paths are placeholders (`path/to/...`), so a tutorial stops with
a "file not found" error until you set them.

Every tutorial reads the NWB file with MATLAB's built-in HDF5 functions. You
do not need MatNWB.

## Requirements

MATLAB R2019b or later. The code was tested on R2019b. Toolboxes by
tutorial:

| Toolbox | 1a | 1b | 1c | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|---|
| Statistics and Machine Learning | yes | yes | yes | yes | yes | yes | yes |
| Signal Processing | yes | yes | yes | yes | | | |
| Image Processing | yes | yes | yes | | optional | | |
| Computer Vision | yes | | | | | | |
| Curve Fitting | | | | 2d | yes | | |
| Optimization | | | | | 3d | | |
| Parallel Computing | | | | | optional (3d) | | |

In tutorial 2, Statistics is used only to print comparisons with the
deposit. In tutorial 3, Image Processing is used only if you set
`useMat2gray = true` in stage 3c.

**Memory.** Tutorials 3 and 4 load whole matrices of 320 cells × 3.2 million
frames into memory, about 4 GB each in single precision. Tutorial 3 holds
about 8 GB (voltage and stimulus). Tutorial 4 holds about 10 GB (voltage, a
spike-subtracted copy and trigger masks). Use a machine with at least 16 GB
of RAM; 32 GB is safer.

---

## Tutorial 1 — from raw movie to spike times

**Folder:** `tutorial1_movie_to_spikes/`
**Data:** `M_YQ0201_29_behavior_rawexcerpt.nwb`, about 360 GB. It holds the
raw camera movie for the first four recording chunks (1,152,000 frames, about
24 minutes) of one spontaneous-activity session, plus the processed products
each stage checks itself against. Sixteen cells are marked as tutorial cells
(`/units/tutorial_cell`). The motion-corrected pixels around each of them are
stored, so the demixing can be run on any of the sixteen exactly as it was
run in the original analysis.

Three stages turn camera frames into spike times.

| Stage | Function | What it does | Checked against |
|---|---|---|---|
| 1a | `tutorial1a_motionCorrection` | dense optical-flow registration, and each cell's motion trace | the stored motion-corrected pixels and motion traces |
| 1b | `tutorial1b_demixing` | separates each cell from its overlapping neighbours by non-negative matrix factorisation | the stored footprints and traces |
| 1c | `tutorial1c_spikeFindingAndErrorRates` | matched-filter spike detection, with false-positive and false-negative rates measured from the data | the stored spike times |

```matlab
f = 'path/to/M_YQ0201_29_behavior_rawexcerpt.nwb';
r1a = tutorial1a_motionCorrection(f);               % 1500 frames, about 100 s
r1b = tutorial1b_demixing(f, [6 235]);              % about 60 s per cell
r1c = tutorial1c_spikeFindingAndErrorRates(f, [6 235]);   % about 30 s per cell
```

Each function draws its diagnostic figures and returns a `results`
structure. Pass `false` as the doPlot argument to skip the figures. With no
cell list, 1b uses the first two tutorial cells and 1c the first five.

**The ideas.**

- **1a.** When the mouse runs, the brain moves under the window. A cell's
  fluorescence then changes for two reasons you cannot tell apart later: it
  fired, or it moved. The stage registers every frame to a still reference by
  dense optical flow, since the tissue bends rather than shifting as one
  piece. The flow field, averaged over each cell, gives that cell's motion
  trace. Later stages use it to pick trustworthy frames. The motion is small:
  under one pixel at peak.
- **1b.** Many cells sit within a blur radius of each other, so a simple pixel
  average picks up the neighbours' spikes. In a study about coupling, leaked
  spikes look exactly like coupling. The movie is sign-flipped, so the
  indicator's downward signal becomes positive and a non-negative
  factorisation fits. It is then split into a fast band and a slow band. Each
  cell's footprint is learned on the fast band, where its own spikes make it
  separable, then frozen and used for both bands.
- **1c.** After the sign flip, real spikes are always positive. So the
  *negative* peaks of the filtered trace are a clean sample of that cell's
  noise. Setting the threshold at a chosen rank among the negative peaks fixes
  the false-positive rate in advance. A two-part model of the detected spike
  heights then gives each spike's chance of being noise, and the fraction of
  real spikes missed below threshold.

**What you should see.** Recomputed motion traces match the stored ones at
|r| ≈ 0.996. The corrected pixels match at r ≈ 0.996. For the nearest
neighbour of cell 6, demixing removes about 90% of the leaked spike signal.
Detected spike times match the published ones for 96–98% of spikes. The
mismatches are events within about 2% of the threshold.

**Things to know.**

- **One deliberate change from the published pipeline.** The published rule
  for choosing still frames for the registration template picks every frame
  of the first chunk, so the template came from frames 1–1000. By default,
  stage 1a instead uses the frames where the treadmill registered no movement
  at all. This moves the template to frames 9678–10677 and doubles the residual
  reduction. To reproduce the published template exactly, call
  `buildRegistrationTemplate(f, struct('stillRule', 'smoothedThreshold'))`.
- **Numbers differ slightly from the full session.** Several quantities are
  estimated over the whole 10-chunk recording in the published pipeline:
  motion-free frames, noise level, threshold. On a 4-chunk excerpt they come
  out slightly different. The false-positive budget is given as a rate, so it
  scales with the length of the recording. On the full session it gives
  exactly the published 500 events per cell.
- MATLAB's `opticalFlowFarneback` silently returns zero motion for images
  outside the range 0 to 1. `opticalFlowMotionCorrection` scales its inputs
  and stops with an error if the flow comes back zero.
- Stage 1a is limited by memory rather than time. Its flow field alone is
  1.8 GB for 1500 frames. Lower `nFrames` first if you run short.
- Residual correlation between neighbours after demixing is not a failure.
  These cells are genuinely coupled, so demixing should not drive their
  correlation to zero.

---

## Tutorial 2 — intrinsic properties measured with light

**Folder:** `tutorial2_intrinsic_properties/`
**Data:** either sparseOpto file, `M_YQ0201_27_sparseOpto.nwb` or
`M_YQ0201_29_sparseOpto.nwb`.

Five numbers per cell, all measured optically:

| Property | Meaning |
|---|---|
| maximum firing rate (Hz) | rate at the top of the 5 s ramp |
| pulse adaptation | early vs late firing rate within one 500 ms step |
| optical rheobase | the voltage the cell had reached when it started firing on the ramp |
| after-depolarization | from the average of the cell's own light-evoked spikes |
| membrane time constant (ms) | exponential fit to the voltage after the light turns off |

These are the five properties behind Figure 2 and Figures 5 and 6.

| Stage | Script | Needs | Time |
|---|---|---|---|
| 2a | `tutorial2a_stimulusProtocol` | stimulus only | about 60 s |
| 2b | `tutorial2b_firingRateProperties` | spikes only | about 30 s |
| 2c | `tutorial2c_rheobase` | voltage | about 90 s |
| 2d | `tutorial2d_adpAndMembraneConstant` | voltage | about 60 s |

```matlab
nwbFile = 'path/to/M_YQ0201_27_sparseOpto.nwb';
tutorial2a_stimulusProtocol
tutorial2b_firingRateProperties
tutorial2c_rheobase
tutorial2d_adpAndMembraneConstant
```

**The idea.** An optical recording has no calibrated input and no calibrated
output. Neither is needed. Every trace is divided by that cell's own spike
height, so voltage is measured in "action potentials". Rheobase is expressed
as a voltage the cell reached, not as a light intensity. So every property is
either a rate or a voltage in spike heights.

**What you should see.** Four of the five properties match the values stored
in the NWB file exactly, in both sessions. The after-depolarization matches
exactly for 97% of cells. The rest differ slightly because the published
pipeline removed a few spikes during large movements, a step the public files
cannot support.

**Things to know.**

- **Read time blocks, not cells.** The voltage is stored in blocks of 10,000
  frames covering all 320 cells. Reading one cell's full trace touches every
  block, about 3 GB for a 10 MB trace. The stages therefore stream through
  time and handle all cells at once. Do not rewrite them as a loop over cells.
- **Rheobase uses plain interpolation, not the stored subthreshold trace.** The
  stored `subthreshold_voltage` is smoothed with a Savitzky–Golay filter,
  which bends the slow ramp this measurement reads. Stage 2c removes spikes by
  linear interpolation instead, as the property was defined. It prints both
  versions so you can see the difference.
- In `M_YQ0201_27_sparseOpto.nwb`, seven cells have a stored membrane time
  constant of exactly 0, and two of them pass QC. These are placeholders, not
  measurements. Stage 2d lists them, leaves them out of its comparison, and
  computes real values for them.

---

## Tutorial 3 — synaptic connectivity and IPSP fitting

**Folder:** `tutorial3_ipsp_connectivity/`
**Data:** `M_YQ0201_27_sparsePulseHad.nwb` (stages 3a and 3c also run on
mouse 29).

Stimulate one cell with brief blue pulses and average another cell's voltage
around those pulses. An inhibitory synapse shows up as a small dip. The
problem is that scattered blue light produces a positive artifact in every
cell. It is much larger than the synaptic dip and strongest at short
distances, where synapses are most likely. So the method cannot just look
for a dip.

| Stage | Script | What it teaches | Time |
|---|---|---|---|
| 3a | `tutorial3a_bluePulseSpikeCount` | how many spikes a pulse evokes when individual spikes cannot be seen | 20 s after loading |
| 3b | `tutorial3b_bluePulseSTA` | the pair waveform and its two baseline corrections | 20 s |
| 3c | `tutorial3c_crosstalkTemplate` | building each cell's light-artifact template | 10 s |
| 3d | `tutorial3d_detectAndFitIpsp` | constrained detection against an empirical null, and IPSP kinetics | about 10 s per pair |

Each script first loads the whole voltage matrix, which takes 25 to 100 s.

**How the detection works.** Each pair's waveform is fitted twice. In one fit
it is a baseline plus the light artifact, which may only be positive. In the
other, a synaptic term, which may only be negative, is added. The sign limits
stop the two from absorbing each other. The test statistic is how much the
synaptic term improves the fit. It is judged against the same postsynaptic
cell averaged around pulses delivered more than 800 µm away, which carry the
same artifact but cannot carry a synapse from this cell. Benjamini–Hochberg
correction is then applied across all pairs within 400 µm.

**What you should see.** The spike-per-pulse estimate, pair waveforms and
artifact templates all match the stored values to within rounding. On the
ten example pairs, stage 3d reproduces the published p-values exactly and the
amplitudes to 5 × 10⁻⁵.

**Things to know.**

- **Stage 3a's "spike count" is not a count.** During the pulse the cell sits
  on a large depolarisation and spikes cannot be detected. The number is
  estimated from the power spectrum of the residual instead. For the same
  reason, the artifact template's rising phase is built from the fitted
  decay constant, not a fitted rise.
- **Normalisation.** The stored pair waveforms are divided by the
  postsynaptic cell's spike height and by the presynaptic cell's spikes per
  pulse. The null set built from `normalized_voltage` carries only the first
  of these, so stage 3d divides it by spikes per pulse too. Leave this out
  and every p-value is wrong.
- **A known imperfection, kept on purpose.** In the published pipeline, the
  pulse times for the observed waveform and for its null set are defined one
  frame apart. Stage 3d keeps this so its results match the published ones.
  To correct it, pass `data.t_blue` instead of `tBlueNull` to
  `generatePairNullSet` in stage 3d. The comment at that line explains. Over
  the full family of 9,124 pairs, correcting it changes 138 of about 950
  connection calls.
- **Only the independent kinetic fits are used.** The fitting code also
  computes a pooled (empirical-Bayes) fit, but no published number comes from
  it.
- **Subsets.** The stages run on a few cells or pairs so that they finish in
  minutes. Benjamini–Hochberg depends on the whole family of pairs, so
  detection calls from a subset are compared with the stored calls, not
  recomputed. A full run of all 9,124 pairs takes about 10 hours on 8 workers.
- `cfg.referenceMat` and `cfg.resultsMat` are optional comparison files from
  the original analysis. They are not in the public release. Leave them as
  `''`.

---

## Tutorial 4 — spikelets and gap junctions

**Folder:** `tutorial4_gap_junction_spikelets/`
**Data:** `M_YQ0201_27_sparsePulseHad.nwb`.

Two cells joined by a gap junction pass each other a **spikelet**: a small,
smoothed copy of one cell's spike that appears in the other within a frame or
two. To make a cell spike you must light it, and that light also reaches its
neighbours. So a deflection in cell *j* after cell *i* spikes could be
coupling, or cell *j* responding to the same light.

The fix is the Hadamard block. In every fourth recording chunk, each cell is
lit on and off in its own pattern, so at any moment about half the field is
lit. For any pair there are many frames where the first cell is lit and the
second is dark. Averaging only over those removes shared light as an
explanation.

| Stage | Script | What it teaches | Time |
|---|---|---|---|
| 4a | `tutorial4a_orthogonalTriggers` | choosing spikes when only the presynaptic cell is lit | about 2 min |
| 4b | `tutorial4b_spikeletWaveform` | the amplitude measure, and the two artifacts that survive | about 2.5 min |
| 4c | `tutorial4c_commonModeCorrection` | removing a distance trend and a field-wide common signal | about 2 min |
| 4d | `tutorial4d_spikeletDetection` | testing against spikes from far-away cells, then Benjamini–Hochberg | about 3 min |

**The key idea (4c).** Two artifacts survive the Hadamard trick. One is a
smooth trend with distance, from the presynaptic cell's own light spreading.
The other is a common signal shared by all cells. They are removed in turn,
and **both fits leave out the peak**. Otherwise a real spikelet in the
nearest distance bin would be fitted as part of the trend and subtracted
away.

**What you should see.** The spike sets match the stored ones with no spikes
missing (0.01% extra, from the movement step the public files cannot
support). The waveforms match to 10⁻⁷, the corrected amplitudes to
r > 0.9999, and the detection decisions 12 of 12.

**Things to know.**

- **Memory, not time, is the limit.** Keep traces in single precision. Do not
  run the pair loop in `parfor`: it copies the whole working set to every
  worker and runs out of memory.
- **The stored `hadamard_cross_spike_sta` stops at 400 µm.** Pairs further
  apart were never computed and are stored as zeros. The common-mode
  correction needs cells beyond 800 µm, so this tutorial recomputes the
  averages it needs instead of reading them.
- `runFullSpikeletDetection(f, outDir)` runs the whole family of 9,124 pairs,
  which takes about 43 hours. `runFullSpikeletDetection(f, outDir, true)` is a
  3-minute test run.

---

## Tutorial 5 — is the IPSP set by the presynaptic or the postsynaptic cell?

**Folder:** `tutorial5_pre_post_variance/`
**Data:** both sparsePulseHad files. Only the small pairwise table is read, so
loading takes about a second.

Tutorial 3 gives each connected pair an IPSP amplitude and a decay time. This
tutorial asks where their variation comes from: the cell sending the IPSP,
the cell receiving it, or the particular pair. This is the analysis behind
Figure 3L and 3O.

Each pair's log amplitude (or log decay) is modelled with fixed effects for
mouse and distance, plus a random offset for the presynaptic cell and one for
the postsynaptic cell. This is a linear mixed model with crossed random
intercepts. It estimates how much the offsets vary across cells: σ²_pre,
σ²_post and a pair-specific remainder. A parametric bootstrap gives 95%
intervals.

| Stage | Script | Response | Result |
|---|---|---|---|
| 5a | `tutorial5a_amplitudeVariance` | IPSP amplitude | σ²_pre 0.144 [0.116, 0.175] > σ²_post 0.058 [0.042, 0.075] |
| 5b | `tutorial5b_decayVariance` | IPSP decay, weighted by fit precision | σ²_pre 0.021 [0, 0.050] < σ²_post 0.111 [0.070, 0.156] |

Amplitude follows mainly the presynaptic cell; decay follows mainly the
postsynaptic cell. In both, none of the 2,000 bootstrap replicates crossed
zero, so p ≤ 0.001.

**What you should see.** These numbers exactly. The tutorial is
bit-identical to the published analysis at every step, down to each bootstrap
replicate. Each 2,000-replicate bootstrap takes about 7 minutes. Set
`cfg.nBootstrap = 200` for a quick look.

**Things to know.**

- The loader keeps cells that passed quality control and drops connections
  shorter than 60 µm, as the published analysis did.
- Some decay fits stopped at a bound (6 ms or 400 ms). They are not dropped;
  their large standard errors give them less weight. Stage 5b shows that the
  conclusion does not depend on the weighting.
- The source code gave its p-values as a typed-in constant. The tutorial
  computes a two-sided bootstrap p from the replicates instead.

---

## Helper functions

Tutorials 1, 3 and 4 call functions from the analysis pipeline. Copies are in
`helpers/`. They are unchanged apart from a header saying where each came
from.

| Folder | Functions | Origin | Used by |
|---|---|---|---|
| `helpers/voltage_imaging/` | `get_sta_mat_raw`, `get_sta_mat_2`, `get_sta_mat_single`, `fitSelfBlueSta`, `whitened_matched_filter`, `adaptive_thresh`, `spike_trace`, `cellSpkFtprntInit2` | the author's voltage-imaging library. `whitened_matched_filter` and `adaptive_thresh` implement the method of VolPy (Cai et al., 2021, *PLoS Comput Biol* 17: e1008806). | 1, 3, 4 |
| `helpers/psp_fitting/` | 16 functions for IPSP detection and fitting, including `benjaminiHochberg` | written for this study | 3, 4 |
| `helpers/caiman_gpl/` | `HALS_temporal`, `HALS_spatial`; `ALS_nonneg`, `ALS_nonneg_BackgroundUpdateOnly` | [CaImAn-MATLAB](https://github.com/flatironinstitute/CaImAn-MATLAB); the ALS files are the author's modifications of it. **GNU GPL v2**, see `license.txt` there. | 1b |
| `helpers/cohen_lab/` | `toimg`, `tovec`; `SeeResiduals_vec` | Cohen Lab, Harvard University; `SeeResiduals_vec` is the author's modification of the lab's `SeeResiduals` | 1b, 3, 4 |
| `helpers/ironclust/` | `fft_clean` | [IronClust](https://github.com/flatironinstitute/ironclust) | 1b, 1c |

## Files in each tutorial folder

- `tutorialNx_<name>.m`: the scripts or functions you run, one per stage.
  Each is written as a narrative, with comments that explain the method.
- camelCase `.m` files: reusable analysis functions. They draw nothing.
- `plot...Diagnostics.m`: the figures. They compute nothing that changes a
  result.
- `load...ForX.m`: reads the NWB file and returns arrays named as in the
  original analysis.

Comments in the code sometimes refer to the original analysis scripts or to
internal notes by file name. Those files are not part of this repository;
the comments are kept as a record of where each step came from.
