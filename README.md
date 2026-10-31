# All-optical electrophysiology of layer 1 interneurons

MATLAB code for a study of layer 1 (L1) interneurons in mouse barrel cortex,
recorded with wide-field, all-optical electrophysiology in awake mice.

In each session about 320 NDNF-positive L1 interneurons were imaged at once
with the voltage indicator Voltron, while a digital micromirror device used the
channelrhodopsin CheRiff to stimulate chosen cells with blue light. The mice
ran freely on a treadmill, and a face camera recorded whisking. From these
recordings the study measures:

- each cell's intrinsic excitability, optically, with no electrode;
- inhibitory synaptic connections and electrical (gap-junction) coupling
  between hundreds of cell pairs at once;
- how intrinsic properties and connectivity shape spontaneous network activity;
- how these differ between the NPY-positive and NPY-negative subtypes, and
  how brain state (quiet versus whisking) changes them.

## What is in this repository

| Folder | What it is for |
|---|---|
| [`figure_reproduction/`](figure_reproduction/README.md) | Rebuild the data panels of manuscript Figures 2 to 6 from the public NWB files. |
| [`tutorials/`](tutorials/README.md) | Five step-by-step tutorials that teach the analysis methods: from raw movie to spike times, intrinsic properties, synaptic connectivity, gap junctions, and pre- versus postsynaptic variance. |

The two parts are independent. Use the figure scripts to check a published
result. Use the tutorials to learn a method, or to apply it to your own data.
Each folder has its own README with instructions.

## Citation

If you use this code or the data, please cite the preprint:

Qi Y, Parot VJ, Wu X, Plutkis SE, Lavis LD, Cohen AE. Large-scale all-optical
electrophysiology links excitability and connectivity to circuit dynamics in
vivo. *bioRxiv* (2026).
[doi:10.64898/2026.10.05.756595](https://doi.org/10.64898/2026.10.05.756595)

```bibtex
@article{Qi2026,
  title     = {Large-scale all-optical electrophysiology links excitability and connectivity to circuit dynamics in vivo},
  author    = {Qi, Yitong and Parot, Vicente J. and Wu, Xiang and Plutkis, Sarah E. and Lavis, Luke D. and Cohen, Adam E.},
  journal   = {bioRxiv},
  year      = {2026},
  doi       = {10.64898/2026.10.05.756595},
  url       = {https://www.biorxiv.org/content/10.64898/2026.10.05.756595v1}
}
```

## Data

https://dandiarchive.org/dandiset/001960

Every script reads the NWB files directly. Nothing else needs to be
downloaded. The files are large, so keep them on a fast local drive.

The READMEs and scripts call the files by the names used during analysis.
On DANDI the same files have the standard DANDI names:

| Name in this repository | Name on DANDI | Size |
|---|---|---|
| `M_YQ0201_27_sparseOpto.nwb` | `sub-M-YQ0201-27/sub-M-YQ0201-27_ses-sparseOpto_behavior+ophys.nwb` | 24.5 GB |
| `M_YQ0201_29_sparseOpto.nwb` | `sub-M-YQ0201-29/sub-M-YQ0201-29_ses-sparseOpto_behavior+ophys.nwb` | 19.9 GB |
| `M_YQ0201_27_sparsePulseHad.nwb` | `sub-M-YQ0201-27/sub-M-YQ0201-27_ses-sparsePulseHad_behavior+ophys.nwb` | 29.7 GB |
| `M_YQ0201_29_sparsePulseHad.nwb` | `sub-M-YQ0201-29/sub-M-YQ0201-29_ses-sparsePulseHad_behavior+ophys.nwb` | 25.5 GB |
| `M_YQ0201_27_behavior.nwb` | `sub-M-YQ0201-27/sub-M-YQ0201-27_ses-behavior_behavior+ophys.nwb` | 26.4 GB |
| `M_YQ0201_29_behavior.nwb` | `sub-M-YQ0201-29/sub-M-YQ0201-29_ses-behavior_behavior+ophys.nwb` | 22.7 GB |
| `M_YQ0201_29_behavior_rawexcerpt.nwb` | `sub-M-YQ0201-29/sub-M-YQ0201-29_ses-behavior-rawexcerpt_behavior+ophys.nwb` | 356.8 GB |

You can keep the DANDI names. Just pass the path to the file you downloaded.

## Getting started

1. Install MATLAB R2019b or later. The code was written and tested on R2019b.
2. Download or clone this repository.
3. In MATLAB, add the whole repository to the path:

   ```matlab
   addpath(genpath('path/to/this/repository'))
   ```

4. Download the NWB files you need. Each README lists which files each
   script reads.
5. Follow the README in `figure_reproduction/` or `tutorials/`.

The NWB files are read with MATLAB's built-in HDF5 functions, so you do not
need MatNWB or any other NWB package.

**Install time.** Installing the code takes a few minutes. There is nothing
to compile: you download it and add it to the MATLAB path. Installing MATLAB
and its toolboxes, if you do not already have them, is separate. Downloading
the data takes longer: one 25 GB session takes about 35 minutes at
100 Mbit/s.

## Requirements

**Software**

- Tested on Windows 11. The code is plain MATLAB with no compiled files, so
  it should also run on macOS and Linux, but this has not been tested.
- MATLAB R2019b or later. All code was tested on R2019b (9.7).
- MathWorks toolboxes. Which ones you need depends on the script; each
  README has a table. Tested with the versions that ship with R2019b:

  | Toolbox | Version tested |
  |---|---|
  | Statistics and Machine Learning | 11.6 |
  | Signal Processing | 8.3 |
  | Image Processing | 11.0 |
  | Curve Fitting | 3.5.10 |
  | Optimization | 8.4 |
  | Computer Vision | 9.1 |
  | Parallel Computing (optional) | 7.1 |

- The [UMAP package for MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)
  from File Exchange, for the Figure 2 and Figure 6 embeddings only. Tested
  with the January 2020 release.

**Hardware**

No special hardware is needed. The code does not use a GPU.

- **Memory:** 16 GB of RAM is enough for most scripts. Tutorials 3 and 4
  hold 8 to 10 GB of data at once, so 32 GB is safer for them.
- **Disk:** 20 to 30 GB for each session file you download, and about 360 GB
  for the raw-movie excerpt used in tutorial 1.
- **Tested on:** a desktop PC with an Intel Core i9-10850K (10 cores),
  128 GB of RAM and Windows 11. The run times in the other READMEs were
  measured on this machine, with the NWB files on a local drive.

## Using the code on your own data

The figure scripts and the tutorial scripts read these NWB files at fixed
internal paths, so as written they run only on this dataset. The method
functions they call take plain MATLAB arrays, so you can use them on other
voltage-imaging data:

| Step | Function | Main input |
|---|---|---|
| Motion correction | `opticalFlowMotionCorrection` (tutorial 1) | a movie (rows × columns × frames) and a reference image |
| Demixing overlapping cells | `demixCellNMF` (tutorial 1) | the motion-corrected pixels around one cell, and starting footprints for it and its neighbours |
| Spike detection | `detectSpikesMatchedFilter` (tutorial 1) | one cell's voltage trace |
| Spike-detection error rates | `fitSpikeHeightMixture` (tutorial 1) | the output of `detectSpikesMatchedFilter` |
| Spike-triggered averages | `get_sta_mat_2`, `get_sta_mat_single` (`tutorials/helpers/voltage_imaging/`) | spike times and voltage traces |
| IPSP detection and fitting | functions in `tutorials/helpers/psp_fitting/` | averaged pair waveforms, stimulus times and cell positions |
| Correcting for many tests | `benjaminiHochberg` (`tutorials/helpers/psp_fitting/`) | a list of p-values |

Type `help` followed by the function name in MATLAB to see its inputs and
options. The easiest way to start is to run the matching tutorial on the
public data first. Then replace the lines where the tutorial loads arrays
from the NWB file with your own arrays.

## Third-party code

The repository includes a few functions written by others. Each file
says where it came from at the top.

- `figure_reproduction/utils/`: `colorcet` (Peter Kovesi) and `scatter_kde`
  (Nils Haentjens).
- `tutorials/helpers/`: code from CaImAn-MATLAB, IronClust and the Cohen Lab
  at Harvard. The CaImAn-derived files are under the GNU GPL v2; its licence
  text is in `tutorials/helpers/caiman_gpl/license.txt`. The IronClust file is
  under the Apache License 2.0; its licence text is in
  `tutorials/helpers/ironclust/LICENSE`. See the tutorials README for the full
  list.

## Licence

The code written for this repository is released under the
[MIT License](LICENSE). You may use, copy, change and share it, including
for commercial use, as long as you keep the copyright notice and licence text.

The third-party files listed above keep their own licences. The
CaImAn-derived files in `tutorials/helpers/caiman_gpl/` are under the GNU GPL
v2, so if you redistribute those files, or code built on them, the GPL terms
apply to that part.
