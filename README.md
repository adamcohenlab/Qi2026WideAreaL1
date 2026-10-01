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

*Under construction.* The manuscript is not yet public. Citation details will
be added here when it is.

## Data

https://dandiarchive.org/dandiset/001960

Every script reads the NWB files directly. Nothing else needs to be
downloaded. The files are large (about 20 to 30 GB per session, and about 360 GB
for the raw-movie excerpt used in tutorial 1), so keep them on a fast local
drive.

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

## Requirements

- MATLAB R2019b or later.
- MathWorks toolboxes. Which ones you need depends on the script; each
  README has a table. The full set is: Statistics and Machine Learning,
  Signal Processing, Image Processing, Curve Fitting, Optimization, and
  Computer Vision. Parallel Computing is optional.
- The [UMAP package for MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)
  from File Exchange, for the Figure 2 and Figure 6 embeddings only.

## Third-party code

The repository includes a few functions written by others. Each file
says where it came from at the top.

- `figure_reproduction/utils/`: `colorcet` (Peter Kovesi) and `scatter_kde`
  (Nils Haentjens).
- `tutorials/helpers/`: code from CaImAn-MATLAB, IronClust and the Cohen Lab
  at Harvard. The CaImAn-derived files are under the GNU GPL v2; its licence
  text is in `tutorials/helpers/caiman_gpl/license.txt`. See the tutorials
  README for the full list.

## Licence

*Under construction.*
