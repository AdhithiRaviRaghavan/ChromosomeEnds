# ChromosomeEnds

Code used for the analysis in:

**Raghavan et al., _"Multiple mechanisms suppress the initiation of meiotic recombination near chromosome ends"_**  
bioRxiv 2025.02.27.640173 — [https://doi.org/10.1101/2025.02.27.640173](https://doi.org/10.1101/2025.02.27.640173)

---

## Overview

This repository contains all scripts and reusable functions used to generate the figures and results in the above manuscript.  
The work focuses on understanding how axis protein distribution and chromatin structure suppress recombination at chromosome ends during meiosis.

---

## Repository Structure

```
ChromosomeEnds/
├── Figures/               # R scripts to generate each figure (e.g., Figure1.R)
├── R/                     # Modular utility functions (e.g., signal normalization, metaplots)
├── input_coverage_pipeline/ # SLURM-ready scripts to generate input bedGraph & region coverage
├── README.md              # This file

```

---

## Setup & Usage

To reproduce figure panels or run individual analyses:

1. Place your file -> `.bdg` or `.bed` files in a local `data/` directory (not uploaded).
2. Run figure scripts from the `Figures/` folder.
3. Utility functions are automatically loaded via `source("R/filename.R")`.

Most analyses depend on:
- [`hwglabr2`](https://github.com/hochwagenlab/hwglabr2)
- `GenomicRanges`, `EnrichedHeatmap`, `ggplot2`, `rtracklayer`, and other Bioconductor tools

---

## Citation

If you use or adapt this code, please cite:

> Raghavan et al., 2025. _Multiple mechanisms suppress the initiation of meiotic recombination near chromosome ends._ bioRxiv. https://doi.org/10.1101/2025.02.27.640173

## Note

This repository is a working version of the code used in the manuscript.  
Figure scripts and associated functions may change as the manuscript goes through peer review.

