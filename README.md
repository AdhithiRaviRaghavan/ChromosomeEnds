# ChromosomeEnds

Selected analysis code for:

**Raghavan AR, *et al.* Distinct chromatin regulators downmodulate meiotic axis protein deposition and DNA break induction at chromosome ends.**
*eLife* Reviewed Preprint (2026). [doi:10.7554/eLife.110189.1](https://doi.org/10.7554/eLife.110189.1) · Earlier version: [bioRxiv](https://doi.org/10.1101/2025.02.27.640173)

> **Note:** This repository contains selected analysis code and is a work in progress. The complete code will be made public after final publication.

## Overview

R and SLURM code that takes multiple sequencing datasets from raw alignments to publication figures, studying how chromatin regulators shape meiotic axis proteins and DNA break formation at chromosome ends in budding yeast.

- **HPC pipeline**: alignments to BAM and bedGraph coverage, quantified across repetitive subtelomeric elements and genome-wide bins
- **Normalisation**: genome-wide mean, spike-in normalisation between mutants, and relative coverage for copy-number comparison across strains
- **Statistics**: bootstrapping with 95% confidence intervals, effect sizes and t-tests
- **Profiling**: telomere-distance profiles, metaplots and genome annotation integration

## Repository structure

| Folder | Contents |
|---|---|
| [`input_coverage_pipeline/`](input_coverage_pipeline/) | SLURM coverage and normalisation pipeline |
| [`R/`](R/) | Reusable, documented R functions |
| [`Figures/`](Figures/) | Scripts reproducing Figures 1–3 |

Each folder has its own README with details.

## Usage

Requires R with Bioconductor (`GenomicRanges`, `rtracklayer`, `EnrichedHeatmap`), [`hwglabr2`](https://github.com/hochwagenlab/hwglabr2), `samtools` and `bedtools`.

1. Generate coverage on a cluster: `sbatch input_coverage_pipeline/normalize_relative_coverage_array.sh`
2. Place `.bdg` and `.bed` files in a local `data/` folder.
3. Run a figure script, for example `source("Figures/Figure1.R")`.

## Citation

Please cite the eLife reviewed preprint above if you use or adapt this code.
