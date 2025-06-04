# Input Coverage Pipeline

This folder contains SLURM-ready scripts used to process input ChIP-seq coverage profiles for the chromosome ends study.

## Overview

The goal is to:
- Convert SAM files to sorted, indexed BAMs
- Generate genome coverage (bedGraph format)
- Quantify average signal in genomic features: X elements, Y' elements, subtelomeric 20 kb, and internal regions
- Normalize regional coverage to genome-wide signal for relative enrichment or copy number estimation

## Folder Contents

| Script | Description |
|--------|-------------|
| `01_convert_sam_to_sorted_bam.sh` | Converts `*_MACS_input.sam` files to sorted and indexed `.bam` format |
| `02_generate_coverage_from_bam.sh` | Creates coverage bedGraph from sorted BAM files |
| `03_coverage_per_region.sh` | Calculates coverage for X, Y′, and subtelomeric regions |
| `04_internal_20kb_bins_coverage.sbatch` | Quantifies average coverage in 20 kb internal genome bins |
| `05_expand_bedgraph_per_base.sbatch` | Expands bedGraph to per-base resolution for fine-scale plots |
| `region_beds/` | Folder of BED files specifying genomic regions of interest |
| `normalize_relative_coverage_array.sh` | Full pipeline: starts from SAM files, converts to BAM, generates bedGraph, and calculates region-wise normalized coverage (relative to genome-wide mean). Useful for comparing signal or copy number across strains.

## Dependencies

- `samtools` ≥ 1.20
- `bedtools` ≥ 2.29.2
- SLURM job submission

## Usage

Run each script sequentially on your HPC cluster using:

For example: 

```bash
sbatch 01_convert_sam_to_sorted_bam.sh
