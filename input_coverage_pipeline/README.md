# Input Coverage Pipeline

This folder contains SLURM-ready scripts used to process input ChIP-seq coverage profiles for the chromosome ends study.

## Overview

The goal is to:
- Convert SAM files to sorted, indexed BAMs
- Generate genome coverage (bedGraph format)
- Quantify average signal in genomic features: X elements, Y' elements, subtelomeric 20 kb, and internal regions

## Folder Contents

| Script | Description |
|--------|-------------|
| `01_convert_sam_to_bam.sbatch` | Converts `*_MACS_input.sam` files to sorted and indexed `.bam` format |
| `02_generate_bedgraph.sh` | Creates coverage bedGraph from sorted BAM files |
| `03_calculate_region_coverage.sbatch` | Calculates coverage for X, Y′, and subtelomeric regions |
| `04_internal_coverage_20kb.sbatch` | Quantifies average coverage in 20 kb internal genome bins |
| `05_per_base_coverage.sbatch` | Expands bedGraph to per-base resolution for fine-scale plots |
| `region_beds/` | Folder of BED files specifying genomic regions of interest |

## Dependencies

- `samtools` ≥ 1.20
- `bedtools` ≥ 2.29.2
- SLURM job submission

## Usage

Run each script sequentially on your HPC cluster using:

```bash
sbatch 01_convert_sam_to_bam.sbatch
