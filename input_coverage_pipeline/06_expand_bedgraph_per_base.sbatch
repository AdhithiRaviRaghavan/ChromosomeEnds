#!/bin/bash
#SBATCH --job-name=expand_bed
#SBATCH --output=expand_bed.out
#SBATCH --error=expand_bed.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G

module load bedtools/intel/2.29.2

BEDGRAPH="combined_input.bedgraph"
OUTDIR="per_base_coverage"
mkdir -p "$OUTDIR"

awk '{
  for (i = $2; i < $3; i++) {
    print $1, i, i+1, $4
  }
}' "$BEDGRAPH" > "$OUTDIR/expanded_input.bed"
