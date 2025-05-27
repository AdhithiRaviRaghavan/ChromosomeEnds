#!/bin/bash
#SBATCH --job-name=per_region_cov
#SBATCH --output=per_region_cov.out
#SBATCH --error=per_region_cov.err
#SBATCH --time=01:00:00

module load bedtools/intel/2.29.2

BEDGRAPH="combined_input.bedgraph"
BED_DIR="region_beds"
OUT_DIR="per_element_coverage"
mkdir -p "$OUT_DIR"
REGIONS=("x_elements.bed" "y_elements.bed" "subtel_20kb.bed")

for BEDFILE in "${REGIONS[@]}"; do
  BASENAME=$(basename "$BEDFILE" .bed)
  BED="$BED_DIR/$BEDFILE"
  OUT="$OUT_DIR/${BASENAME}_coverage.tsv"
  echo -e "Region_ID\tChrom\tStart\tEnd\tLength\tMapped_Bases\tCoverage_X" > "$OUT"

  bedtools intersect -a "$BEDGRAPH" -b "$BED" -wa -wb | \
  awk '{
    len = $3 - $2;
    region = $5 ":" $6 "-" $7;
    cov = len * $4;
    region_len[region] = $7 - $6;
    chrom[region] = $5;
    start[region] = $6;
    end[region] = $7;
    mapped[region] += cov;
  }
  END {
    for (r in region_len) {
      cov_x = mapped[r] / region_len[r];
      print r, chrom[r], start[r], end[r], region_len[r], mapped[r], cov_x;
    }
  }' OFS='\t' | sort -k2,2 -k3,3n >> "$OUT"
done
