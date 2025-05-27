#!/bin/bash
#SBATCH --job-name=internal_bins
#SBATCH --output=internal_bins.out
#SBATCH --error=internal_bins.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

module load samtools/intel/1.20
module load bedtools/intel/2.29.2

FASTA="Yue.SK1.genome.fa"
BEDGRAPH="combined_input.bedgraph"
BED_DIR="region_beds"
OUTDIR="per_element_coverage"

samtools faidx "$FASTA"
cut -f1,2 "${FASTA}.fai" > SK1Yue.genome

cat "$BED_DIR"/*.bed | sort -k1,1 -k2,2n | bedtools merge -i - > merged_special.bed

bedtools makewindows -g SK1Yue.genome -w 20000 > genome_bins_20kb.bed
bedtools subtract -a genome_bins_20kb.bed -b merged_special.bed > internal_20kb_bins.bed
bedtools map -a internal_20kb_bins.bed -b "$BEDGRAPH" -c 4 -o mean > "$OUTDIR/internal_20kb_cov.bed"

awk 'BEGIN{OFS="\t"} {print $1 ":" $2 "-" $3, $1, $2, $3, $4}' "$OUTDIR/internal_20kb_cov.bed" > "$OUTDIR/internal_20kb_coverage.tsv"
