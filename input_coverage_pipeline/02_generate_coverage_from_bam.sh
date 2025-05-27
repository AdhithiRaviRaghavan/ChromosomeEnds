#!/bin/bash
#SBATCH --job-name=generate_cov
#SBATCH --output=generate_cov.out
#SBATCH --error=generate_cov.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1

module load samtools/intel/1.20
module load bedtools/intel/2.29.2

mkdir -p bedgraphs

for BAM in *.sorted.bam; do
  PREFIX=${BAM%%.sorted.bam}
  samtools view -b -F 4 "$BAM" | bedtools genomecov -ibam stdin -bg > "bedgraphs/${PREFIX}_coverage.bedgraph"
done
