#!/bin/bash
#SBATCH --job-name=sam_to_bam
#SBATCH --output=sam_to_bam.out
#SBATCH --error=sam_to_bam.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1

module load samtools/intel/1.20

for INPUT_SAM in *MACS_input.sam; do
  PREFIX=$(basename "$INPUT_SAM" .sam)
  samtools view -bS "$INPUT_SAM" > "${PREFIX}.bam"
  samtools sort "${PREFIX}.bam" -o "${PREFIX}.sorted.bam"
  samtools index "${PREFIX}.sorted.bam"
done
