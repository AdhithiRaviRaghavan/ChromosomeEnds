#!/bin/bash
#SBATCH --job-name=relcov_array                        # Job name
#SBATCH --output=relcov_array_%A_%a.out                # Standard output
#SBATCH --error=relcov_array_%A_%a.err                 # Standard error
#SBATCH --array=0-7                                    # Array job index (update based on # of input files - 1)
#SBATCH --time=02:00:00                                # Time limit
#SBATCH --cpus-per-task=2                              # CPUs per task

# Load required modules
module load samtools/intel/1.20
module load bedtools/intel/2.29.2

# Define directories
SAM_DIR="/scratch/arr609/Manuscript_Sir3/MultipleMapping/SingleEndReads/AllSams"             # Input SAM files
BEDGRAPH_DIR="${SAM_DIR}/AllBedGraphs"                                                       # Output directory for bedGraphs
BED_DIR="/scratch/arr609/Manuscript_Sir3/MultipleMapping/ReadsCounts/region_beds"            # BED files for Y', X, and subtelomeric regions
OUT_DIR="${SAM_DIR}/AllTsv"                                                                  # Output directory for TSV results

# Create output directories if they don't exist
mkdir -p "$BEDGRAPH_DIR" "$OUT_DIR"

# Get the list of SAM files and select one based on array index
SAMS=($(ls "$SAM_DIR"/*.sam))
SAM=${SAMS[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$SAM" .sam)
echo "🔄 Processing task $SLURM_ARRAY_TASK_ID: $BASENAME"

# Step 1: Convert SAM to sorted BAM and index
BAM_SORTED="${BEDGRAPH_DIR}/${BASENAME}.sorted.bam"
samtools view -bS "$SAM" | samtools sort -o "$BAM_SORTED"
samtools index "$BAM_SORTED"

# Step 2: Generate a genome-wide coverage bedGraph from mapped reads only (-F 4 filters out unmapped reads)
BEDGRAPH_FILE="${BEDGRAPH_DIR}/${BASENAME}_MACS_input_coverage.bedgraph"
samtools view -b -F 4 "$BAM_SORTED" | \
  bedtools genomecov -ibam - -bg > "$BEDGRAPH_FILE"
echo "✅ BedGraph generated: $BEDGRAPH_FILE"

# Step 3: Compute genome-wide mean coverage (used for normalization)
# Weighted average: sum(coverage * region length) / total length
GENOME_MEAN=$(awk '{sum+=($3-$2)*$4; len+=($3-$2)} END {print sum/len}' "$BEDGRAPH_FILE")
echo "📈 Genome-wide mean coverage: $GENOME_MEAN"

# Step 4: For each genomic region (Y′ elements, X elements, subtelomeric 20 kb), compute:
# - Mean input coverage from bedGraph
# - Relative coverage = region_mean / genome_mean
for REGION in y_elements x_elements subtel_20kb; do
  REGION_BED="${BED_DIR}/${REGION}.bed"
  OUT_FILE="${OUT_DIR}/${BASENAME}_${REGION}_mean_coverage.tsv"

  # bedtools map: extracts mean coverage over specified regions
  # awk: compute relative coverage, set to 0 if no coverage
  bedtools map -a "$REGION_BED" -b "$BEDGRAPH_FILE" -c 4 -o mean | \
    awk -v GMEAN="$GENOME_MEAN" 'BEGIN{OFS="\t"} 
      {rel_cov = ($4=="." || $4=="") ? 0 : $4/GMEAN; 
       print $1, $2, $3, $4, rel_cov}' > "$OUT_FILE"

  echo "✅ Saved: $OUT_FILE"
done

echo "🎉 Task $SLURM_ARRAY_TASK_ID complete!"
