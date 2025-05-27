################################################################################
# Figure2.R — Full Panel Script for Figure 2
# Manuscript: "Multiple mechanisms suppress the initiation of meiotic recombination near chromosome ends"
# Author: Adhithi R. Raghavan
# Created: January 2025
#
# Panels:
# - 2A: Telomeric signal depletion at X-with-Y′ vs. X-only ends
# - 2B: Signal in synthetic fusion chromosomes (ChrIV-I, Cen1D/Cen4D/WT)
# - 2C: Coding density vs. telomere distance
#
# Notes:
# - This script assumes local .bdg and genome annotation files in user-defined paths
# - Utility functions like `gendiv()` should be sourced from R/normalize_signal.R
# - This script is provided as a working reference — subject to updates during manuscript review
################################################################################

#Figure 2:Multiple cis-acting features contribute to axis protein depletion near telomeres

#Load necessary libraries:
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
library(reshape2)
library(tibble)
library(dplyr)
library(magrittr)
library(pbapply)
library(readr)
library(cowplot)

# Define local input file paths
hop1_file <- "data/Hop1_WT.bdg"
red1_file <- "data/Red1_WT.bdg"
rec8_file <- "data/Rec8_WT.bdg"
cen1D_fusion_file <- "data/Cen1D_Fusion.bdg"
red1_fusion_wt_file <- "data/Red1_Fusion_WT.bdg"
cen4D_fusion_file <- "data/Cen4D_Fusion.bdg"

# Import ChIP signal data
Hop1_WT <- import_bedGraph(hop1_file)
Red1_WT <- import_bedGraph(red1_file)
Rec8_WT <- import_bedGraph(rec8_file)

cen1D_Fusion <- import_bedGraph(cen1D_fusion_file)
Red1_Fusion_WT <- import_bedGraph(red1_fusion_wt_file)
cen4D_Fusion <- import_bedGraph(cen4D_fusion_file)


######################################################################################
# Figure 2A: Distance from telomeres seperated by X only and XwithY' ends
######################################################################################
# Function: Normalize signal to genome average
source("R/normalize_signal.R")

# Normalize the genome signal
Hop1_WT <- gendiv(Hop1_WT)
Red1_WT <- gendiv(Red1_WT)
Rec8_WT <- gendiv(Rec8_WT)

# Function: Extract telomeric signal from normalized input
# Used in Figure 2A, assumes input already normalized using gendiv().
teloSeqSignal_simple <- function(Sample, length_to_collect = 120000, Genome = 'SK1Yue') {
  Sample_telo <- hwglabr2::signal_from_telomeres2(Sample, length_to_collect, genome = Genome)
  return(Sample_telo)
}

# Extract telomere signals for each sample
hop1_telo <- teloSeqSignal_simple(Hop1_WT)
red1_telo <- teloSeqSignal_simple(Red1_WT)
rec8_telo <- teloSeqSignal_simple(Rec8_WT)

######################################################################################
# X-with-Y' Ends: Distance from Telomeres
######################################################################################
XY_Ends_DistanceFromTelo <- function(hop1_telo, gff_data = 'SK1Yue', length_to_collect = 120000) {
  # Load GFF data
  gff <- get_gff(gff_data)
  xelement <- gff[grep('X_element', gff$type)]
  yelement <- gff[grep('Y_prime_element', gff$type)]
  
  # Identify X elements with Y' elements
  xelement_df <- data.frame(xelement)
  yelement_df <- data.frame(yelement)
  XwithYrows_all <- merge(
    xelement_df[xelement_df$type == "X_element", ],
    yelement_df[yelement_df$type %in% c("Y_prime_element", "Y_prime_element_partial"), ],
    by = c("seqnames", "strand")
  )
  
  # Adjust strand information
  XwithYrows_all$strand <- as.character(XwithYrows_all$strand)
  XwithYrows_all$strand[is.na(XwithYrows_all$strand)] <- "+"
  XwithYrows_all$strand[XwithYrows_all$strand == "-"] <- "L"
  XwithYrows_all$strand[XwithYrows_all$strand == "+"] <- "R"
  
  # Filter telomere data based on strand and chromosome
  chromosomes_strand <- XwithYrows_all[, c("seqnames", "strand")]
  hop1_telo_df <- data.frame(seqnames = hop1_telo$chr, strand = hop1_telo$arm)
  Sample_telo2 <- merge(hop1_telo, chromosomes_strand, by.x = c("chr", "arm"), by.y = c("seqnames", "strand"))
  
  # Average signal from different chromosomes
  Sample_telo1 <- Sample_telo2[, 4:ncol(Sample_telo2)]
  Sample_teloGA <- colMeans(Sample_telo1, na.rm = TRUE)
  
  # Compress and smooth data
  Sample_teloGADF <- data.frame(position = seq(1, length_to_collect), value = Sample_teloGA)
  Sample_teloCompressed <- hwglabr2::compress_signal_track(Sample_teloGADF, window_size = 200)
  Sample_telo <- ksmooth(x = Sample_teloCompressed$position, y = Sample_teloCompressed$window_mean, bandwidth = 25000)
  
  return(list(Sample_telo, chromosomes_strand))
}

# Process X-with-Y' ends
XY_hop1_telo <- XY_Ends_DistanceFromTelo(hop1_telo)
XY_red1_telo <- XY_Ends_DistanceFromTelo(red1_telo)
XY_rec8_telo <- XY_Ends_DistanceFromTelo(rec8_telo)

# Plot X-with-Y' ends
plot_XY <- function(data1, data2, data3) {
  plot(data1[[1]], col = "forestgreen", xlab = "Distance to chr end", ylab = "Average ChIP signal",
       main = "XY' ends", type = "l", bty = "n", lwd = 5, ylim = c(0.4, 1.6), yaxp = c(0.4, 1.6, 6))
  lines(data2[[1]], type = "l", lwd = 5, col = "red")
  lines(data3[[1]], type = "l", lwd = 5, col = "purple")
  abline(h = 1, col = "gray60", lwd = 3, lty = 3)
  legend(50000, 1.79, c("Red1", "Hop1", "Rec8"), text.col = c("forestgreen", "red", "purple"), bty = "n")
}

p1 <- as.ggplot(~plot_XY(XY_red1_telo, XY_hop1_telo, XY_rec8_telo))
p1
######################################################################################
# X-Only Ends: Distance from Telomeres
######################################################################################
XOnly_Ends_DistanceFromTelo <- function(hop1_telo, gff_data = 'SK1Yue', length_to_collect = 120000) {
  gff <- get_gff(gff_data)
  xelement <- gff[grep('X_element', gff$type)]
  xelement_df <- data.frame(xelement)
  
  # Adjust strand information
  xelement_df$strand <- as.character(xelement_df$strand)
  xelement_df$strand[is.na(xelement_df$strand)] <- "+"
  xelement_df$strand[xelement_df$strand == "-"] <- "L"
  xelement_df$strand[xelement_df$strand == "+"] <- "R"
  
  chromosomes_strand <- setdiff(xelement_df[, c("seqnames", "strand")], XY_hop1_telo[[2]])
  hop1_telo_df <- data.frame(seqnames = hop1_telo$chr, strand = hop1_telo$arm)
  Sample_telo2 <- merge(hop1_telo, chromosomes_strand, by.x = c("chr", "arm"), by.y = c("seqnames", "strand"))
  
  Sample_telo1 <- Sample_telo2[, 4:ncol(Sample_telo2)]
  Sample_teloGA <- colMeans(Sample_telo1, na.rm = TRUE)
  
  Sample_teloGADF <- data.frame(position = seq(1, length_to_collect), value = Sample_teloGA)
  Sample_teloCompressed <- hwglabr2::compress_signal_track(Sample_teloGADF, window_size = 200)
  Sample_telo <- ksmooth(x = Sample_teloCompressed$position, y = Sample_teloCompressed$window_mean, bandwidth = 25000)
  
  return(list(Sample_telo, chromosomes_strand))
}

# Process X-only ends
XOnly_hop1_telo <- XOnly_Ends_DistanceFromTelo(hop1_telo)
XOnly_red1_telo <- XOnly_Ends_DistanceFromTelo(red1_telo)
XOnly_rec8_telo <- XOnly_Ends_DistanceFromTelo(rec8_telo)

# Plot X-only ends
plot_XOnly <- function(data1, data2, data3) {
  plot(data1[[1]], col = "forestgreen", xlab = "Distance to chr end", ylab = "Average ChIP signal",
       main = "X-only ends", type = "l", bty = "n", lwd = 5, ylim = c(0.4, 1.6), yaxp = c(0.4, 1.6, 6))
  lines(data2[[1]], type = "l", lwd = 5, col = "red")
  lines(data3[[1]], type = "l", lwd = 5, col = "purple")
  abline(h = 1, col = "gray60", lwd = 3, lty = 3)
  legend(50000, 1.79, c("Red1", "Hop1", "Rec8"), text.col = c("forestgreen", "red", "purple"), bty = "n")
}

p2 <- as.ggplot(~plot_XOnly(XOnly_red1_telo, XOnly_hop1_telo, XOnly_rec8_telo))

#Combining the plots, and plotting together
combined_plot <- grid.arrange(p1, p2, ncol = 2)
combined_plot
#ggsave("Fig2A.pdf", plot = combined_plot, width = 10, height = 6, dpi = 300, units = "in")


######################################################################################
# Figure 2B: Anaysis of Fusion chromosomes
######################################################################################

# Function to normalize bedgraph data by genome average
source("R/normalize_signal.R")
# Normalize input datasets
cen1D_Fusion_divided <- gendiv(cen1D_Fusion)  # Normalize cen1D_Fusion data
Red1_WT_divided <- gendiv(Red1_WT)  # Normalize Red1_WT data
cen4D_Fusion_divided <- gendiv(cen4D_Fusion)  # Normalize cen4D_Fusion data

# Function to sort bedgraph data and assign sequence lengths
sort_and_seqlength <- function(genome_name, divided_data) {
  genome_info <- hwglabr2::get_chr_coordinates(genome = genome_name)  # Retrieve genome coordinates
  divided_data <- sort(GenomeInfoDb::sortSeqlevels(divided_data))  # Sort sequences by levels
  genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)  # Sort genome information
  GenomeInfoDb::seqlengths(divided_data) <- GenomeInfoDb::seqlengths(genome_info)  # Assign sequence lengths to the data
  return(divided_data)
}

# Apply sorting and sequence length assignment
cen1D_Fusion_divided <- sort_and_seqlength("SK1_S288CYue", cen1D_Fusion_divided)  # Process cen1D_Fusion
Red1_WT_divided <- sort_and_seqlength("SK1_S288CYue", Red1_WT_divided)  # Process Red1_WT
cen4D_Fusion_divided <- sort_and_seqlength("SK1_S288CYue", cen4D_Fusion_divided)  # Process cen4D_Fusion

# Extract signal from telomeres
length_to_collect <- 160000  # Length of sequence to collect from telomeres
Red1_WT_divided_Telo <- signal_from_telomeres2(Red1_WT_divided, length_to_collect = length_to_collect, averaging_window = 1, genome = "SK1_S288CYue")
cen1D_Fusion_divided_Telo <- signal_from_telomeres2(cen1D_Fusion_divided, length_to_collect = length_to_collect, averaging_window = 1, genome = "SK1_S288CYue")
cen4D_Fusion_divided_Telo <- signal_from_telomeres2(cen4D_Fusion_divided, length_to_collect = length_to_collect, averaging_window = 1, genome = "SK1_S288CYue")

# Function to add L and R labels to chromosomes
add_LR_to_chromosome <- function(divided_data) {
  half <- nrow(divided_data) / 2  # Divide data into two halves
  df_L <- divided_data[1:half, ]  # Left half of the data
  df_R <- divided_data[(half + 1):(2 * half), ]  # Right half of the data
  df_L$chr <- paste0(df_L$chr, "_L")  # Add "_L" suffix to Left half
  df_R$chr <- paste0(df_R$chr, "_R")  # Add "_R" suffix to Right half
  divided_data_LR <- rbind(df_L, df_R)  # Combine Left and Right halves
  divided_data_LR <- divided_data_LR[, !names(divided_data_LR) %in% c("arm", "size_cat")]  # Remove unnecessary columns
  return(divided_data_LR)
}

# Apply L and R labeling
Red1_WT_divided_Telo <- add_LR_to_chromosome(Red1_WT_divided_Telo)
cen1D_Fusion_divided_Telo <- add_LR_to_chromosome(cen1D_Fusion_divided_Telo)
cen4D_Fusion_divided_Telo <- add_LR_to_chromosome(cen4D_Fusion_divided_Telo)


# Function to filter signal for specific chromosomes
FilterForChromosome <- function(signal, chromosomename, length_to_collect) {
  # Filter signal for the specified chromosome
  signal <- signal[signal$chr == chromosomename, , drop = FALSE]  # Subset data based on chromosome name
  signal <- signal[, !names(signal) %in% c("chr")]  # Drop the chromosome name column
  signal <- t(signal)  # Transpose the data
  signal_compressed <- data.frame(position = seq(1, length_to_collect), value = signal)  # Create compressed signal data frame
  signal_compressed <- hwglabr2::compress_signal_track(signal_compressed, window_size = 100)  # Compress signal with a defined window size
  Sample_telo <- ksmooth(x = signal_compressed$position, y = signal_compressed$window_mean, bandwidth = 1000)  # Smooth the signal using kernel smoothing
  return(data.frame(Sample_telo))  # Return the smoothed data as a data frame
}

# Filter signal for specific chromosomes
Red1_WT_divided_Telo_chrIL_value_S288C <- FilterForChromosome(Red1_WT_divided_Telo, "chrI_S288C_L", length_to_collect = length_to_collect)  # Filter chrI Left arm
Fusion_divided_Telo_chrIL_value_S288C <- FilterForChromosome(cen1D_Fusion_divided_Telo, "chrI_S288C_L", length_to_collect = length_to_collect)
Fusion4D_divided_Telo_chrIL_value_S288C <- FilterForChromosome(cen4D_Fusion_divided_Telo, "chrI_S288C_L", length_to_collect = length_to_collect)

# Need only 120 kb of ChrIV, removing the rest by filtering
Red1_WT_divided_Telo_chr4R_value_S288C <- FilterForChromosome(Red1_WT_divided_Telo, "chrIV_S288C_R", length_to_collect = length_to_collect)  # Filter chrIV Right arm
Fusion_divided_Telo_chr4R_value_S288C <- FilterForChromosome(cen1D_Fusion_divided_Telo, "chrIV_S288C_R", length_to_collect = length_to_collect)
Fusion4D_divided_Telo_chr4R_value_S288C <- FilterForChromosome(cen4D_Fusion_divided_Telo, "chrIV_S288C_R", length_to_collect = length_to_collect)
Red1_WT_divided_Telo_chr4R_value_S288C <- Red1_WT_divided_Telo_chr4R_value_S288C[Red1_WT_divided_Telo_chr4R_value_S288C$x < 120000, ]  # Retain data within 120 kb
Fusion_divided_Telo_chr4R_value_S288C <- Fusion_divided_Telo_chr4R_value_S288C[Fusion_divided_Telo_chr4R_value_S288C$x < 120000, ]
Fusion4D_divided_Telo_chr4R_value_S288C <- Fusion4D_divided_Telo_chr4R_value_S288C[Fusion4D_divided_Telo_chr4R_value_S288C$x < 120000, ]


# Function to clean and pivot data for plotting
clean_pivoted_dataframes <- function(df1, df2, df3) {
  df <- cbind(df1, df2[, 2], df3[, 2])  # Combine three data frames
  colnames(df) <- c("Position", "WT", "ChrIV-I(Cen1D)", "ChrIV-I(Cen4D)")  # Rename columns
  df <- df[complete.cases(df[, 2:4]) & apply(df[, 2:4], 1, function(x) all(x >= 0.4)), ]  # Filter valid rows and remove low values
  df_long <- tidyr::pivot_longer(df, cols = c("WT", "ChrIV-I(Cen1D)", "ChrIV-I(Cen4D)"), names_to = "Type", values_to = "Value")  # Transform to long format
  return(df_long[complete.cases(df_long), ])  # Return cleaned long-format data
}

chrI_S288C_L <- clean_pivoted_dataframes(Red1_WT_divided_Telo_chrIL_value_S288C,Fusion_divided_Telo_chrIL_value_S288C,Fusion4D_divided_Telo_chrIL_value_S288C)
chrIV_S288C_R <- clean_pivoted_dataframes(Red1_WT_divided_Telo_chr4R_value_S288C,Fusion_divided_Telo_chr4R_value_S288C, Fusion4D_divided_Telo_chr4R_value_S288C)


span = 1

plot_two_dataframes <- function(df1, df2, chromosomename1, chromosomename2, title) {
  
  df1$Position_flipped <- 120000  - df1$Position
  
  # Create a data frame for the point
  point_df <- data.frame(x = ifelse(title == "S288C", 151967, 154687),
                         y = 0,
                         Type = NA)  # Add a dummy Type column to match aesthetics
  
  
  original_plot <- ggplot(df1, aes(x = Position_flipped, y = Value, color = Type, group = Type)) +
    geom_line(alpha = 0.2, linetype = "dashed", size = 0.3) +  
    geom_smooth(method = "loess", se = FALSE, span = span, size = 1.5) +
    labs(x = chromosomename1, y = "Average ChIP signal (ChIP/Input)", color = "Type") +  
    scale_color_manual(values = c("#724099", "#CB3366","#E5A31F")) +
    scale_x_continuous(breaks = seq(120000 , 0, by = -20000), labels = seq(0, 120000 , by = 20000)) +  
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")) +
    scale_y_continuous(position = "right", limits = c(0, 3), breaks = seq(0, 3, by = 1))+
    ggtitle(title)
  
  point_df <- data.frame(x = ifelse(title == "S288C", 151967, 154687),
                         y = 0,
                         Type = NA)  # Add a dummy Type column to match aesthetics
  
  
  p <- ggplot(df2, aes(x = Position, y = Value, color = Type, group = Type)) +
    geom_line(alpha = 0.2, linetype = "dashed", size = 0.3) +  
    geom_smooth(method = "loess", se = FALSE, span =span, size = 1.5) +  
    labs(x = chromosomename2, y = "Avergae ChIP signal (ChIP/Input)", color = "Type") +  
    scale_color_manual(values = c("#724099", "#CB3366","#E5A31F")) +  
    scale_x_continuous(limits = c(0,160000),breaks = seq(0, 160000 , by = 20000)) +  
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")) +
    scale_y_continuous(limits = c(0, 10) , breaks = seq(0, 10, by = 1)) +
    ggtitle(title) + geom_point(data = point_df, aes(x = x, y = y, color = Type), size = 3) # Create a data frame for the point
  
  
  # Combine both plots into one
  combined_plot <- cowplot::plot_grid(original_plot, p, ncol = 2, align = 'v', axis = 'l', rel_heights = c(1, 1))
  
  return(combined_plot)
}

plot_S288C <- plot_two_dataframes(chrIV_S288C_R, chrI_S288C_L, chromosomename1 = "chrIV-R", chromosomename2 ="chrI-L",  "S288C")
plot_S288C

######################################################################################
# Figure 2C: Coding Distance Versus Distance from Telomeres
######################################################################################

# Load GFF annotations and extract mRNA and Y' elements
gff <- get_gff('SK1Yue')  # Load GFF file for the SK1Yue genome
RNA <- gff[grep('mRNA', gff$type)]  # Extract mRNA annotations
Y_prime_element <- gff[grep("Y_prime_element", gff$type)]  # Extract Y' elements

# Combine mRNA and Y' elements
RNA <- c(RNA, Y_prime_element)

# Get chromosome lengths
chr_length <- get_chr_coordinates("SK1Yue")@seqinfo@seqlengths

# Define bin size for telomere regions
binsize <- 120000

# Create bins for telomere regions for each chromosome
bins <- lapply(seq_along(chr_length), function(i) {
  chr_names <- paste0("chr", as.roman(1:16))  # Chromosome names in Roman numerals
  chr <- chr_names[i]
  length <- chr_length[i]
  
  # Define bins for the start and end of each chromosome
  bin_start <- data.frame(seqnames = chr, start = 1, end = binsize)
  bin_end <- data.frame(seqnames = chr, start = length - binsize + 1, end = length)
  
  # Combine bins for the current chromosome
  rbind(bin_start, bin_end)
})

# Combine bins into a single data frame
bins_df <- do.call(rbind, bins)

# Convert bins to GRanges object
telomere_granges <- GRanges(
  seqnames = bins_df$seqnames,
  ranges = IRanges(start = bins_df$start, end = bins_df$end)
)

# Create sub-bins of 10 kb within the telomere regions
bin_size <- 10000
sub_bins <- unlist(tile(telomere_granges, width = bin_size))

# Find overlaps between sub-bins and RNA elements
overlaps <- findOverlaps(sub_bins, RNA)

# Calculate overlap widths between bins and RNA
overlapping_ranges <- pintersect(sub_bins[queryHits(overlaps)], RNA[subjectHits(overlaps)])
overlap_lengths <- tapply(width(overlapping_ranges), queryHits(overlaps), sum)

# Initialize coding density for all bins
coding_density <- rep(0, length(sub_bins))

# Assign overlap lengths to the corresponding bins
coding_density[as.numeric(names(overlap_lengths))] <- overlap_lengths

# Normalize coding density by bin width
coding_density <- coding_density / width(sub_bins)

# Create a data frame with coding density and bin information
results <- data.frame(
  seqnames = as.character(seqnames(sub_bins)),
  start = start(sub_bins),
  end = end(sub_bins),
  coding_density = coding_density
)

# Rearrange bins to label and flip second halves
rearranged_results <- data.frame()
for (chromosome in unique(results$seqnames)) {
  chr_data <- results[results$seqnames == chromosome, ]
  
  # First 12 bins
  first_part <- chr_data[1:12, ]
  first_part$bin <- paste0("Bin", 1:12)
  
  # Reverse bins after 12 for the second part
  second_part <- chr_data[13:nrow(chr_data), ]
  reversed_bins <- rev(1:nrow(second_part)) + 12
  second_part$bin <- paste0("Bin", reversed_bins)
  
  # Combine and append
  rearranged_results <- rbind(rearranged_results, rbind(first_part, second_part))
}

# Reorder columns for better readability
rearranged_results <- rearranged_results[, c("seqnames", "start", "end", "bin", "coding_density")]

# Function to reshape and flip bins for alternate rows
reshape_and_flip_bins <- function(data, bins_per_row = 12) {
  reshaped <- data %>%
    group_by(seqnames) %>%
    mutate(row_group = ceiling(row_number() / bins_per_row)) %>%
    group_by(seqnames, row_group) %>%
    summarise(bin_values = list(coding_density), .groups = "drop") %>%
    mutate(bin_values = if_else(row_group %% 2 == 0, rev(bin_values), bin_values)) %>%
    ungroup() %>%
    unnest_wider(bin_values, names_sep = "_")
  
  colnames(reshaped) <- c("Chromosome", paste0("Bin", 1:bins_per_row))
  return(reshaped)
}

# Apply the reshape and flip function
flipped_results <- reshape_and_flip_bins(results)

# Remove chromosome columns and calculate mean coding density
flipped_results2 <- flipped_results[, -c(1:2)]
Sample_teloGA <- colMeans(flipped_results2, na.rm = TRUE)
Sample_teloGA <- as.data.frame(Sample_teloGA)
rownames(Sample_teloGA) <- paste0("Bin", 1:12)

# Compute standard deviations for coding density
Sample_teloGA_sds <- apply(flipped_results2, 2, sd, na.rm = TRUE)

# Combine mean and standard deviation into a single data frame
Sample_teloGA_stats <- data.frame(
  Bin = paste0("Bin", 1:12),
  Mean_Coding_Density = Sample_teloGA,
  Std_Dev_Coding_Density = Sample_teloGA_sds
)

# Add distance ranges from the telomere
bin_size <- 10000
Sample_teloGA_stats <- Sample_teloGA_stats %>%
  mutate(
    Distance_from_Telomere = seq(bin_size, bin_size * nrow(Sample_teloGA_stats), by = bin_size)
  )

# Define range labels for the x-axis
range_labels <- paste0(
  seq(0, max(Sample_teloGA_stats$Distance_from_Telomere) - 10000, by = 10000) / 1000, 
  "-", 
  seq(10000, max(Sample_teloGA_stats$Distance_from_Telomere), by = 10000) / 1000, 
  " kb"
)

# Plot Mean Coding Density vs. Distance from Telomere
p <- ggplot(Sample_teloGA_stats, aes(x = Distance_from_Telomere, y = Sample_teloGA)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "darkblue", linetype = "solid", size = 1) +
  labs(
    title = "Mean Coding Density vs. Distance from Telomere (with Y')",
    x = "Distance from Telomere (kb)",
    y = "Mean Coding Density"
  ) +
  scale_x_continuous(
    breaks = Sample_teloGA_stats$Distance_from_Telomere,
    labels = range_labels
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for clarity
  )

p
