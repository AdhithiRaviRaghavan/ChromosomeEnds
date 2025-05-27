################################################################################
# Figure1.R — Full Panel Script (B–F)
# Raghavan et al., bioRxiv 2025. https://doi.org/10.1101/2025.02.27.640173
#
# Author: Adhithi R. Raghavan
# Created: January 2025
#
# Description:
# This script reproduces all panels of Figure 1 from the manuscript:
# "Multiple mechanisms suppress the initiation of meiotic recombination near
# chromosome ends."
#
# Panels:
# - 1A: Schematic (not generated here; illustrated externally)
# - 1B: Average telomeric ChIP signal (Hop1, Red1, Rec8)
# - 1C: Bootstrapped genome-wide signal + average at subtelomeric 20 kb
# - 1D: Meta-signal in X elements with vs. without Y′ ends
# - 1E: Meta-signal in Y′ elements
# - 1F: Barplots comparing signal at X elements with/without Y′ (Cohen's d, t-test)
#
# Required utilities (loaded manually or via `source()`):
# - telo_signal.R
# - bootstrap_distribution.R

################################################################################


#Figure 1: Meiotic axis extends to the telomere-associated sequences

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


# Define local input file paths 
hop1_file <- "data/Hop1_WT.bdg"
red1_file <- "data/Red1_WT.bdg"
rec8_file <- "data/Rec8_WT.bdg"

# Import ChIP signal data
Hop1_WT <- import_bedGraph(hop1_file)
Red1_WT <- import_bedGraph(red1_file)
Rec8_WT <- import_bedGraph(rec8_file)

######################################################################################
# Figure 1A: Schematic representation of the distinct regions present in the chromosome ends 
# Illustrator figure
######################################################################################


######################################################################################
# Figure 1B: Distance from telomeres 
######################################################################################

#This function calculates the distance from the telomeres and produces an output that can be directly plotted.
#The function will set genome average to 1
# Function will also perform spike in normalization, set normalization factor to 1, if you dont want the spike in values.

Red1_WT <- sort(GenomeInfoDb::sortSeqlevels(Red1_WT))

source("R/telo_signal.R")

hop1_telo <- teloSeqSignal(Hop1_WT, normalizationfactor = 1)
red1_telo <- teloSeqSignal(Red1_WT, normalizationfactor = 1)
rec8_telo <- teloSeqSignal(Rec8_WT, normalizationfactor = 1)


myplot <- function(data1,data2,data3) {
  plot(data1[[1]], col="forestgreen", xlab='Distance to chr end', ylab='Average ChIP signal', type='l', bty='n', lwd=5, ylim = c(0.6,1.6),  yaxp = c(0.6, 1.6, 5))
  lines(data2[[1]] , type='l', lwd=5, col="red") 
  lines(data3[[1]], type='l', lwd=5, col="purple") 
  abline(h = 1, col = "gray60", lwd=3, lty=3) 
  legend(50000,1.79, c("Red1", "Hop1", "Rec8"), text.col = c("forestgreen", "red", "purple"), bty = "n")
  
}

library(ggplotify)
p1 <- as.ggplot(~myplot(red1_telo,hop1_telo,rec8_telo ))
p1 

#ggsave("Fig1B.pdf", plot = p1, width = 5, height = 4, dpi = 300, units = "in")

######################################################################################
# Figure 1C: Axis Protein Distribution with Bootstrapping & Mean Signal at 20kb Ends
######################################################################################

# Function: Bootstrapped Distribution of Signal Across the Genome
# This function calculates genome-wide signal distribution using bootstrapping to estimate
# variability in signal intensity across the genome. It performs the following:
# 1. Normalizes signal to the genome average.
# 2. Divides the genome into non-overlapping bins of fixed size (e.g., 20kb).
# 3. Uses bootstrap resampling to calculate mean, median, and confidence intervals (95% CI).

source("R/bootstrap_distribution.R")

### Analyze Hop1, Red1, and Rec8 Signals Using the Function ###
hop1test <- AxisDistributionBootstrapped(Hop1_WT, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32)
red1test <- AxisDistributionBootstrapped(Red1_WT, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32)
rec8test <- AxisDistributionBootstrapped(Rec8_WT, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32)

### Combine Results into Data Frames for Visualization ###
group1s_gg <- data.frame(Data = "Rec8", Score = rec8test[[1]])
group2s_gg <- data.frame(Data = "Hop1", Score = hop1test[[1]])
group3s_gg <- data.frame(Data = "Red1", Score = red1test[[1]])

# Merge all groups
groups_all <- rbind(group1s_gg, group2s_gg, group3s_gg)
groups_all$Strain <- factor(groups_all$Data, levels = c("Rec8", "Hop1", "Red1"))

# Medians and confidence intervals
median_all <- data.frame(
  Strain = c("Rec8", "Hop1", "Red1"),
  Median = c(rec8test[[2]], hop1test[[2]], red1test[[2]])
)
confidenceinterval_all_lower <- data.frame(
  Lower = c(rec8test[[3]][1], hop1test[[3]][1], red1test[[3]][1])
)
confidenceinterval_all_upper <- data.frame(
  Upper = c(rec8test[[3]][2], hop1test[[3]][2], red1test[[3]][2])
)
median_ci_all <- cbind(median_all, confidenceinterval_all_lower, confidenceinterval_all_upper)

### Plotting Signal Distribution with Bootstrap ###
b <- ggplot() +
  geom_violin(data = groups_all, mapping = aes(Strain, Score, fill = Strain)) +
  scale_fill_manual(values = c("#E69F00", "#0072B2", "#000099")) +
  geom_point(mapping = aes(Strain, Median), data = median_ci_all, color = "black") +
  geom_errorbar(mapping = aes(x = Strain, ymin = Lower, ymax = Upper), data = median_ci_all, width = 0.2) +
  ylim(0.6, 1.6) +
  theme_classic() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.5)

### Calculate Mean Signal in Last 20 kb Regions ###
# Define last 20 kb regions
genome_info <- hwglabr2::get_chr_coordinates(genome = "SK1Yue")
seqlength_all <- data.frame(melt(GenomeInfoDb::seqlengths(genome_info)))
seqlength_all <- rownames_to_column(seqlength_all)
colnames(seqlength_all) <- c("Chr", "length")

# Define the last 20 kb regions
binsize <- 20000
coordinates_left <- data.frame()
coordinates_right <- data.frame()
coordinates_all <- data.frame()

for (i in 1:nrow(seqlength_all)) {
  coordinates_left <- cbind(seqlength_all$Chr[i], 0, binsize)
  colnames(coordinates_left) <- c("Chromosome", "Start", "End")
  coordinates_right <- cbind(seqlength_all$Chr[i], seqlength_all$length[i] - binsize, seqlength_all$length[i])
  colnames(coordinates_right) <- c("Chromosome", "Start", "End")
  coordinates_all <- rbind(coordinates_all, coordinates_left, coordinates_right)
}

# Create and import bed file for last 20 kb regions
write.table(coordinates_all, "Last20kb.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
Last20kb_bed <- rtracklayer::import.bed("Last20kb.bed")

# Calculate mean signal in last 20 kb
calculate_mean_last20kb <- function(test_results) {
  signal_matrix <- EnrichedHeatmap::normalizeToMatrix(
    test_results[[4]], Last20kb_bed,
    value_column = "binned_score",
    mean_mode = "weighted", extend = 0, k = 1, empty_value = NA, smooth = FALSE, target_ratio = 1
  )
  return(mean(signal_matrix, na.rm = TRUE))
}



Hop1_mean_last20kb <- calculate_mean_last20kb(hop1test)
Red1_mean_last20kb <- calculate_mean_last20kb(red1test)
Rec8_mean_last20kb <- calculate_mean_last20kb(rec8test)

# Add means to the plot
Last20kb_Means <- data.frame(
  Strain = c("Rec8", "Hop1", "Red1"),
  Mean = c(Rec8_mean_last20kb, Hop1_mean_last20kb, Red1_mean_last20kb)
)

b <- b +
  geom_point(mapping = aes(Strain, Mean), data = Last20kb_Means, color = "#CC0000")

b
# Save and display the plot
#ggsave("Figure1C.pdf", plot = b, width = 5, height = 4, dpi = 300, units = "in")





######################################################################################
# Figure 1D: Meta-Signal Distribution in X Elements with and without Y' Ends
######################################################################################

# Load Genome Information and GFF Data
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')  # Load chromosome info
gff <- get_gff('SK1Yue')  # Load GFF file for SK1Yue genome

# Extract X and Y' elements from GFF data
xelement <- gff[grep('X_element', gff$type)]
yelement <- gff[grep('Y_prime_element', gff$type)]
xelement_df <- data.frame(xelement)
yelement_df <- data.frame(yelement)

######################################################################################
# Part A: X Elements with Y' Ends
######################################################################################

source("R/metaplots.R")
XwithY_gr <- get_XwithY_GRanges("SK1Yue")

# Run metaplot for X elements with Y' ends
Hop1_signalatmetaORF <- MetaplotXwithY(Hop1_WT,1,regions = XwithY_gr )
Red1_signalatmetaORF <- MetaplotXwithY(Red1_WT,1, regions = XwithY_gr)
Rec8_signalatmetaORF <- MetaplotXwithY(Rec8_WT,1, regions = XwithY_gr)

# Combine results for plotting
n_windows = floor(1000/10)
group1s_gg <- data.frame(Data = "Hop1", Position = seq(1, n_windows), Hop1_signalatmetaORF)
group2s_gg <- data.frame(Data = "Red1", Position = seq(1, n_windows), Red1_signalatmetaORF)
group3s_gg <- data.frame(Data = "Rec8", Position = seq(1, n_windows), Rec8_signalatmetaORF)
plotc <- rbind(group1s_gg, group2s_gg, group3s_gg)

# Add upstream (Cen) and downstream (Tel) boundaries
quarter1 <- round(n_windows / 3)
quarter2 <- quarter1 * 2

# Plot X with Y' ends
library(ggplot2)
c_a <- ggplot(plotc, aes(x = Position, y = Mean, group = Data, fill = Data, colour = Data)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line() +
  geom_vline(xintercept = c(quarter1, quarter2), lty = 3) +  # Add Cen/Tel lines
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.5) +
  scale_x_continuous(
    breaks = c(quarter1, quarter2),
    labels = c("Cen", "Tel")  # Update labels to Cen and Tel
  ) +
  ylim(0, 2.8) +
  theme_classic() +
  ggtitle("X with Y' ends")
c_a 

######################################################################################
# Part B: X Elements without Y' Ends
######################################################################################


source("R/metaplots.R")
Xonlygr <- get_Xonly_GRanges("SK1Yue")



# Run metaplot for X-only elements
Hop1_signalatmetaORF_Xonly <- MetaplotXonly(Hop1_WT,1, regions = Xonlygr)
Red1_signalatmetaORF_Xonly <- MetaplotXonly(Red1_WT,1, regions = Xonlygr)
Rec8_signalatmetaORF_Xonly <- MetaplotXonly(Rec8_WT,1, regions = Xonlygr)

# Combine results for plotting
group1s_gg_xonly <- data.frame(Data = "Hop1", Position = seq(1, n_windows), Hop1_signalatmetaORF_Xonly)
group2s_gg_xonly <- data.frame(Data = "Red1", Position = seq(1, n_windows), Red1_signalatmetaORF_Xonly)
group3s_gg_xonly <- data.frame(Data = "Rec8", Position = seq(1, n_windows), Rec8_signalatmetaORF_Xonly)
plotc_b <- rbind(group1s_gg_xonly, group2s_gg_xonly, group3s_gg_xonly)

# Plot X-only elements
c_b <- ggplot(plotc_b, aes(x = Position, y = Mean, group = Data, fill = Data, colour = Data)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line() +
  geom_vline(xintercept = c(quarter1, quarter2), lty = 3) +  # Add Cen/Tel lines
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.5) +
  scale_x_continuous(
    breaks = c(quarter1, quarter2),
    labels = c("Cen", "Tel")  # Update labels to Cen and Tel
  ) +
  ylim(0, 2.8) +
  theme_classic() +
  ggtitle("X only ends")
c_b 

######################################################################################
# Combine Plots
######################################################################################

# Combine the plots side by side
library(gridExtra)
c <- grid.arrange(c_a, c_b, ncol = 2)

# Save the combined plot
#ggsave("Figure1D.pdf", plot = c, width = 10, height = 6, dpi = 300, units = "in")

######################################################################################
# Figure 1E: Meta-Signal Distribution in Y' Elements
######################################################################################
source("R/metaplots.R")
Yprimegr <- get_Yprime_GRanges("SK1Yue")

# Run the function for different datasets
Hop1_signalatmetaORF_Y <- MetaplotY(Hop1_WT)
Red1_signalatmetaORF_Y <- MetaplotY(Red1_WT)
Rec8_signalatmetaORF_Y <- MetaplotY(Rec8_WT)

# Prepare data for plotting
n_windows <- floor(1000 / 10)  # Define number of windows
group1s_gg_Y <- data.frame(Data = "Hop1", Position = seq(1, n_windows), Hop1_signalatmetaORF_Y)
group2s_gg_Y <- data.frame(Data = "Red1", Position = seq(1, n_windows), Red1_signalatmetaORF_Y)
group3s_gg_Y <- data.frame(Data = "Rec8", Position = seq(1, n_windows), Rec8_signalatmetaORF_Y)

# Combine all datasets into one for plotting
plotd <- rbind(group1s_gg_Y, group2s_gg_Y, group3s_gg_Y)



# Plot the data
library(ggplot2)
d <- ggplot(plotd, aes(x = Position, y = Mean, group = Data, fill = Data, colour = Data)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line() +
  geom_vline(xintercept = c(25, 75), lty = 3) +  # Vertical lines for upstream and downstream regions
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.5) +
  scale_x_continuous(
    breaks = c(25, 75),
    labels = c("Cen", "Tel")
  ) +
  ylim(0, 2.8) +
  ylab("Mean") +
  ggtitle("Y' element") +
  theme_classic() +
  theme(legend.position = "bottom")

# Display the plot
d

#ggsave("Figure1E.pdf",plot = d, width = 5, height = 4, dpi = 300, units = "in")

######################################################################################
# Figure 1F: Comparison of Signal Between X Elements with and without Y' Ends
######################################################################################

# Reuse `XwithYrows_all_Grange` and `Xonly_all_Grange` from earlier script

# Normalize to Genome Average and Calculate Mean Signal for X Elements
## This function calculates the normalized signal and generates a data frame of means per X element.
MetaplotX_MeanPerX <- function(signal, regions) {
  # Normalize signal to genome average
  gendiv <- function(bdg) {
    gavg <- hwglabr2::average_chr_signal(bdg)$genome_avrg
    bdg$score <- bdg$score / gavg
    return(bdg)
  }
  signal_divided <- gendiv(signal)
  
  # Get mean signal for specified regions
  signal_matrix <- EnrichedHeatmap::normalizeToMatrix(
    signal_divided, regions,
    value_column = "score", mean_mode = "absolute",
    extend = 0, k = 1, empty_value = NA, smooth = FALSE, target_ratio = 1
  )
  return(data.frame(signal_matrix))
}

# Run the function for X with Y' ends and X-only ends
Hop1_signal_XwithY <- MetaplotX_MeanPerX(Hop1_WT, XwithYrows_all_Grange)
Red1_signal_XwithY <- MetaplotX_MeanPerX(Red1_WT, XwithYrows_all_Grange)
Rec8_signal_XwithY <- MetaplotX_MeanPerX(Rec8_WT, XwithYrows_all_Grange)

Hop1_signal_Xonly <- MetaplotX_MeanPerX(Hop1_WT, Xonly_all_Grange)
Red1_signal_Xonly <- MetaplotX_MeanPerX(Red1_WT, Xonly_all_Grange)
Rec8_signal_Xonly <- MetaplotX_MeanPerX(Rec8_WT, Xonly_all_Grange)

# Function to perform t-test and return p-value
perform_t_test <- function(df1, df2) {
  combined_data <- bind_rows(
    mutate(df1, Group = "Group 1"),
    mutate(df2, Group = "Group 2")
  )
  
  # Perform t-test and return p-value
  t_test_result <- t.test(t1 ~ Group, data = combined_data)
  return(t_test_result$p.value)
}

# Apply t-test for all comparisons and get p-values
hop1_pvalue <- perform_t_test(Hop1_signal_XwithY, Hop1_signal_Xonly)
red1_pvalue <- perform_t_test(Red1_signal_XwithY, Red1_signal_Xonly)
rec8_pvalue <- perform_t_test(Rec8_signal_XwithY, Rec8_signal_Xonly)

# Collect all p-values into a vector
p_values <- c(hop1_pvalue, red1_pvalue, rec8_pvalue)

# Apply Benjamini-Hochberg correction
adjusted_pvalues <- p.adjust(p_values, method = "BH")

# Create Bar Plots for Mean Signal Comparison
## This function creates bar plots with means, error bars, and significance stars.
makebarplots <- function(df1, df2, plot_title = "Comparison of Means", custom_colors = NULL, corrected_pvalue) {   
  combined_data <- bind_rows(     
    mutate(df1, Group = "Group 1"),     
    mutate(df2, Group = "Group 2")   
  )     
  
  # Calculate means and standard deviations   
  group_summary <- combined_data %>%     
    group_by(Group) %>%     
    summarise(mean_t1 = mean(t1), sd_t1 = sd(t1))     
  
  # Create the bar plot   
  bar_plot <- ggplot(group_summary, aes(x = Group, y = mean_t1)) +     
    geom_bar(stat = "identity", position = "dodge", fill = custom_colors[1]) +     
    geom_errorbar(aes(ymin = mean_t1 - sd_t1, ymax = mean_t1 + sd_t1),                   
                  position = position_dodge(0.9), width = 0.2, color = "black", linewidth = 0.7) +     
    labs(title = plot_title, x = "Group", y = "Mean Value") +     
    theme_classic() +
    scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5)) # Set Y-axis limits
  
  # Add stars based on adjusted p-value   
  significance_label <- ifelse(corrected_pvalue <= 0.001, "***",                                
                               ifelse(corrected_pvalue <= 0.01, "**",                                       
                                      ifelse(corrected_pvalue <= 0.05, "*", "n.s.")))      
  
  # Define a fixed Y position for the significance label
  fixed_y_position <- 1.75  # Adjust this value based on your specific data
  
  # Add significance stars to the plot   
  bar_plot <- bar_plot +      
    annotate("text", x = 1.5, y = fixed_y_position, label = significance_label, size = 8, vjust = -0.5)  # Adjust vjust if needed
  
  return(bar_plot) 
}

# Generate plots using adjusted p-values
plot1 <-makebarplots(Hop1_signal_Xonly,Hop1_signal_XwithY, plot_title = "Hop1", custom_colors = c("#80D1EB", "#80D1EB"), corrected_pvalue = adjusted_pvalues[1])
plot2 <- makebarplots(Red1_signal_Xonly, Red1_signal_XwithY, plot_title = "Red1", custom_colors = c("#FFA07A", "#FFA07A"), corrected_pvalue = adjusted_pvalues[2])
plot3 <- makebarplots(Rec8_signal_Xonly, Rec8_signal_XwithY, plot_title = "Rec8", custom_colors = c("#00FA9A", "#00FA9A"), corrected_pvalue = adjusted_pvalues[3])

library(patchwork)
# Combine plots
combined_plot <- plot3 + plot1 + plot2
combined_plot <- combined_plot + plot_layout(ncol = 3)
print(combined_plot)

library(effsize)
cohen_d1 <- cohen.d(Rec8_signal_Xonly$t1, Rec8_signal_XwithY$t1)
cohen_d1 
cohen_d2 <- cohen.d(Hop1_signal_Xonly$t1, Hop1_signal_XwithY$t1)
cohen_d2
cohen_d3 <- cohen.d(Red1_signal_Xonly$t1, Hop1_signal_XwithY$t1)
cohen_d3
cohen_d <- c(cohen_d1, cohen_d2 ,cohen_d3)

cohen_d
#ggsave("Figure1F.pdf",plot = combined_plot, width = 5, height = 4, dpi = 300, units = "in")