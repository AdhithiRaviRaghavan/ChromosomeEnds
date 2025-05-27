################################################################################
# Figure3.R — Full Panel Script (A–D)
# Raghavan et al., bioRxiv 2025 — https://doi.org/10.1101/2025.02.27.640173
#
# Author: Adhithi R. Raghavan
# Created: January 2025
#
# Description:
# This script reproduces all panels for Figure 3 from the manuscript:
# "Multiple mechanisms suppress the initiation of meiotic recombination near chromosome ends."
# It focuses on Red1 ChIP-seq signal changes across mutants, normalized using SNP-ChIP spike-in.
#
# Panels:
# - 3A: Red1 ChIP signal as a function of distance from telomeres
# - 3B: Genome-wide bootstrapping of Red1 ChIP signal
# - 3C: Metaplot of Red1 signal at X elements (with/without Y′ ends)
# - 3D: Metaplot of Red1 signal at Y′ elements
#
# Utility functions required (must be sourced separately):
# - R/metaplots.R
# - R/normalize_signal.R
#
# Notes:
# - This script uses spike-in normalization factors calculated via `hwglabr2::spikein_normalization_factor`.
# - This is a working research script and may evolve as the manuscript undergoes review.
################################################################################

#Figure 3:Differential recruitment of Red1 at chromosome ends by Rec8-dependent and independent pathways.

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

######################################
# Define local input file paths
######################################

red1_file <- "data/Red1_WT.bdg"
red1_rec8_file <- "data/Red1_Rec8.bdg"
red1_phd_file <- "data/Red1_Phd.bdg"
red1_phd_rec8_file <- "data/Red1_Phd_Rec8.bdg"

######################################
# Import ChIP-seq signal
######################################

Red1_WT <- import_bedGraph(red1_file)
Red1_Rec8 <- import_bedGraph(red1_rec8_file)
Red1_Phd <- import_bedGraph(red1_phd_file)
Red1_Phd_Rec8 <- import_bedGraph(red1_phd_rec8_file)

######################################
# Spike-in Normalization Factors
######################################

ref_input_7797 <- "spikein/ref_input_7797.txt"
ref_chip_7797 <- "spikein/ref_chip_7797.txt"
test_input_5187 <- "spikein/test_input_5187.txt"
test_chip_5187 <- "spikein/test_chip_5187.txt"
test_input_9251 <- "spikein/test_input_9251.txt"
test_chip_9251 <- "spikein/test_chip_9251.txt"
test_input_10517 <- "spikein/test_input_10517.txt"
test_chip_10517 <- "spikein/test_chip_10517.txt"

#Calculating Spike in:
normalizationfactor <- hwglabr2::spikein_normalization_factor(
  ref_input_counts = ref_input_counts,
  ref_chip_counts = ref_chip_counts,
  test_input_counts = test_input_5187,
  test_chip_counts = test_chip_5187
)

normalizationfactor_phd <- hwglabr2::spikein_normalization_factor(
  ref_input_counts = ref_input_counts,
  ref_chip_counts = ref_chip_counts,
  test_input_counts = test_input_9251,
  test_chip_counts = test_chip_9251
)

normalizationfactor_phd_rec8 <- hwglabr2::spikein_normalization_factor(
  ref_input_counts = ref_input_counts,
  ref_chip_counts = ref_chip_counts,
  test_input_counts = test_input_10517,
  test_chip_counts = test_chip_10517
)

################################################################################
# Figure 3A: Red1 ChIP signal as a function of distance from telomeres
################################################################################

# This function calculates average ChIP signal as a function of distance from telomeres.
# It sets genome-wide average to 1 and applies spike-in normalization.
# To disable spike-in, set `normalizationfactor = 1`.

source("R/telo_signal.R")
Red1_Rec8_telo <- teloSeqSignal(Red1_Rec8, normalizationfactor = normalizationfactor)
Red1_Phd_telo <- teloSeqSignal(Red1_Phd, normalizationfactor = normalizationfactor_phd)
Red1_Phd_Rec8_telo <- teloSeqSignal(Red1_Phd_Rec8, normalizationfactor = normalizationfactor_phd_rec8)

myplot <- function(data1,data2,data3,data4) {
  plot(data1[[1]],   col="#CC3366", xlab='Distance to chr end', ylab='Average ChIP signal', type='l', bty='n', lwd=5,ylim = c(0.2,1.2),  yaxp = c(0.2, 1.2, 5), xaxp = c(0,200000,10))
  lines(data2[[1]] , type='l', lwd=5, col="#FFB000") 
  lines(data3[[1]] , type='l', lwd=5, col="#619CFF") 
  lines(data4[[1]] , type='l', lwd=5, col="#6900Ae") 
  abline(h = data1[[2]], col = "#CC3366", lwd=3, lty=3) 
  abline(h = data2[[2]], col = "#FFB000", lwd=3, lty=3) 
  abline(h = data3[[2]], col = "#619CFF", lwd=3, lty=3) 
  abline(h = data4[[2]], col = "#6900Ae", lwd=3, lty=3) 
  
  legend(80000,0.6, c("Red1",  "Red1-Rec8Delta", "Red1-PHD", "Red1-PHDRec8"), text.col = c("#CC3366", "#FFB000", "#619CFF","#6900Ae"), bty = "n")
  
}

library(ggplotify)
p1 <- as.ggplot(~myplot(red1_telo,Red1_Rec8_telo,Red1_Phd_telo,Red1_Phd_Rec8_telo))
p1 

#ggsave("Fig3A.pdf", plot = p1, width = 5, height = 4, dpi = 300, units = "in")


##############################################################################################################
# Figure 3B:Axis factors bootstrapped
###############################################################################################################

# Function: Bootstrapped Distribution of Signal Across the Genome
# This function calculates genome-wide signal distribution using bootstrapping to estimate
# variability in signal intensity across the genome. It performs the following:
# 1. Normalizes signal to the genome average.
# 2. Divides the genome into non-overlapping bins of fixed size (e.g., 20kb).
# 3. Uses bootstrap resampling to calculate mean, median, and confidence intervals (95% CI).

source("R/bootstrap_distribution.R")

### Analyze Hop1, Red1, and Rec8 Signals Using the Function ###
red1test <- AxisDistributionBootstrapped(signal = Red1_WT, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32, normalizationfactor=1)
red1_rec8test <- AxisDistributionBootstrapped(signal = Red1_Rec8, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32, normalizationfactor = normalizationfactor)
red1_phdtest <- AxisDistributionBootstrapped(signal = Red1_Phd, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32, normalizationfactor = normalizationfactor_phd)
red1_phdrec8test <- AxisDistributionBootstrapped(signal = Red1_Phd_Rec8, tilesize = 20000, bootstrappedNo = 1000, samplePerBootstrap = 32, normalizationfactor = normalizationfactor_phd_rec8)

### Combine Results into Data Frames for Visualization ###
group1s_gg <- data.frame(Data="Red1",red1test[[1]])
colnames(group1s_gg) <- c("Strain", "Score")

group2s_gg <- data.frame(Data="Red1_Rec8",red1_rec8test[[1]])
colnames(group2s_gg) <- c("Strain", "Score")

group3s_gg <- data.frame(Data="Red1_Phd",red1_phdtest[[1]])
colnames(group3s_gg) <- c("Strain", "Score")

group4s_gg <- data.frame(Data="Red1_Phd_Rec8",red1_phdrec8test[[1]])
colnames(group4s_gg) <- c("Strain", "Score")

# Merge all groups
groups_all = rbind(group1s_gg,group2s_gg,group3s_gg,group4s_gg)
groups_all$Strain <- factor(groups_all$Strain, levels = c( "Red1", "Red1_Rec8","Red1_Phd","Red1_Phd_Rec8"))

# Medians and confidence intervals
median_all <-data.frame(c(red1test[[2]],red1_rec8test[[2]],red1_phdtest[[2]],red1_phdrec8test[[2]]))
row.names(median_all) <- c( "Red1", "Red1_Rec8","Red1_Phd","Red1_Phd_Rec8")
median_all <- rownames_to_column(median_all)
colnames(median_all) <- c("Strain", "Median")
median_all$Strain <- factor(median_all$Strain, levels = c("Red1", "Red1_Rec8","Red1_Phd","Red1_Phd_Rec8"))

confidenceinterval_all_lower <-data.frame(c(red1test[[3]][1],red1_rec8test[[3]][1],red1_phdtest[[3]][1],red1_phdrec8test[[3]][1] ))
confidenceinterval_all_upper <-data.frame(c(red1test[[3]][2],red1_rec8test[[3]][2],red1_phdtest[[3]][2],red1_phdrec8test[[3]][2] ))
median_ci_all <-cbind(median_all, confidenceinterval_all_lower,confidenceinterval_all_upper)
colnames(median_ci_all) <- c("Strain", "Median", "Lower", "Upper")

### Calculate Mean Signal in Last 20 kb Regions ###
# Define last 20 kb regions
genome_info <- hwglabr2::get_chr_coordinates(genome = "SK1Yue")
seqlength_all <- data.frame(melt(GenomeInfoDb::seqlengths(genome_info)))
seqlength_all <- rownames_to_column(seqlength_all)
colnames(seqlength_all) <- c("Chr", "length")

#Getting the list of the last 20 kb
spreadingfromcoordinate = data.frame()
dataframe_left = data.frame()

binsize = 20000
coordinates_left = data.frame()
coordinates_right = data.frame()
coordinates_all = data.frame()
for(i in 1:nrow(seqlength_all)){
  coordinates_left <- cbind(seqlength_all$Chr[i],0, binsize)
  colnames( coordinates_left) <- c("Chromosome", "Start", "End")
  coordinates_right <- cbind(seqlength_all$Chr[i],seqlength_all$length[i] - binsize,seqlength_all$length[i] )
  colnames( coordinates_right) <- c("Chromosome", "Start", "End")
  coordinates_all <- rbind(coordinates_all, coordinates_left,coordinates_right )
}


#Making a bed file
write.table(coordinates_all,'Last20kb.bed',quote = F,row.names = F,col.names = F,sep='\t')
Last20kb_bed = rtracklayer::import.bed('Last20kb.bed')

# Calculate mean signal in last 20 kb
calculate_mean_last20kb <- function(test_results) {
  signal_matrix <- EnrichedHeatmap::normalizeToMatrix(
    test_results[[4]], Last20kb_bed,
    value_column = "binned_score",
    mean_mode = "weighted", extend = 0, k = 1, empty_value = NA, smooth = FALSE, target_ratio = 1
  )
  return(mean(signal_matrix, na.rm = TRUE))
}

Red1_WT_last20kbmean <- calculate_mean_last20kb(red1test)
Red1_Rec8_last20kbmean <- calculate_mean_last20kb(red1_rec8test)
Red1_Phd_last20kbmean  <- calculate_mean_last20kb(red1_phdtest)
Red1_Phd_Rec8_last20kbmean   <- calculate_mean_last20kb(red1_phdrec8test)

Last20kb_Means <- data.frame(
  Strain = c("Red1", "Red1_Rec8","Red1_Phd","Red1_Phd_Rec8"),
  Mean = c(Red1_WT_last20kbmean,Red1_Rec8_last20kbmean,Red1_Phd_last20kbmean,Red1_Phd_Rec8_last20kbmean )
)

                       
####Plotting###
b <- ggplot() + 
  geom_violin(data = groups_all[groups_all$Strain != "Red1_Phd_Rec8", ], 
              mapping = aes(Strain, Score, fill = Strain), 
              width = 0.8, position = position_dodge(width = 0.5), adjust = 1) +
  geom_violin(data = groups_all[groups_all$Strain == "Red1_Phd_Rec8", ], 
              mapping = aes(Strain, Score, fill = Strain), 
              width = 0.8, position = position_dodge(width = 0.5), adjust = 0.5) +
  scale_fill_manual(values = c("#CC3366", "#FFB000", "#619CFF", "#6900Ae")) +
  geom_point(mapping = aes(Strain, Median), data = median_ci_all, position = position_dodge(width = 0.8)) +
  geom_point(mapping = aes(Strain, Mean), data = Last20kb_Means, colour = "#CC0000", position = position_dodge(width = 0.8)) +
  geom_errorbar(mapping = aes(x = Strain, ymin = Lower, ymax = Upper), 
                data = median_ci_all, width = 0.2, position = position_dodge(width = 0.8)) +
  scale_y_continuous(limits = c(0.15, 1.3), breaks = seq(0, 1.4, by = 0.2)) +
  theme_classic() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#CC3366", size = 0.5) +
  geom_hline(yintercept = normalizationfactor, linetype = "dashed", color = "#FFB000", size = 0.5) +
  geom_hline(yintercept = normalizationfactor_phd, linetype = "dashed", color = "#619CFF", size = 0.5) +
  geom_hline(yintercept = normalizationfactor_phd_rec8, linetype = "dashed", color = "#6900Ae", size = 0.5) +
  guides(fill = "none")

print(b)

#ggsave("Fig3B.pdf", plot = p1, width = 5, height = 4, dpi = 300, units = "in")

######################################################################################
# Figure 3C: Meta-Signal Distribution in X Elements with and without Y' Ends
######################################################################################


source("R/metaplots.R")
XwithY_gr <- get_XwithY_GRanges("SK1Yue")

# Calculate meta-signals
Red1_signal <- MetaplotXwithY(Red1_WT, normalizationfactor = 1,regions = XwithY_gr )
Red1_Rec8_signal <- MetaplotXwithY(Red1_Rec8, normalizationfactor,regions = XwithY_gr )
Red1_Phd_signal <- MetaplotXwithY(Red1_Phd, normalizationfactor_phd,regions = XwithY_gr )
Red1_PhdRec8_signal <- MetaplotXwithY(Red1_Phd_Rec8, normalizationfactor_phd_rec8, regions = XwithY_gr )

# Combine data for plotting
n_windows <- floor(1000 / 10)
plot_a <- rbind(
  data.frame(Data = "Red1", Position = seq(1, n_windows), Red1_signal),
  data.frame(Data = "Red1-Rec8Delta", Position = seq(1, n_windows), Red1_Rec8_signal),
  data.frame(Data = "Red1-PHD", Position = seq(1, n_windows), Red1_Phd_signal),
  data.frame(Data = "Red1-PHDRec8", Position = seq(1, n_windows), Red1_PhdRec8_signal)
)

c_a <- ggplot(allxgroups, aes(x = Position, y = Mean, group = Data, fill = Data, colour = Data)) +
  theme_minimal() +
  geom_vline(xintercept = c(33, 66), lty = 3) +
  scale_x_continuous(breaks = c(33, 66), labels = c('Cen', 'Tel')) +
  ylim(0, 2.8) +
  ggtitle("X with Y' ends") +
   geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Data), alpha = 0.3, color = NA) +
  geom_line(aes(color = Data), size = 1) +
  scale_fill_manual(values = c("Red1" = "#CC3366", "Red1-Rec8Delta" = "#FFB000", "Red1-PHD" = "#619CFF","Red1-PHDRec8"="#6900Ae")) +
  scale_color_manual(values = c("Red1" = "#CC3366", "Red1-Rec8Delta" = "#FFB000", "Red1-PHD" = "#619CFF","Red1-PHDRec8"="#6900Ae"))+ theme_classic()

c_a

######################################################################################
# Part B: X Elements without Y' Ends
######################################################################################

source("R/metaplots.R")
Xonlygr <- get_Xonly_GRanges("SK1Yue")

#Runing the function
Red1_signalatmetaORF_Xonly <- MetaplotXonly(Red1_WT, 1, regions = Xonlygr)
Red1_Rec8_signalatmetaORF_Xonly <- MetaplotXonly(Red1_Rec8, normalizationfactor,regions = Xonlygr )
Red1_Phd_signalatmetaORF_Xonly <- MetaplotXonly(Red1_Phd, normalizationfactor_phd, regions = Xonlygr)
Red1_Phd_Rec8_signalatmetaORF_Xonly <- MetaplotXonly(Red1_Phd_Rec8, normalizationfactor_phd_rec8, regions = Xonlygr)

# Combine data for plotting
plot_b <- rbind(
  data.frame(Data = "Red1", Position = seq(1, n_windows),Red1_signalatmetaORF_Xonly),
  data.frame(Data = "Red1-Rec8Delta", Position = seq(1, n_windows),Red1_Rec8_signalatmetaORF_Xonly),
  data.frame(Data = "Red1-PHD", Position = seq(1, n_windows),Red1_Phd_signalatmetaORF_Xonly),
  data.frame(Data = "Red1-PHDRec8", Position = seq(1, n_windows),Red1_Phd_Rec8_signalatmetaORF_Xonly)
)

c_b <- ggplot(plot_b, aes(x = Position, y = Mean, group = Data, fill = Data, colour = Data)) +
  theme_minimal() +
  geom_vline(xintercept = c(33, 66), lty = 3) +
  scale_x_continuous(breaks = c(33, 66), labels = c('Cen', 'Tel')) +
  ylim(0, 2.8) +
 
  ggtitle("X only ends") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Data), alpha = 0.3, color = NA) +
  geom_line(aes(color = Data), size = 1) +
  scale_fill_manual(values = c("Red1" = "#CC3366", "Red1-Rec8Delta" = "#FFB000", "Red1-PHD" = "#619CFF","Red1-PHDRec8"="#6900Ae")) +
  scale_color_manual(values = c("Red1" = "#CC3366", "Red1-Rec8Delta" = "#FFB000", "Red1-PHD" = "#619CFF","Red1-PHDRec8"="#6900Ae")) +
  theme(legend.position = "none") +
  theme_classic()

c_b

library(gridExtra)
x_all <- grid.arrange(
  c_b, c_a, ncol = 2,
  widths = c(4.08, 4.08),
  heights = unit(5.5, "cm")  )

#ggsave("Figure3C.pdf", plot = x_all, width = 10, height = 6, dpi = 300, units = "in")

######################################################################################
# Figure 3D: Meta-Signal Distribution in Y' Elements
######################################################################################
source("R/metaplots.R")
Yprimegr <- get_Yprime_GRanges("SK1Yue")

#Runing the function
Red1_signalatmetaORF_Y <- MetaplotY(Red1_WT,1, regions = Yprimegr)
Red1_Rec8_signalatmetaORF_Y <-MetaplotY(Red1_Rec8, normalizationfactor, regions = Yprimegr)
Red1_Phd_signalatmetaORF_Y <- MetaplotXonly(Red1_Phd, normalizationfactor_phd, regions = Yprimegr)
Red1_Phd_Rec8_signalatmetaORF_Y <- MetaplotXonly(Red1_Phd_Rec8, normalizationfactor_phd_rec8, regions = Yprimegr)

# Prepare data for plotting
n_windows = floor(1000/10)
plot_b <- rbind(
  data.frame(Data = "Red1", Position = seq(1, n_windows),Red1_signalatmetaORF_Y),
  data.frame(Data = "Red1-Rec8Delta", Position = seq(1, n_windows),Red1_Rec8_signalatmetaORF_Y),
  data.frame(Data = "Red1-PHD", Position = seq(1, n_windows),Red1_Phd_signalatmetaORF_Y),
  data.frame(Data = "Red1-PHDRec8", Position = seq(1, n_windows),Red1_Phd_Rec8_signalatmetaORF_Y)
)


d <- ggplot(plot_b, aes(x = Position, y = Mean, group = Data, fill = Data, colour = Data)) +
  theme_minimal() +
  geom_vline(xintercept = c(25,75), lty = 3) +
  scale_x_continuous(breaks = c(25,75), labels = c('upstream', 'downstream')) +
  ylim(0, 2.8) +
  ggtitle("Y' element") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Data), alpha = 0.3, color = NA) +
  geom_line(aes(color = Data), size = 1) +
  scale_fill_manual(values = c("Red1" = "#CC3366", "Red1-Rec8Delta" = "#FFB000", "Red1-PHD" = "#619CFF","Red1-PHDRec8"="#6900Ae")) +
  scale_color_manual(values = c("Red1" = "#CC3366", "Red1-Rec8Delta" = "#FFB000", "Red1-PHD" = "#619CFF","Red1-PHDRec8"="#6900Ae"))

d
d<- d + 
  theme_classic()+ theme(legend.position="none")
d


y_all <- grid.arrange(
  d, ncol = 1,
  widths = c(4.08),
  heights = unit(5.5, "cm") )
#ggsave("Figure3D.pdf", plot = Y_all, width = 10, height = 6, dpi = 300, units = "in")


