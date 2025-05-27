#' Calculate Signal from Telomeres
#'
#' This function calculates the signal distribution from chromosome ends,
#' normalizing to genome average and applying optional spike-in normalization.
#'
#' @param Sample A GenomicRanges object.
#' @param remove_cen Logical. Remove centromere regions?
#' @param length_to_collect How many bp from each telomere end to collect.
#' @param Genome Genome assembly. Default = "SK1Yue".
#' @param normalizationfactor Optional normalization multiplier.
#' @return A list containing smoothed telomeric signal and genome average.
teloSeqSignal <- function(Sample, remove_cen = FALSE, length_to_collect = 120000, Genome = 'SK1Yue', normalizationfactor = 1) {
  library(hwglabr2)
  library(dplyr)
  library(magrittr)
  library(pbapply)
  library(readr)
  
  if (Genome == 'SK1Yue') {
    cen <- hwglabr2::get_chr_coordinates('SK1Yue')  
  }
  
  if (remove_cen) {
    half_length <- floor(cen_region_length / 2)
    offset <- floor(GenomicRanges::width(cen) / 2)
    GenomicRanges::start(cen) <- (GenomicRanges::start(cen) + offset - half_length)
    GenomicRanges::end(cen) <- (GenomicRanges::end(cen) - offset + half_length)
    Sample <- Sample[!IRanges::overlapsAny(Sample, cen)]
  }
  
  Sample_genomeAvg <- hwglabr2::average_chr_signal(Sample)[[2]]
  GenomicRanges::score(Sample) <- GenomicRanges::score(Sample) / Sample_genomeAvg
  GenomicRanges::score(Sample) <- GenomicRanges::score(Sample) * normalizationfactor
  
  Sample_telo <- hwglabr2::signal_from_telomeres2(Sample, length_to_collect, genome = Genome)
  Sample_telo1 <- Sample_telo[, 4:ncol(Sample_telo)]
  Sample_teloGA <- colMeans(Sample_telo1, na.rm = TRUE)
  
  Sample_teloGADF <- data.frame(position = seq(1, length_to_collect), value = Sample_teloGA)
  Sample_teloCompressed <- hwglabr2::compress_signal_track(Sample_teloGADF, window_size = 200)
  Sample_teloSmooth <- ksmooth(x = Sample_teloCompressed$position, y = Sample_teloCompressed$window_mean, bandwidth = 25000)
  
  list(Sample_teloSmooth, hwglabr2::average_chr_signal(Sample)[[2]])
}
