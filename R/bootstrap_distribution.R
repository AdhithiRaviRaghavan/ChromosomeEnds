#' Bootstrapped Distribution of Signal Across Genome
#'
#' Estimate signal variability by bootstrapping across genome bins.
#'
#' @param signal A GenomicRanges object with signal.
#' @param genome Genome assembly.
#' @param tilesize Tile size in bp.
#' @param bootstrappedNo Number of bootstrap replicates.
#' @param samplePerBootstrap Number of bins per replicate.
#' @param replacement Sample with replacement?
#' @param normalizationfactor Spike-in normalization factor.
#' @return List of bootstrap results, median, CI, and binned signal.
AxisDistributionBootstrapped <- function(signal, genome = 'SK1Yue', tilesize, bootstrappedNo = 1000, samplePerBootstrap, replacement = TRUE, normalizationfactor = 1) {
  signal_divided <- gendiv(signal)
  GenomicRanges::score(signal_divided) <- GenomicRanges::score(signal_divided) * normalizationfactor
  genome_info <- hwglabr2::get_chr_coordinates(genome = genome)
  signal_divided_sort <- sort(GenomeInfoDb::sortSeqlevels(signal_divided))
  GenomeInfoDb::seqlengths(signal_divided_sort) <- GenomeInfoDb::seqlengths(genome_info)
  
  bins_signal_divided <- GenomicRanges::tileGenome(
    GenomeInfoDb::seqlengths(signal_divided_sort),
    tilewidth = tilesize,
    cut.last.tile.in.chrom = TRUE
  )
  
  score_signal_divided <- GenomicRanges::coverage(signal_divided_sort, weight = "score")
  bins_signal_divided <- GenomicRanges::binnedAverage(bins_signal_divided, score_signal_divided, "binned_score")
  matrix_signal <- matrix(bins_signal_divided$binned_score)
  
  bootstrap <- replicate(bootstrappedNo, sample(x = matrix_signal, samplePerBootstrap, replace = replacement), simplify = FALSE)
  Meanof32_bootstrap <- sapply(bootstrap, mean, simplify = TRUE)
  median_bootstrap <- median(Meanof32_bootstrap)
  confidenceinterval_bootstrap <- matrix(quantile(Meanof32_bootstrap, c(0.025, 0.975)))
  
  list(Meanof32_bootstrap, median_bootstrap, confidenceinterval_bootstrap, bins_signal_divided)
}
