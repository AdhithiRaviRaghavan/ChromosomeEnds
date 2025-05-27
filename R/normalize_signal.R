#' Normalize Signal to Genome Average
#'
#' This function normalizes signal scores in a GenomicRanges object by dividing 
#' them by the genome-wide average signal.
#'
#' @param bdg A GenomicRanges object containing signal scores.
#' @return A GenomicRanges object with normalized signal scores.
#' @examples
#' normalized_signal <- gendiv(my_bedgraph)
gendiv <- function(bdg) {
  gavg <- hwglabr2::average_chr_signal(bdg)$genome_avrg
  bdg$score <- bdg$score / gavg
  return(bdg)
}
