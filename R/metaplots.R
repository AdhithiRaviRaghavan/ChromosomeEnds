#' Get GRanges for X elements with adjacent Y′ elements
#' Matches by seqnames, strand, and valid Y′ types
get_XwithY_GRanges <- function(genome = "SK1Yue") {
  gff <- get_gff(genome)
  xelement_df <- data.frame(gff[grep("X_element", gff$type)])
  yelement_df <- data.frame(gff[grep("Y_prime_element|Y_prime_element_partial", gff$type)])

  matched_rows <- purrr::map_dfr(1:nrow(yelement_df), function(i) {
    y <- yelement_df[i, ]
    dplyr::filter(
      xelement_df,
      seqnames == y$seqnames,
      strand == y$strand,
      type == "X_element"
    )
  })

  matched_unique <- dplyr::distinct(matched_rows)
  gr <- GenomicRanges::makeGRangesFromDataFrame(matched_unique)
  flank <- floor(GenomicRanges::width(gr))
  GenomicRanges::start(gr) <- GenomicRanges::start(gr) - flank
  GenomicRanges::end(gr) <- GenomicRanges::end(gr) + flank
  GenomicRanges::trim(gr)
}

#' Compute metaplot signal around X elements with Y′ ends
#'
#' Normalizes signal to genome-wide average, rescales with spike-in factor,
#' and extracts meta-signal with 95% CI.
#'
#' @param signal GenomicRanges object with score column
#' @param normalizationfactor Numeric (default = 1)
#' @param regions GRanges for regions of interest
#' @return Data frame with mean and 95% CI
#' @export
MetaplotXwithY <- function(signal, normalizationfactor = 1, regions) {
  gendiv <- function(bdg) {
    gavg <- hwglabr2::average_chr_signal(bdg)$genome_avrg
    message("Genome-wide average signal: ", round(gavg, 3))
    bdg$score <- bdg$score / gavg
    return(bdg)
  }
  signal_divided <- gendiv(signal)
  GenomicRanges::score(signal_divided) <- GenomicRanges::score(signal_divided) * normalizationfactor
  n_windows <- floor(1000 / 10)
  message("Computing signal matrix over provided GRanges...")
  signal_matrix <- EnrichedHeatmap::normalizeToMatrix(
    signal_divided, regions,
    value_column = "score", mean_mode = "absolute",
    extend = 0, k = n_windows, empty_value = NA,
    smooth = FALSE, target_ratio = 1
  )
  hwglabr2::signal_mean_and_ci(signal_matrix, ci = 0.95, rep_bootstrap = 1000, na_rm = TRUE)
}

#' Get GRanges for X elements without adjacent Y′ elements
get_Xonly_GRanges <- function(genome = "SK1Yue") {
  gff <- get_gff(genome)
  xelement_df <- data.frame(gff[grep("X_element", gff$type)])
  yelement_df <- data.frame(gff[grep("Y_prime_element|Y_prime_element_partial", gff$type)])

  matched_rows <- purrr::map_dfr(1:nrow(yelement_df), function(i) {
    y <- yelement_df[i, ]
    dplyr::filter(
      xelement_df,
      seqnames == y$seqnames,
      strand == y$strand,
      type == "X_element"
    )
  })

  unmatched <- dplyr::anti_join(xelement_df, matched_rows, by = "ID")
  unmatched_gr <- GenomicRanges::makeGRangesFromDataFrame(unmatched)
  flank <- floor(GenomicRanges::width(unmatched_gr))
  GenomicRanges::start(unmatched_gr) <- GenomicRanges::start(unmatched_gr) - flank
  GenomicRanges::end(unmatched_gr) <- GenomicRanges::end(unmatched_gr) + flank
  GenomicRanges::trim(unmatched_gr)
}

#' Compute metaplot signal around X-only elements
#'
#' Normalizes signal to genome-wide average, rescales with spike-in factor,
#' and extracts meta-signal with 95% CI.
#'
#' @param signal GenomicRanges object with score column
#' @param normalizationfactor Numeric (default = 1)
#' @param regions GRanges for regions of interest
#' @return Data frame with mean and 95% CI
#' @export
MetaplotXonly <- function(signal, normalizationfactor = 1, regions) {
  gendiv <- function(bdg) {
    gavg <- hwglabr2::average_chr_signal(bdg)$genome_avrg
    message("Genome-wide average signal: ", round(gavg, 3))
    bdg$score <- bdg$score / gavg
    return(bdg)
  }
  signal_divided <- gendiv(signal)
  GenomicRanges::score(signal_divided) <- GenomicRanges::score(signal_divided) * normalizationfactor
  n_windows <- floor(1000 / 10)
  message("Computing signal matrix over X-only regions...")
  signal_matrix <- EnrichedHeatmap::normalizeToMatrix(
    signal_divided, regions,
    value_column = "score", mean_mode = "absolute",
    extend = 0, k = n_windows, empty_value = NA,
    smooth = FALSE, target_ratio = 1
  )
  hwglabr2::signal_mean_and_ci(signal_matrix, ci = 0.95, rep_bootstrap = 1000, na_rm = TRUE)
}

#' Get GRanges for Y′ elements (with flanking ±0.5x width)
get_Yprime_GRanges <- function(genome = "SK1Yue") {
  gff <- get_gff(genome)
  yelements <- gff[grep("Y_prime_element", gff$type)]
  flank <- floor(GenomicRanges::width(yelements) / 2)
  GenomicRanges::start(yelements) <- GenomicRanges::start(yelements) - flank
  GenomicRanges::end(yelements) <- GenomicRanges::end(yelements) + flank
  GenomicRanges::trim(yelements)
}

#' Compute metaplot signal across Y′ elements
#'
#' Normalizes signal to genome-wide average and computes metaplot with 95% CI
#'
#' @param signal GenomicRanges object with score
#' @param regions GRanges of Y′ elements (e.g., from get_Yprime_GRanges())
#' @return Data frame with meta-signal mean and CI
#' @export
MetaplotY <- function(signal, regions) {
  gendiv <- function(bdg) {
    gavg <- hwglabr2::average_chr_signal(bdg)$genome_avrg
    bdg$score <- bdg$score / gavg
    return(bdg)
  }

  signal_divided <- gendiv(signal)
  n_windows <- floor(1000 / 10)

  message("Calculating signal in Y′ elements...")

  signal_matrix <- EnrichedHeatmap::normalizeToMatrix(
    signal_divided, regions,
    value_column = "score", mean_mode = "absolute",
    extend = 0, k = n_windows, empty_value = NA,
    smooth = FALSE, target_ratio = 1
  )

  hwglabr2::signal_mean_and_ci(signal_matrix, ci = 0.95, rep_bootstrap = 1000, na_rm = TRUE)
}
