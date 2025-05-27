#' Annotate Genomic Ranges
#'
#' Labels windows as centromere/telomere/other based on proximity thresholds.
#'
#' @param df Data frame with genomic coordinates.
#' @param cen_df Data frame with centromere coordinates.
#' @param gen_chr_len Data frame with chromosome lengths.
#' @param binsize Centromere proximity in bp.
#' @param Telo_binsize Telomere proximity in bp.
#' @return Annotated data frame.
annotate_ranges <- function(df, cen_df, gen_chr_len, binsize = 10000, Telo_binsize = 20000) {
  df$c <- NA
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(cen_df)) {
      if (df$seqnames[i] == cen_df$seqnames[j] &&
          df$start[i] >= cen_df$start[j] - binsize &&
          df$start[i] <= cen_df$start[j] + binsize) {
        df$c[i] <- "centromere"
        break
      }
    }
    
    if (is.na(df$c[i])) {
      for (j in 1:nrow(gen_chr_len)) {
        if (df$seqnames[i] == gen_chr_len$chr[j] &&
            (df$start[i] <= Telo_binsize || df$start[i] >= (gen_chr_len$len[j] - Telo_binsize))) {
          df$c[i] <- "telomere"
          break
        }
      }
    }
    
    if (is.na(df$c[i])) {
      df$c[i] <- "non-centromere and non-telomere"
    }
  }
  df
}
