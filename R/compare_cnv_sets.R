
#' Compare two CNV sets
#'
#' Return all pairs of overlapping segments for a given sample with the same (exact
#' ) CN from two separate sets, e.g. two genotyping run in the same collection.
#'
#' @param cnvsA CNVs table
#' @param cnvsB CNVs table
#'
#' @export
#'
#' @import data.table

# compare two CNV sets, return all pairs of overlapping segments for each sample
# whit the same (exact) CN
compare_cnv_sets <- function(cnvsA, cnvsB) {
  a <- copy(cnvsA[, .(sample_ID, chr, start, end, CN)])
  b <- copy(cnvsB[, .(sample_ID, chr, start, end, CN)])
  
  dt <- data.table()
  for (s in unique(c(cnvsA$sample_ID, cnvsB$sample_ID))) {
    aa <- a[sample_ID == s, ]
    bb <- b[sample_ID == s, ]
    if (aa[, .N] == 0 | bb[, .N] == 0) next
    setkey(bb, start, end)
    tmp <- foverlaps(aa, bb)[CN == i.CN]
    tmp[, ':=' (i.sample_ID = NULL, i.CN = NULL, i.chr = NULL)]
    dt <- rbind(dt, tmp)
  }
  colnames(dt) <- c('sample_ID', 'chr', 'startB', 'endB', 'CN', 'startA', 'endA')
  dt[, iou := (pmin(endA, endB) - pmax(startA, startB)) /
       (pmax(endA, endB) - pmin(startA, startB))]
  return(dt)
}
