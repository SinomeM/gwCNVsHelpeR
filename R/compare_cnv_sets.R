
#' Compare two CNV sets
#'
#' Return all pairs of overlapping segments for a given sample with the same (exact
#' ) CN from two separate sets, e.g. two genotyping run in the same collection.
#'
#' @param cnvsA CNVs table
#' @param cnvsB CNVs table
#' @param gt_cn choose to check the exct CN or just the GT
#'
#' @export
#'
#' @import data.table

# compare two CNV sets, return all pairs of overlapping segments for each sample
# whit the same (exact) CN
compare_cnv_sets <- function(cnvsA, cnvsB, gt_cn = 'CN') {
  if (gt_cn == 'CN') {
    a <- copy(cnvsA[, .(sample_ID, chr, start, end, CN)])
    b <- copy(cnvsB[, .(sample_ID, chr, start, end, CN)])
  }
  if (gt_cn == 'GT') {
    a <- copy(cnvsA[, .(sample_ID, chr, start, end, GT)])
    b <- copy(cnvsB[, .(sample_ID, chr, start, end, GT)])
  }
  
  dt <- data.table()
  for (s in unique(c(cnvsA$sample_ID, cnvsB$sample_ID))) {
    aa <- a[sample_ID == s, ]
    bb <- b[sample_ID == s, ]
    if (aa[, .N] == 0 | bb[, .N] == 0) next
    setkey(bb, start, end)
    if (gt_cn == 'CN') {
      tmp <- foverlaps(aa, bb)[CN == i.CN & chr == i.chr, ]
      tmp[, ':=' (i.sample_ID = NULL, i.CN = NULL, i.chr = NULL)]
    }
    if (gt_cn == 'GT') {
      tmp <- foverlaps(aa, bb)[GT == i.GT & chr == i.chr, ]
      tmp[, ':=' (i.sample_ID = NULL, i.GT = NULL, i.chr = NULL)]
    }
    dt <- rbind(dt, tmp)
  }
  if (gt_cn == 'CN')
    colnames(dt) <- c('sample_ID', 'chr', 'startB', 'endB', 'CN', 'startA', 'endA')
  if (gt_cn == 'GT')
    colnames(dt) <- c('sample_ID', 'chr', 'startB', 'endB', 'GT', 'startA', 'endA')
  dt[, iou := (pmin(endA, endB) - pmax(startA, startB)) /
       (pmax(endA, endB) - pmin(startA, startB))]
  return(dt)
}
