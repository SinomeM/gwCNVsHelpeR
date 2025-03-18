

#                   #
# Combine CNVR sets #
#                   #

combine_CNVRs_across <- function(cnvrs1, cnvrs2, chr_arms, min_iou = 0.9) {
  a <- copy(cnvrs1[[2]])
  b <- copy(cnvrs2[[2]])
  a[, ix := 1:.N]
  b[, ix := 1:.N]

  acnvs <- copy(cnvrs1[[1]])
  bcnvs <- copy(cnvrs2[[1]])

  for (i in 1:chr_arms[, .N]) {
    # one chr arm at the time reduce computation
    cc <- chr_arms[i]
    aa <- a[chr == cc$chr & start <= cc$end & end >= cc$start, ]
    bb <- b[chr == cc$chr & start <= cc$end & end >= cc$start, ]
    if (nrow(aa) == 0 | nrow(bb) == 0) next
    # compare with the other, there is no avoiding using one as the "reference"
    # but it should not have a big impact since we do only one round and
    # select the best match
    for (ii in 1:aa[, .N]) {
      aaa <- aa[ii]
      bbb <- bb[start <= aaa$end & end >= aaa$start, ]
      bbb[, iou := (pmin(end, aaa$end) - pmax(start, aaa$start) + 1) /
          (pmax(end, aaa$end) - pmin(start, aaa$start) + 1)]
      bbb <- bbb[iou >= min_iou, ]
      # let's you know if multiple matches have been found
      if (bbb[, .N] > 1) message(bbb[, .N])
      # this should be safe even if there is more than one match
      if (bbb[, .N] != 0)
        a[ix == aaa$ix, match_ix := bbb[iou == max(iou), ix]]
    }
  }
  
  # combine the CNVRs
  a[, new_name := paste0('a_', ix)]
  b <- merge(b, a[, .(match_ix, new_name)], by.x = 'ix', by.y = 'match_ix', all.x = T)
  b[is.na(new_name), new_name := paste0('b_', ix)]
  
  # update the CNVs
  acnvs <- merge(acnvs, a[, .(CNVR, new_name)], by = 'CNVR')
  bcnvs <- merge(bcnvs, b[, .(CNVR, new_name)], by = 'CNVR')
  cnvs <- rbind (acnvs, bcnvs, fill = T)
  cnvs[, CNVR := new_name]
  
  # recreate the CVNRs
  cnvrs <- cnvs[, .N, by = CNVR]
  setnames(cnvrs, 'N', 'totN')
  cnvrs <- merge(cnvrs, cnvs[GT == 1, .N, by = CNVR], all.x = T)
  cnvrs[is.na(N), N := 0]
  setnames(cnvrs, 'N', 'N_dels')
  cnvrs <- merge(cnvrs, cnvs[GT == 2, .N, by = CNVR], all.x = T)
  cnvrs[is.na(N), N := 0]
  setnames(cnvrs, 'N', 'N_dups')
  cnvrs <- merge(cnvrs, cnvs[, median(chr), by = CNVR], all.x = T)
  setnames(cnvrs, 'V1', 'chr')
  cnvrs <- merge(cnvrs, cnvs[, median(start), by = CNVR], all.x = T)
  setnames(cnvrs, 'V1', 'start')
  cnvrs <- merge(cnvrs, cnvs[, median(end), by = CNVR], all.x = T)
  setnames(cnvrs, 'V1', 'end')
  cnvrs <- merge(cnvrs, cnvs[, median(numsnp), by = CNVR], all.x = T)
  setnames(cnvrs, 'V1', 'numsnp')
  cnvrs <- merge(cnvrs, cnvs[, mean(p_false), by = CNVR], all.x = T)
  setnames(cnvrs, 'V1', 'mean_Pfalse')
  cnvrs <- merge(cnvrs, cnvs[, median(p_false), by = CNVR], all.x = T)
  setnames(cnvrs, 'V1', 'median_Pfalse')

  setnames(cnvrs, 'totN', 'N')

  return(list(cnvs, cnvrs))
}



# add separate count for dels and dups and mean p_false
recreate_CNVRs <- function(cnvrs, cnvs) {
  a <- copy(cnvrs[[1]])
  b <- copy(cnvrs[[2]])
  setnames(b, 'n', 'tot_N')
  a <- merge(a[, .(sample_ID, chr, start, end, numsnp, CNVR)],
             cnvs[, .(sample_ID, chr, start, end, GT, pred, pred_prob, p_false)])

  b <- merge(b, a[GT == 1, .N, by = CNVR], all.x = T)
  setnames(b, 'N', 'N_dels')
  b <- merge(b, a[GT == 2, .N, by = CNVR], all.x = T)
  setnames(b, 'N', 'N_dups')
  b <- merge(b, a[, mean(p_false), by = CNVR], all.x = T)
  setnames(b, 'V1', 'mean_Pfalse')
  b <- merge(b, a[, median(p_false), by = CNVR], all.x = T)
  setnames(b, 'V1', 'median_Pfalse')

  setnames(b, 'tot_N', 'N')

  return(b)
}

