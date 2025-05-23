
# Functions that do not belong in a specific category and/or do not
# do much

# This code requires data.table


# Prepare CNV table for plotting etc
prepare_cnv_table <- function(cnvs, chr_arms) {
  if ('GT' %in% colnames(cnvs)) {
    cnvs[GT == 1, CNV := 'Deletions']
    cnvs[GT == 2, CNV := 'Duplications']
  }
  if ('pred' %in% colnames(cnvs)) {
    cnvs[pred == 1, prediction := 'false']
    cnvs[pred == 2, prediction := 'true_del']
    cnvs[pred == 3, prediction := 'true_dup']
  }
  if ('real_numsnp' %in% colnames(cnvs))
    cnvs[, ':=' (numsnp = real_numsnp, real_numsnp = NULL, ix = 1:.N)]
  cnvs[, length := end - start + 1]
  cnvs[, centre := round(start + length / 2)]
  for (i in 1:chr_arms[, .N]) {
    ca <- chr_arms[i]
    cnvs[chr == ca$chr & start <= ca$end & end >= ca$start, chr_arm := ca$arm_ID]
  }
}



# count how many markers are captured by the marker set only by markers with at
# least 5 carriers
markers_counts <- function(markers, cnvs) {
  a <- copy(markers[[1]])
  b <- copy(markers[[2]])

  good_m <- b[N >= 5, ix]
  capt <- a[ix %in% good_m, length(unique(paste0(sample_ID, chr, i.start)))]
  tot <- cnvs[, .N]

  message('Markers set capture ', capt, ' CNVs out of ', tot,
          ' (', round(capt/tot, 2), ')')
}



prepare_marker <- function(dt, chr_arms, samples) {

  dt[GT == 1, CNV := 'Deletions'][GT == 2, CNV := 'Duplications']
  dt[, centre := round(start + (end - start + 1)/2)]
  dt[, frequency := N / samples[, .N]]

  for (i in 1:chr_arms[, .N]) {
    ca <- chr_arms[i]
    dt[chr == ca$chr & start <= ca$end & end >= ca$start, chr_arm := ca$arm_ID]
  }
  return(dt)
}





#                      #
# Misc helper (deCODE) #
#                      #

# Prepare everything for export
process_markers_table <- function(dt) {
  dt[between(N, 1, 4), N := NA]

  if ('N_dels' %in% colnames(dt))
    dt[between(N_dels, 1, 4) | between(N_dups, 1, 4), ':=' (N_dels = NA, N_dups = NA)]

  return(dt)
}

# process_bins_gc_table <- function(dt) {
#   dt[N_dels < 5, ':=' (N_dels = NA, N_dups = NA)]
#   dt[N_dups < 5, ':=' (N_dups = NA, N_dels = NA)]
#   dt[N_all < 5, ':=' (N_all = NA, median_all = NA)]
# }
process_bins_gc_table <- process_markers_table


process_bins_gc_table_genes <- function(dt) {
  dt[between(N_dels, 1, 4) | between(N_dups, 1, 4),
       ':=' (N_dels = NA, N_dups = NA, N_dels_complete = NA, N_dups_complete = NA,
             N_dels_partial = NA, N_dups_partial = NA)]
  dt[between(N_all, 1, 4), ':=' (N_all = NA, median_all = NA)]
}

