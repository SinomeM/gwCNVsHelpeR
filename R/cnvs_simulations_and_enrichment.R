


#                   #
# CNVRs simulations #
#                   #

inputs_cnvrs_simulation <- function(cnvs, cnvrs, snps, chr_arms,
                                    min_snps_chr_arm = NULL) {
  # restrict to 1:22 always
  chr_arms <- chr_arms[chr %in% 1:22, ]

  # order CNVs and SNPs by start, per chromosome
  setorder(cnvs, chr, start)
  setorder(snps, Chr, Position)
  setorder(cnvs, chr, start)
  if ('real_numsnp' %in% colnames(cnvs))
      cnvs[, numsnp := real_numsnp]

  # add an index per chromosome arm to the SNPs, and the arm of each SNP and CNV
  for (i in 1:chr_arms[, .N]) {
    dt <- chr_arms[i]
    snps[Chr == dt$chr & between(Position, dt$start, dt$end),
           ':=' (a_ix = 1:.N, arm = dt$arm_ID)]
    cnvs[chr == dt$chr & between(start, dt$start, dt$end), arm := dt$arm_ID]
    cnvrs[chr == dt$chr & between(start, dt$start, dt$end), arm := dt$arm_ID]
  }

  # add an index to the CNV table and the full SNP table
  cnvs[, ix := 1:.N]
  snps[, ix := 1:.N]

  cnvrs[, ix := 1:.N]

  if (!is.null(min_snps_chr_arm))
    chr_arms <- chr_arms[arm_ID %in% snps[, .N, by = arm][
                                          N >= min_snps_chr_arm, arm], ]

  # add max SNP ix for each arm
  for (i in 1:chr_arms[, .N])
    chr_arms[i, ':=' (max_ix = snps[arm == chr_arms[i, arm_ID], max(ix)],
                      min_ix = snps[arm == chr_arms[i, arm_ID], min(ix)])]

  return(list(cnvs[, .(sample_ID, chr, arm, numsnp, GT, CN, CNVR)],
              cnvrs[, .(CNVR, chr, arm, start, end, n, numsnp)],
              snps[, !('Name')], chr_arms))
}

cnvrs_simulation <- function(cnvs, cnvrs, snps, chr_arms, length_OK_prop = 0.20,
                             verbose = F) {

  cnvrs[, length := end- start + 1]

  cnvrs <- new_start_end_SNPs(cnvrs, snps, chr_arms, length_OK_prop)

  i <- 1
  cnvrs[new_length_OK == T, length_OK_iterations := i-1]
  # reprocess those with length very different from before until they fit,
  # slow but fair (statistically speaking). 1% of CNVRs are allowed to have
  # 'bad' new_length to speed things up
  while (cnvrs[new_length_OK == F, .N] > cnvrs[, .N] * 0.01) {
    if (verbose) message(i, ': ', cnvrs[new_length_OK == F, .N], '\n')
    i <- i+1
    # if stuck increase by 25% every 50 iterations
    if (i %% 50 == 0) length_OK_prop <- length_OK_prop + length_OK_prop * 0.25
    cnvrs <- rbind(cnvrs[new_length_OK == T, ],
                   new_start_end_SNPs(cnvrs[new_length_OK == F, ],
                                      snps, chr_arms, length_OK_prop))
    cnvrs[new_length_OK == T & is.na(length_OK_iterations),
            length_OK_iterations := i-1]
  }

  cnvrs[is.na(length_OK_iterations), length_OK_iterations := i]
  # Now simulate CNVs based on the new CNVR position
  # They are allowed to move a little bit
  a <- cnvs[sample(1:.N, round(.N/2)), ]
  b <- fsetdiff(cnvs, a)
  # half will inherit the start, half the end +/- 0 to 10% numsnp wiggling
  cnvrs[, snps_wiggle :=
        round(sample(seq(0, 0.1, length.out = 10), .N, replace = T) * numsnp)]
  a <- merge(a[, !('arm')],cnvrs[, .(CNVR, arm, new_start_SNP,
                                     snps_wiggle, min_ix, max_ix)], by = 'CNVR')
  a[, new_start_SNP := new_start_SNP + (snps_wiggle * sample(c(1, -1), .N,
                                                             replace = T))]
  a[, new_end_SNP := new_start_SNP + numsnp]
  b <- merge(b[, !('arm')], cnvrs[, .(CNVR, arm, new_end_SNP,
                                      snps_wiggle, min_ix, max_ix)], by = 'CNVR')
  b[, new_end_SNP := new_end_SNP + (snps_wiggle * sample(c(1, -1), .N, replace = T))]
  b[, new_start_SNP := new_end_SNP - numsnp]

  # add step to check if any CNV fell out of the chromosome when CNVR is
  # close to an end
  cnvs <- check_boundaries(rbind(a, b))

  cnvs <- snp_to_coord(cnvs, snps)

  # extract chr from the arm
  cnvs[, new_chr := as.integer(gsub('\\D', '', arm))]
  cnvrs[, new_chr := as.integer(gsub('\\D', '', arm))]

  # return the new CNV table
  cnvs <- cnvs[, .(CNVR, new_chr, arm, new_start, new_end, sample_ID, GT, CN, numsnp)]
  colnames(cnvs) <- c('CNVR_origin', 'chr', 'arm', 'start', 'end', 'sample_ID',
                      'GT', 'CN', 'numsnp')
  cnvrs <- cnvrs[, .(CNVR, new_chr, arm, new_start, new_end, numsnp,
                     n, length_OK_iterations)]
  colnames(cnvrs) <- c('CNVR_origin', 'chr', 'arm', 'start', 'end', 'numsnp',
                       'N', 'length_OK_iterations')

  return(list(cnvs, cnvrs))

}


new_start_end_SNPs <- function(dt, snps, chr_arms, length_OK_prop) {
  # if in the while there is already min_ix and max_ix
  if ('min_ix' %in% colnames(dt))
    dt <- dt[, !c('min_ix', 'max_ix', 'new_start', 'new_end')]
  # sample new end from the entire SNPs table, assign chr amr based on that
  dt[, new_end_SNP := sample(snps$ix, .N, replace = T)]
  dt <- merge(dt[, !('arm')], snps[, .(ix, arm)],
                by.x = 'new_end_SNP', by.y = 'ix')
  # add max SNP ix of the arm to each new CNV
  dt <- merge(dt, chr_arms[, .(arm_ID, max_ix, min_ix)],
                by.x = 'arm', by.y = 'arm_ID')

  dt <- check_boundaries(dt)

  dt <- snp_to_coord(dt, snps)

  # new length must be comparable
  dt[, ':=' (new_length = new_end - new_start + 1, new_length_OK = F,
                min_length = length - (length_OK_prop * length),
                max_length = length + (length_OK_prop * length))]
  dt[between(new_length, min_length, max_length), new_length_OK := T]

  return(dt)
}


check_boundaries <- function(dt) {
  # check if the end is outside the arm, if it is move it back
  dt[, over := new_end_SNP - max_ix]
  dt[over > 0, new_end_SNP := new_end_SNP - over]
  dt[, new_start_SNP := new_end_SNP - numsnp]
  # check if the start is outside the arm, if it is move it foreword
  dt[, under := min_ix - new_start_SNP]
  dt[under > 0, new_start_SNP := new_start_SNP + under]
  dt[, new_end_SNP := new_start_SNP + numsnp]
}

snp_to_coord <- function(dt, snps) {
  # now start and end SNP are OK, move to coordinates
  dt <- merge(dt, snps[, .(ix, Position)],
                by.x = 'new_start_SNP', by.y = 'ix')
  dt <- merge(dt, snps[, .(ix, Position)],
                by.x = 'new_end_SNP', by.y = 'ix')
  setnames(dt, c('Position.x', 'Position.y'), c('new_start', 'new_end'))
  return(dt)
}


#                 #
# CNVs enrichment #
#                 #


count_annotation_full_cnvrs <- function(x, base_path, bins) {
  cnv_files <- list.files(base_path)
  cnvs <- readRDS(paste0(base_path, '/', cnv_files[x]))[[1]]
  tmp <- binned_cnvs(cnvs, format = 'count', bins = bins)
  
  # add back bins with N = 0 if necessary
  dt <- tmp[[2]][!is.na(ix), .(ix, GT, N)]
  tmp <- bins[, .(ix)]
  a <- copy(tmp)
  b <- copy(tmp)
  tmp <- rbind(a[, GT := 1], b[, GT := 2])
  dt <- merge(dt, tmp, all.y = T)
  dt[is.na(N), N := 0]

  return(dt)
}


save_baseline <- function(pt, bins, save_pt, bins_name) {
  len <- length(list.files(pt))

  baseline <- BiocParallel::bplapply(seq_len(len), function(x)
                                     count_annotation_full_cnvrs(x, pt, bins))

  dt <- unique(do.call(rbind, baseline)[, lapply(.SD, median), by = c('ix', 'GT')])
  colnames(dt) <- c(bins_name, 'GT', 'median_N')
  fwrite(dt, save_pt, sep = '\t')
}


get_enrichment_table <- function(cnvs, baseline_pt, bins, bins_name) {

  bins <- copy(bins)
  dt <- binned_cnvs(cnvs, format = 'count', bins = bins)[[2]][
          !is.na(ix), .(ix, GT, N)]
  setnames(dt, 'ix', bins_name)
  baseline <- fread(baseline_pt)
  dt <- merge(baseline, dt, by = c(bins_name, 'GT'), all.x = T)
  dt[is.na(N), N := 0]

  # combine dels and dups
  dt <- merge(dt[GT == 1, !c('GT', 'N_simul')],
              dt[GT == 2, !c('GT', 'N_simul')], by = bins_name, all = T)
  colnames(dt) <- c(bins_name, 'median_dels', 'N_dels',
                    'median_dups', 'N_dups')
  dt[is.na(median_dels), median_dels := 0]
  dt[is.na(median_dups), median_dups := 0]
  dt[is.na(N_dels), N_dels := 0]
  dt[is.na(N_dups), N_dups := 0]
  dt[, ':=' (median_all = median_dels + median_dups, N_all = N_dels + N_dups)]

  # fisher test pval
  n_cnvs <- cnvs[, .N]
  n_dels <- cnvs[GT == 1, .N]
  n_dups <- cnvs[GT == 2, .N]
  for (i in 1:dt[, .N]) {
  #for (i in 1:100) {
    a <- dt[i, matrix(c(N_all, n_cnvs - N_all, round(median_all),
                        n_cnvs - round(median_all)), ncol = 2, byrow = T)]
    b <- dt[i, matrix(c(N_dels, n_dels - N_dels, round(median_dels),
                        n_dels - round(median_dels)), ncol = 2, byrow = T)]
    c <- dt[i, matrix(c(N_dups, n_dups - N_dups, round(median_dups),
                        n_dups - round(median_dups)), ncol = 2, byrow = T)]
    dt[i, ':=' (p_all = fisher.test(a)$p.value, p_dels = fisher.test(b)$p.value,
                p_dups = fisher.test(c)$p.value)]
  }

  dt[, ':=' (log2FC_all = log((N_all+1) / (median_all+1), 2),
             log2FC_dels = log((N_dels+1) / (median_dels+1), 2),
             log2FC_dups = log((N_dups+1) / (median_dups+1), 2),
             padj_all = p.adjust(p_all, 'fdr'),
             padj_dels = p.adjust(p_dels, 'fdr'),
             padj_dups = p.adjust(p_dups, 'fdr'))]
  dt[, ':=' (log10_padj_all = log(padj_all, 10), log10_padj_dels = log(padj_dels, 10),
             log10_padj_dups = log(padj_dups, 10))]

  return(dt)
}



# Version for genes (also considers exonic or not)

count_annotation_genes_cnvrs_bplapply <- function(x, base_path, genes, exons) {
  cnv_files <- list.files(base_path)
  cnvs <- readRDS(paste0(base_path, '/', cnv_files[x]))[[1]]
  
  return(count_annotation_genes_cnvrs(cnvs, genes, exons))
}

count_annotation_genes_cnvrs <- function(cnvs, genes, exons) {
  genes[, ix := ens_ID]
  exons[, ':=' (ix = exon_ID, start = exon_start, end = exon_end)]
  # counts for genes and exons
  g_tmp <- binned_cnvs(cnvs, format = 'count', bins = genes)
  e_tmp <- binned_cnvs(cnvs, format = 'count', bins = exons)
  # collapse exons on a signle ens_ID
  e_tmp <- e_tmp[[1]][!is.na(exon_ID) & !is.na(sample_ID), .(ens_ID, sample_ID, GT)]
  e_tmp <- unique(e_tmp)

  # complete or non complete (on the genes)
  a <- g_tmp[[1]]
  a[i.start <= start & i.end >= end, complete := T]
  a[is.na(complete), complete := F]
  # exonic or not
  b <- unique(e_tmp)
  b[, exonic := T]
  b[is.na(exonic), exonic := F]
  
  # combine the two
  dt <- merge(a[, .(sample_ID, ens_ID, GT, complete)], b,
              by = c('sample_ID', 'ens_ID', 'GT'), all = T)
  dt <- dt[!is.na(sample_ID) & !is.na(ens_ID), ]
  dt[is.na(exonic), exonic := F]
  
  dt <- dt[, .N, by = c('ens_ID', 'GT', 'complete', 'exonic')]
  
  # add back lines with N = 0 if necessary
  tmp <- genes[, .(ens_ID)]
  t1 <- copy(tmp); t1[, ':=' (GT = 1, exonic = T, complete = T)]
  t2 <- copy(tmp); t2[, ':=' (GT = 2, exonic = T, complete = T)]
  t3 <- copy(tmp); t3[, ':=' (GT = 1, exonic = F, complete = T)]
  t4 <- copy(tmp); t4[, ':=' (GT = 2, exonic = F, complete = T)]
  t5 <- copy(tmp); t5[, ':=' (GT = 1, exonic = T, complete = F)]
  t6 <- copy(tmp); t6[, ':=' (GT = 2, exonic = T, complete = F)]
  t7 <- copy(tmp); t7[, ':=' (GT = 1, exonic = F, complete = F)]
  t8 <- copy(tmp); t8[, ':=' (GT = 2, exonic = F, complete = F)]
  tmp <- rbind(t1, t2, t3, t4, t5, t6, t7, t8)
  dt <- merge(dt, tmp, all.y = T)
  dt[is.na(N), N := 0]
  
  return(dt)
}


save_baseline_genes <- function(pt, genes, exons, save_pt, bins_name) {
  len <- length(list.files(pt))
  
  baseline <- BiocParallel::bplapply(seq_len(len), function(x)
                                     count_annotation_genes_cnvrs_bplapply(x, pt, genes,
                                                                           exons))

  saveRDS(baseline, paste0(save_pt, '_tmp_baseline.rds'))
  
  dt <- unique(do.call(rbind, baseline)[, lapply(.SD, median), 
                                          by = c('ens_ID', 'GT', 'complete', 'exonic')])
  colnames(dt) <- c(bins_name, 'GT', 'complete', 'exonic', 'median_N')
  dt[, N_simul := len]
  fwrite(dt, save_pt, sep = '\t')
  # if succesfull delete the tmp file 
  unlink(paste0(save_pt, '_tmp_baseline.rds'))
}


get_enrichment_table_genes <- function(cnvs, baseline_pt, genes, exons) {

  dt <- count_annotation_genes_cnvrs(cnvs, genes, exons)
  baseline <- fread(baseline_pt)  
  dt <- merge(baseline, dt, by = c('ens_ID', 'GT', 'complete', 'exonic'), all.x = T)
  dt[is.na(N), N := 0]

  dele <- dt[GT == 1 & exonic == T, .(ens_ID, median_N, N)]
  dupe <- dt[GT == 2 & exonic == T, .(ens_ID, median_N, N)]
  delc <- dt[GT == 1 & exonic == T & complete == T, .(ens_ID, median_N, N)]
  dupc <- dt[GT == 2 & exonic == T & complete == T, .(ens_ID, median_N, N)]
  delp <- dt[GT == 1 & exonic == T & complete == F, .(ens_ID, median_N, N)]
  dupp <- dt[GT == 2 & exonic == T & complete == F, .(ens_ID, median_N, N)]
  
  dele <- merge(dele[,sum(median_N), by = 'ens_ID'], 
                dele[, sum(N), by = 'ens_ID'], by = 'ens_ID')
  dupe <- merge(dupe[,sum(median_N), by = 'ens_ID'],
                dupe[, sum(N), by = 'ens_ID'], by = 'ens_ID')
  delc <- merge(delc[,sum(median_N), by = 'ens_ID'],
                delc[, sum(N), by = 'ens_ID'], by = 'ens_ID')
  dupc <- merge(dupc[,sum(median_N), by = 'ens_ID'],
                dupc[, sum(N), by = 'ens_ID'], by = 'ens_ID')
  delp <- merge(delp[,sum(median_N), by = 'ens_ID'],
                delp[, sum(N), by = 'ens_ID'], by = 'ens_ID')
  dupp <- merge(dupp[,sum(median_N), by = 'ens_ID'],
                dupp[, sum(N), by = 'ens_ID'], by = 'ens_ID')

  setkey(dele, 'ens_ID'); setkey(dupe, 'ens_ID')
  setkey(delc, 'ens_ID'); setkey(dupc, 'ens_ID')
  setkey(delp, 'ens_ID'); setkey(dupp, 'ens_ID')
  
  dt <- dele[dupe,][delc,][dupc,][delp,][dupp,]

  # combine dels and dups
  colnames(dt) <- c('ens_ID', 'median_dels', 'N_dels',
                    'median_dups', 'N_dups', 
                    'median_dels_complete', 'N_dels_complete',
                    'median_dups_complete', 'N_dups_complete', 
                    'median_dels_partial', 'N_dels_partial',
                    'median_dups_partial', 'N_dups_partial')
  dt[is.na(median_dels), median_dels := 0]
  dt[is.na(median_dups), median_dups := 0]
  dt[is.na(N_dels), N_dels := 0]
  dt[is.na(N_dups), N_dups := 0]
  dt[, ':=' (median_all = median_dels + median_dups, N_all = N_dels + N_dups)]
  
  dt[is.na(median_dels_complete), median_dels_complete := 0]
  dt[is.na(median_dups_complete), median_dups_complete := 0]
  dt[is.na(N_dels_complete), N_dels_complete := 0]
  dt[is.na(N_dups_complete), N_dups_complete := 0]
  dt[, ':=' (median_all_complete = median_dels_complete + median_dups_complete,
             N_all_complete = N_dels_complete + N_dups_complete)]
  
  dt[is.na(median_dels_partial), median_dels_partial := 0]
  dt[is.na(median_dups_partial), median_dups_partial := 0]
  dt[is.na(N_dels_partial), N_dels_partial := 0]
  dt[is.na(N_dups_partial), N_dups_partial := 0]
  dt[, ':=' (median_all_partial = median_dels_partial + median_dups_partial,
             N_all_partial = N_dels_partial + N_dups_partial)]

  # fisher test pval
  n_cnvs <- cnvs[, .N]
  n_dels <- cnvs[GT == 1, .N]
  n_dups <- cnvs[GT == 2, .N]
  for (i in 1:dt[, .N]) {
    a <- dt[i, matrix(c(N_all, n_cnvs - N_all, round(median_all),
                        n_cnvs - round(median_all)), ncol = 2, byrow = T)]
    b <- dt[i, matrix(c(N_dels, n_dels - N_dels, round(median_dels),
                        n_dels - round(median_dels)), ncol = 2, byrow = T)]
    c <- dt[i, matrix(c(N_dups, n_dups - N_dups, round(median_dups),
                        n_dups - round(median_dups)), ncol = 2, byrow = T)]
    dt[i, ':=' (p_all = fisher.test(a)$p.value, p_dels = fisher.test(b)$p.value,
                p_dups = fisher.test(c)$p.value)]
                       
    a <- dt[i, matrix(c(N_all_complete, n_cnvs - N_all_complete,
                        round(median_all_complete),
                        n_cnvs - round(median_all_complete)), ncol = 2, byrow = T)]
    b <- dt[i, matrix(c(N_dels_complete, n_dels - N_dels_complete,
                        round(median_dels_complete),
                        n_dels - round(median_dels_complete)), ncol = 2, byrow = T)]
    c <- dt[i, matrix(c(N_dups_complete, n_dups - N_dups_complete,
                        round(median_dups_complete),
                        n_dups - round(median_dups_complete)), ncol = 2, byrow = T)]
    dt[i, ':=' (p_all_complete = fisher.test(a)$p.value,
                p_dels_complete = fisher.test(b)$p.value,
                p_dups_complete = fisher.test(c)$p.value)]
                          
    a <- dt[i, matrix(c(N_all_partial, n_cnvs - N_all_partial,
                        round(median_all_partial),
                        n_cnvs - round(median_all_partial)), ncol = 2, byrow = T)]
    b <- dt[i, matrix(c(N_dels_partial, n_dels - N_dels_partial,
                        round(median_dels_partial),
                        n_dels - round(median_dels_partial)), ncol = 2, byrow = T)]
    c <- dt[i, matrix(c(N_dups_partial, n_dups - N_dups_partial,
                        round(median_dups_partial),
                        n_dups - round(median_dups_partial)), ncol = 2, byrow = T)]
    dt[i, ':=' (p_all_partial = fisher.test(a)$p.value,
                p_dels_partial = fisher.test(b)$p.value,
                p_dups_partial = fisher.test(c)$p.value)]
  }


  dt[, ':=' (log2FC_all = log((N_all+1) / (median_all+1), 2),
             log2FC_dels = log((N_dels+1) / (median_dels+1), 2),
             log2FC_dups = log((N_dups+1) / (median_dups+1), 2),
             
             log2FC_all_complete = log((N_all_complete+1)/(median_all_complete+1), 2),
             log2FC_dels_complete = log((N_dels_complete+1)/(median_dels_complete+1), 2),
             log2FC_dups_complete = log((N_dups_complete+1)/(median_dups_complete+1), 2),
             
             log2FC_all_partial = log((N_all_partial+1)/(median_all_partial+1), 2),
             log2FC_dels_partial = log((N_dels_partial+1)/(median_dels_partial+1), 2),
             log2FC_dups_partial = log((N_dups_partial+1)/(median_dups_partial+1), 2))]

  return(dt)
}


#                 #
# Genes Shuffling #
#                 #


shuffle_genes <- function(x, genes, chr_arms, snps) {

  set.seed(x)

  genes[, length := end - start + 1]
  chr_arms <- copy(chr_arms)
  chr_arms[, ix := 1:.N]
  #indexes for fisrt and last SNP on each arm
  for (i in chr_arms[, ix]) {
    cc <- chr_arms[ix == i, chr]
    st <- chr_arms[ix == i, start]
    en <- chr_arms[ix == i, end]
    s_tmp <- snps[Chr == cc & between(Position, st, en), ]
    if (s_tmp[, .N] == 0) next
    chr_arms[ix == i, ':=' (first_snp = s_tmp[, min(Position)],
                            last_snp = s_tmp[, max(Position)])]
  }
  # length covered by SNPs
  chr_arms[, covered_len := last_snp - first_snp + 1]
  chr_arms <- chr_arms[!is.na(covered_len), ]
  # index start end
  for (i in chr_arms[, ix]) {
    if (i == 1) {
      chr_arms[ix == i, ':=' (st_ix = as.numeric(0), 
                              en_ix = as.numeric(covered_len),
                              cum_len = covered_len)]
      prev <- 1
    }
    else {
      chr_arms[ix == i, st_ix := as.numeric(chr_arms[ix == prev, en_ix+1])]
      chr_arms[ix == i, en_ix := as.numeric(st_ix + covered_len)]
      chr_arms[ix == i, cum_len := as.numeric(chr_arms[ix == prev, cum_len]) +
                 covered_len]
      prev <- i
    }
  }

  dt <- copy(genes)
  max_ix <- chr_arms[, max(cum_len)]
  dt[, new_st_ix := sample(0:max_ix, .N)]
  # for each chromosomal arm, get the chr and go back to the correct dimension
  res <- data.table()
  for (i in chr_arms[, ix]) {
    cc <- chr_arms[ix == i, chr]
    cl <- chr_arms[ix == i, cum_len]
    st_ix <- chr_arms[ix == i, st_ix]
    en_ix <- chr_arms[ix == i, en_ix]
    st <- chr_arms[ix == i, first_snp]
    en <- chr_arms[ix == i, last_snp]
    
    tmp <- dt[between(new_st_ix, st_ix, en_ix)]
    
    tmp[, ':=' (new_chr = cc, new_start = new_st_ix - st_ix + st)]
    tmp[, new_end := new_start + length]
    # in the unlikely case a gene is "on the border" pull it in
    tmp[new_end >= st & new_start < st,
        ':=' (new_start = new_start + length, new_end = new_end + length)]
    tmp[new_end >= en & new_start < en,
        ':=' (new_start = new_start - length, new_end = new_end - length)]
    res <- rbind(res, tmp)
  }
  return(res[, .(ens_ID, new_chr, new_start, new_end)])
}


# Baseline function, modified for shuffled genes
save_baseline_shuffled_genes <- function(cnvs, pt_genes, genes, exons, save_pt) {
  genes_l <- readRDS(pt_genes)
  len <- length(genes_l)

  baseline <-
    BiocParallel::bplapply(seq_len(len), function(x)
                           count_annotation_genes_shuff_bplapply(x, cnvs, genes_l[[x]],
                                                                 genes, exons))

  saveRDS(baseline, paste0(save_pt, '_tmp_baseline.rds'))

  dt <- unique(do.call(rbind, baseline)[, lapply(.SD, median), 
                                          by = c('ens_ID', 'GT', 'complete', 'exonic')])
  colnames(dt) <- c('ens_ID', 'GT', 'complete', 'exonic', 'median_N')

  fwrite(dt, save_pt, sep = '\t')
  # if succesfull delete the tmp file 
  unlink(paste0(save_pt, '_tmp_baseline.rds'))
}

count_annotation_genes_shuff_bplapply <- function(x, cnvs, genes_x, genes, exons) {

  genes <- copy(genes)
  exons <- copy(unique(exons))

  # retrieve original boundaries
  genes_x <- merge(genes_x, genes[, .(ens_ID, start, end)], by = 'ens_ID')
  # compute distance
  genes_x[, dist := new_start - start]
  # Update exons
  exons_x <- merge(unique(exons), genes_x[, .(ens_ID, new_chr, dist)], by = 'ens_ID')
  exons_x[, ':=' (chr = as.integer(new_chr), 
                  exon_start = exon_start + dist, exon_end = exon_end + dist)]

  # update colnames
  genes_x[, ':=' (chr = as.integer(new_chr), start = new_start, end = new_end)]

  return(count_annotation_genes_cnvrs(cnvs, genes_x, exons_x))
}




