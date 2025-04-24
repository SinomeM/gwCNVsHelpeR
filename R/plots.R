
#' Genomic distribution of CNVs
#'
#' Plot the genomic distribution of deletion and duplication
#' genome-wide (chromosomes 1 to 22). The y axes is the rescaled
#' log frequency.
#'
#' @param cnvs CNVs table
#'
#' @param chr_st_en As usual
#'
#' @param n_samples For the frequency calculation
#'
#' @param bin_size Size of the unit
#'
#' @param return_pl_list return all chromosomes individually
#'   instead of four rows
#'
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @import ggpubr
#' @import patchwork

# This function creates the huge plot of the CNV frequencies in chrs 1:22 as
# in the Genome wide CNV draft in the thesis. Positive values for DELs and negative
# values for DUPs. The y axes is a rescaled log frequency, might need more parameters
cnvs_distr_big <- function(cnvs, chr_st_en, n_samples, bin_size = 250000, return_pl_list = F) {

  colours <- c(Deletions = '#E69F00', Duplications = '#56B4E9',
             True = '#009E73', False = '#D55E00', Unknown = '#F0E442')

  bins <- binned_cnvs(cnvs, format = 'count',
                      CNValidatron:::binned_genome(chr_st_en, bin_size))[[2]]
  bins[, centre := round(start + (end - start + 1)/2)]
  bins[GT == 1, CNV := 'Deletions'][GT == 2, CNV := 'Duplications']

  # max and min from non-empty bins
  merge(bins[N != 0, min(start), by = chr],
        bins[N != 0, max(end), by = chr], by = 'chr') -> chr_boundaries
  colnames(chr_boundaries) <- c('chr', 'start', 'end')
  chr_boundaries[, chr := as.character(chr)]
  chr_boundaries <- merge(chr_boundaries, chr_st_en[chr %in% 1:22, ], by = 'chr')
  chr_boundaries[, chr := as.integer(chr)]
  # add one bin
  chr_boundaries[start.x != start.y, start.x := start.x - 250000]
  chr_boundaries[end.x != end.y, end.x := end.x + 250000]
  chr_boundaries <- chr_boundaries[, .(chr, start.x, end.x, centromere)]
  colnames(chr_boundaries) <- c('chr', 'start', 'end', 'centromere')

  # four rows of ~ 700 units in width
  chromosome_groups <- list(1:3, 4:7, 8:12, 13:22)
  sides <- list(c(255, 245, 200),
                c(190, 180, 170, 160),
                c(147, 142, 137, 137, 137),
                c(98, 88, 78, 88, 83, 78, 58, 63, 33, 33))

  # y-axis scale and labels
  y_labs <- c(1/100, 1/1000, 1/10000, 1/100000)
  y_breaks <- log(y_labs) - log(1/100000)
  y_labs <- c(y_labs, y_labs[3:1])
  y_breaks <- c(y_breaks, -y_breaks[3:1])

  # Scaled log frequency
  bins[, freq := N / n_samples]
  bins[, scaled_freq := log(freq) + (-log(1/100000))]
  bins[N == 0, scaled_freq := 0]
  bins[GT == 2, scaled_freq := -scaled_freq]

  # create all plots
  pl_list <- list()
  for (i in 1:22) {
    tmp <- bins[chr == i, ]
    pl <- ggplot(tmp) +
      geom_bar(aes(x = centre/1000000, y = scaled_freq, fill = CNV),
               stat = 'identity', show.legend = F) +
      geom_vline(xintercept = chr_boundaries[chr == i, centromere/1000000],
                 linetype = "dashed", color = "black") + theme_bw() +
      scale_y_continuous(breaks = y_breaks, labels = y_labs * 1000, limits=c(-9,9)) +
      scale_fill_manual(values = colours) +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlim(chr_boundaries[chr == i, (start-1)/1000000], chr_boundaries[chr == i, (end + 1)/1000000 + 1])
    # remove y tick from all but the first
    if (!i %in% c(1, 4, 8, 13))
      pl <- pl + theme(axis.text.y = element_blank())

    pl_list[[i]] <- pl
  }

  if (return_pl_list) return(pl_list)

  # compose rows
  row1 <- pl_list[[1]] + pl_list[[2]] + pl_list[[3]] +
    plot_layout(widths = sides[[1]], nrow=1)
  row2 <- pl_list[[4]] + pl_list[[5]] + pl_list[[6]] + pl_list[[7]] +
    plot_layout(widths = sides[[2]], nrow=1)
  row3 <- pl_list[[8]] + pl_list[[9]] + pl_list[[10]] + pl_list[[11]] + pl_list[[12]] +
    plot_layout(widths = sides[[3]], nrow=1)
  row4 <- pl_list[[13]] + pl_list[[14]] + pl_list[[15]] + pl_list[[16]] + pl_list[[17]] +
    pl_list[[18]] + pl_list[[19]] + pl_list[[20]] + pl_list[[21]] + pl_list[[22]] +
    plot_layout(widths = sides[[4]], nrow=1, ncol = 10)

  # return list of the four rows
  return(list(row1, row2, row3, row4))
}

#' Genomic distribution of CNVs (T/F)
#'
#' Plot the genomic distribution of deletion and duplication
#' genome-wide (chromosomes 1 to 22). The y axes is the rescaled
#' log frequency.
#'
#' @param cnvs CNVs table
#'
#' @param chr_st_en As usual
#'
#'
#' @param bin_size Size of the unit
#'
#' @param gt 1 for deletions, 2 fro duplications
#'
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @import ggpubr
#' @import patchwork
#' @import CNValidatron

# Same as previous but need to be run twice (DELs and DUPs),
# this version is intended to show differences between validated
# CNVs (positive values) and false CNVs (negative values).
# Now the y axes scale in the absolute N value but I might move it to log
cnvs_distr_big_TF <- function(cnvs, chr_st_en, bin_size = 250000, gt) {

  a <- binned_cnvs(cnvs[GT == gt & Visual_Validation == 'True', ], format = 'count',
                   CNValidatron:::binned_genome(chr_st_en, bin_size))[[2]]
  b <- binned_cnvs(cnvs[GT == gt & Visual_Validation == 'False', ], format = 'count',
                   CNValidatron:::binned_genome(chr_st_en, bin_size))[[2]]
  a[, Visual_Validation := 'True']
  b[, ':=' (Visual_Validation = 'False', N = -N)]
  bins <- rbind(a, b)

  # log scale
  bins[N > 0, logN := log(N, 10)]
  bins[N < 0, logN := -log(-N, 10)]
  bins[N == 0, logN := 0]

  bins[, centre := round(start + (end - start + 1)/2)]
  bins[GT == 1, CNV := 'Deletions'][GT == 2, CNV := 'Duplications']

  # boundaries for each chr
  merge(bins[N != 0, min(start), by = chr],
        bins[N != 0, max(end), by = chr], by = 'chr') -> chr_boundaries
  colnames(chr_boundaries) <- c('chr', 'start', 'end')
  chr_boundaries[, chr := as.character(chr)]
  chr_boundaries <- merge(chr_boundaries, chr_st_en[chr %in% 1:22, ], by = 'chr')
  chr_boundaries[start.x > start.y, start.x := start.x - 250000][start.x < 0, start.x := 0]
  chr_boundaries[end.x < end.y, end.x := end.x + 250000][end.x > end.y, end.x := end.y]
  chr_boundaries <- chr_boundaries[, .(chr, start.x, end.x, centromere)]
  colnames(chr_boundaries) <- c('chr', 'start', 'end', 'centromere')

  # four rows of ~ 700 units in width
  chromosome_groups <- list(1:3, 4:7, 8:12, 13:22)
  sides <- list(c(255, 245, 200),
                c(190, 180, 170, 160),
                c(147, 142, 137, 137, 137),
                c(98, 88, 78, 88, 83, 78, 58, 63, 33, 33))
  ylims <- bins[, c(max(logN), min(logN))+0.1]

  pl_list <- list()
  for (i in 1:22) {
    pl <- ggplot(bins[chr == i, ]) +
      geom_bar(aes(x = centre/1000000, y = logN, fill = Visual_Validation),
               stat = 'identity', show.legend = F) +
      geom_vline(xintercept = chr_boundaries[chr == i, centromere/1000000],
                 linetype = "dashed", color = "black") + theme_bw() +
      scale_fill_manual(values = colours) +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlim(chr_boundaries[chr == i, (start-1)/1000000],
           chr_boundaries[chr == i, (end + 1)/1000000 + 1]) +
      ylim(ylims[2], ylims[1])

    if (!i %in% c(1, 4, 8, 13)) pl <- pl + theme(axis.text.y = element_blank())
    pl_list[[i]] <- pl
  }

  # compose and return the four rows
  return(list(
    pl_list[[1]] + pl_list[[2]] + pl_list[[3]] +
      plot_layout(widths = sides[[1]], nrow=1),
    pl_list[[4]] + pl_list[[5]] + pl_list[[6]] + pl_list[[7]] +
      plot_layout(widths = sides[[2]], nrow=1),
    pl_list[[8]] + pl_list[[9]] + pl_list[[10]] + pl_list[[11]] + pl_list[[12]] +
      plot_layout(widths = sides[[3]], nrow=1),
    pl_list[[13]] + pl_list[[14]] + pl_list[[15]] + pl_list[[16]] + pl_list[[17]] +
      pl_list[[18]] + pl_list[[19]] + pl_list[[20]] + pl_list[[21]] + pl_list[[22]] +
      plot_layout(widths = sides[[4]], nrow=1, ncol = 10)))
}


#' Prepare performance metrics table
#'
#' @param dt validated CNVs table
#'
#' @param unk treat unclear calls as unclear (3, excluded),
#'   true (1), or false (2)
#'
#' @export
#'
#' @import data.table

# prepare table for performance metrics plot

get_metrics_per_prob2 <- function(dt, unk = 3) {
  dt <- copy(dt)

  # drop UNK if wanted
  if (unk != 3) message('Unclear CNVs marked as: ', unk)
  dt[vo == 3, vo := unk]
  dt <- dt[vo %in% 1:2, ]

  if (!'prob' %in% colnames(dt)) {
    dt[GT == 1, prob := p_true_del]
    dt[GT == 2, prob := p_true_dup]
  }
  dt[vo == 1, human := T]
  dt[vo == 2, human := F]
  out <- data.table()

  for (i in seq(from = 0, to = 1, by = 0.01)) {
    tp <- dt[prob >= i & human == T, .N]
    tn <- dt[human == F & prob < i, .N]
    fp <- dt[prob >= i & human == F, .N]
    fn <- dt[human == T & prob < i, .N]
    out <- rbind(out, data.table(minprob = i, tp = tp, tn = tn, fp = fp, fn = fn,
                                 precision = signif(tp / (tp + fp), 4),
                                 recall = signif(tp / (tp + fn), 4)))
  }
  return(out)
}

#' Performance metrics plot
#'
#' @param dt output of `get_metrics_per_prob2()`
#'
#'
#' @export
#'
#' @import data.table
#' @import ggplot2

plot_performance <- function(dt) {
  a <- dt[, .(minprob, precision)]
  colnames(a) <- c('minprob', 'performance')
  b <- dt[, .(minprob, recall)]
  colnames(b) <- c('minprob', 'performance')
  a[, Metric := 'Precision / PPV']
  b[, Metric := 'Recall / Sensitivity']
  dt <- rbind(a, b)

  ggplot(dt) +
    geom_point(aes(x = minprob, y = performance, colour = Metric), size = 0.3) +
    geom_line(aes(x = minprob, y = performance, colour = Metric), alpha = 0.75) +
    scale_colour_manual(values = c('Precision / PPV' = '#0000CC', 'Recall / Sensitivity' = '#0066FF')) +
    theme_bw() + xlim(0, 1) + ylab('Performance') + xlab('Prediction Probability Cutoff') + 
    theme(legend.position = c(0.25, 0.15))

}
