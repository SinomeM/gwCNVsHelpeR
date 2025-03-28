% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{cnvs_distr_big}
\alias{cnvs_distr_big}
\title{Genomic distribution of CNVs}
\usage{
cnvs_distr_big(
  cnvs,
  chr_st_en,
  n_samples,
  bin_size = 250000,
  return_pl_list = F
)
}
\arguments{
\item{cnvs}{CNVs table}

\item{chr_st_en}{As usual}

\item{n_samples}{For the frequency calculation}

\item{bin_size}{Size of the unit}

\item{return_pl_list}{return all chromosomes individually
instead of four rows}
}
\description{
Plot the genomic distribution of deletion and duplication
genome-wide (chromosomes 1 to 22). The y axes is the rescaled
log frequency.
}
