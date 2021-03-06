% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motif_analysis.R
\name{LoadFastaData}
\alias{LoadFastaData}
\title{Load the SNP data from fasta files.}
\usage{
LoadFastaData(ref.filename = NULL, snp.filename = NULL,
  ref.urlname = NULL, snp.urlname = NULL, snpids = NULL,
  default.par = FALSE)
}
\arguments{
\item{ref.filename}{a fastq file name for the reference allele sequences.}

\item{snp.filename}{a fastq file name for the SNP allele sequences.}

\item{ref.urlname}{URL of a fastq file for the reference allele sequences.}

\item{snp.urlname}{URL of a fastq file for the SNP allele sequences.}

\item{snpids}{SNP IDs}

\item{default.par}{A boolean for whether using the default Markov parameters.
Default: FALSE.}
}
\value{
A list object containing the following components:
\tabular{ll}{
sequence_matrix \tab A list of integer vectors representing the deroxyribose 
sequence around each SNP.\cr
a1 \tab An integer vector for the deroxyribose at the SNP location on the 
reference genome.\cr
a2 \tab An integer vector for the deroxyribose at the SNP location on the SNP
 genome.\cr
}
The results are coded as: "A"-1, "C"-2, "G"-3, "T"-4.
}
\description{
Load SNP data.
}
\examples{
LoadFastaData(
ref.urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_1.fasta",
snp.urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_2.fasta")
}
\author{
Sunyoung Shin \email{sunyoung.shin@utdallas.edu}, Chandler Zuo 
\email{chandler.c.zuo@gmail.com}
}
