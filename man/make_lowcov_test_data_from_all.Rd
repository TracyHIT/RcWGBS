% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methylation_DNA_coding.R
\name{make_lowcov_test_data_from_all}
\alias{make_lowcov_test_data_from_all}
\title{Extract all sites' features with low coverage from cov file}
\usage{
make_lowcov_test_data_from_all(
  cov_filename = "data/Methyl/GM12878.0.5.strand_plus.cov",
  flank_size = 50,
  output_prefix = "data/Train_Test_Data/GM12878.0.5.lowcov",
  pre_DNA_inter_dir = "data/All_hg38_DNA/",
  file_type = ".gz",
  cov_high_threshold = 1000
)
}
\arguments{
\item{cov_filename}{The file contain "chr", "start", "end", "methylation", "m.read" and "um.read"}

\item{flank_size}{Side chain length, total window width 2*flank_size}

\item{output_prefix}{prefix of output file}

\item{pre_DNA_inter_dir}{The folder of output file}

\item{file_type}{cov_file type ".gz" or ".txt"}

\item{cov_high_threshold}{There are extremely abnormal high coverage sites during sequencing.
A fixed value can be set according to prior knowledge to filter out these sites.}
}
\value{
A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
}
\description{
Extract all sites' features with low coverage from cov file
}
