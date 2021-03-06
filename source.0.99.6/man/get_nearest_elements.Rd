% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMI_met_cis_network.R
\name{get_nearest_elements}
\alias{get_nearest_elements}
\title{get_nearest_elements
This function is to get the neighbour elements on genome for target gene}
\usage{
get_nearest_elements(gene_list, distance = 3e+05, gene_number = 100,
  ref_gene_bed, ref_CpGs_bed, data_node_list)
}
\arguments{
\item{gene_list}{a character vector containing the list of target genes which to be tested}

\item{distance}{a numeric value determining the genome range of potential cis-acting network}

\item{gene_number}{a numeric value determining the max number of neighbour genes in this range}

\item{ref_gene_bed}{a data.frame containing reference gene coorinate with five columns named "name", "chr", "start", "end" and "strand". The coordinates of genes in exp_data_matrix are required to be included in this data.frame.}

\item{ref_CpGs_bed}{a data.frame containing reference CpGS coorinate with four columns names "name", "chr", "start" and "end". The coordinates of CpGs/probes in met_data_matrix are required to be included in this data.frame.}

\item{data_node_list}{a list of node names in data matrix. The return elements which are not in the node list will be excluded.}
}
\value{
a data frame containing the neighbour genes and neighbour CpGs
}
\description{
get_nearest_elements
This function is to get the neighbour elements on genome for target gene
}
\author{
Tong Yin
}
\keyword{internal}
