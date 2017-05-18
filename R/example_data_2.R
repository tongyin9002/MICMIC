#' methylation and expression data from TCGA HNSC dataset
#'
#' @docType data
#' @keywords datasets
#' @name TCGA_HNSC_data
#' @usage data(TCGA_HNSC_data)
#' @format matrix, data.frame
#' @return The methylation beta matrix, and expression log2 matrix
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network. "Comprehensive molecular characterization of gastric adenocarcinoma." Nature 513.7517 (2014): 202-209.


#' @docType data
#' @keywords datasets
#' @rdname  TCGA_HNSC_data
#' @details
#' HNSC_met_data_matrix: methylation values of CpGs in HNSC samples.
#' A data matrix containing the methylation beta values of CpGs in gastric cancer samples from TCGA
#' @usage HNSC_met_data_matrix
"HNSC_met_data_matrix"


#' @docType data
#' @keywords datasets
#' @rdname TCGA_HNSC_data
#' @details
#' HNSC_exp_data_matrix: expression values of genes in HNSC samples.
#' A data matrix containing the log2 transformed expression levels of genes in
#' gastric cancer samples from TCGA
#' @usage HNSC_exp_data_matrix
"HNSC_exp_data_matrix"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_HNSC_data
#' @details
#' HNSC_control_id: the sample ids of normal samples in TCGA HNSC dataset.
#' A vector containing sample ids
#' @usage HNSC_control_id
"HNSC_control_id"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_HNSC_data
#' @details
#' HNSC_tumor_id: the sample ids of tumor samples in TCGA HNSC dataset.
#' A vector containing sample ids
#' @usage HNSC_tumor_id
"HNSC_tumor_id"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_HNSC_data
#' @details
#' HNSC_ref_gene_bed: the gene coordinates in bed format.
#' A data.frame containing the chromosome coordinate information of all the genes
#' nearby the target genes
#' @usage HNSC_ref_gene_bed
#' @format A data.frame
#' \describe{
#'   \item{name}{name of genes}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of genes}
#'   \item{end}{the ending coordinate of genes}
#'   \item{strand}{the strand of genes}
#' }
#' @source \url{https://genome.ucsc.edu}
"HNSC_ref_gene_bed"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_HNSC_data
#' @details
#' HNSC_ref_CpGs_bed: the CpGs coordinates in bed format.
#' A data.frame containing the chromosome coordinate information of all the CpGs
#' nearby the target genes
#' @usage HNSC_ref_CpGs_bed
#' @format A data.frame
#' \describe{
#'   \item{name}{name of CpGs}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of CpGs}
#'   \item{end}{the ending coordinate of CpGs}
#' }
"HNSC_ref_CpGs_bed"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_HNSC_data
#' @details
#' HNSC_sample_class:
#' a data.frame containing the classes(control,tumor,stages or subtypes) information of samples
#' @usage HNSC_sample_class
#' @format A data.frame
#' \describe{
#'   \item{sample_id}{name of samples, should be exactly the same as the colnames of data matrix}
#'   \item{class}{control, tumor, stages, or tumor subtypes}
#' }
"HNSC_sample_class"


