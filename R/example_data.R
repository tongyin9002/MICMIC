#' methylation and expression data from TCGA STAD dataset
#'
#' @docType data
#' @keywords datasets
#' @name TCGA_STAD_data
#' @usage data(TCGA_STAD_data)
#' @format matrix, data.frame
#' @return The methylation beta matrix, and expression log2 matrix
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network. "Comprehensive molecular characterization of gastric adenocarcinoma." Nature 513.7517 (2014): 202-209.


#' @docType data
#' @keywords datasets
#' @rdname  TCGA_STAD_data
#' @details
#' STAD_met_data_matrix: methylation values of CpGs in STAD samples.
#' A data matrix containing the methylation beta values of CpGs in gastric cancer samples from TCGA
#' @usage STAD_met_data_matrix
"STAD_met_data_matrix"


#' @docType data
#' @keywords datasets
#' @rdname TCGA_STAD_data
#' @details
#' STAD_exp_data_matrix: expression values of genes in STAD samples.
#' A data matrix containing the log2 transformed expression levels of genes in
#' gastric cancer samples from TCGA
#' @usage STAD_exp_data_matrix
"STAD_exp_data_matrix"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_STAD_data
#' @details
#' STAD_control_id: the sample ids of normal samples in TCGA STAD dataset.
#' A vector containing sample ids
#' @usage STAD_control_id
"STAD_control_id"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_STAD_data
#' @details
#' STAD_tumor_id: the sample ids of tumor samples in TCGA STAD dataset.
#' A vector containing sample ids
#' @usage STAD_tumor_id
"STAD_tumor_id"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_STAD_data
#' @details
#' STAD_ref_gene_bed: the gene coordinates in bed format.
#' A data.frame containing the chromosome coordinate information of all the genes
#' nearby the target genes
#' @usage STAD_ref_gene_bed
#' @format A data.frame
#' \describe{
#'   \item{name}{name of genes}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of genes}
#'   \item{end}{the ending coordinate of genes}
#'   \item{strand}{the strand of genes}
#' }
#' @source \url{https://genome.ucsc.edu}
"STAD_ref_gene_bed"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_STAD_data
#' @details
#' STAD_ref_CpGs_bed: the CpGs coordinates in bed format.
#' A data.frame containing the chromosome coordinate information of all the CpGs
#' nearby the target genes
#' @usage STAD_ref_CpGs_bed
#' @format A data.frame
#' \describe{
#'   \item{name}{name of CpGs}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of CpGs}
#'   \item{end}{the ending coordinate of CpGs}
#' }
"STAD_ref_CpGs_bed"

#' @docType data
#' @keywords datasets
#' @rdname TCGA_STAD_data
#' @details
#' STAD_sample_class:
#' a data.frame containing the classes(control,tumor,stages or subtypes) information of samples
#' @usage STAD_sample_class
#' @format A data.frame
#' \describe{
#'   \item{sample_id}{name of samples, should be exactly the same as the colnames of data matrix}
#'   \item{class}{control, tumor, stages, or tumor subtypes}
#' }
"STAD_sample_class"















