#' methylation values of CpGs in STAD samples
#'
#' A data matrix containing the methylation beta values of CpGs in gastric cancer
#' samples from TCGA
#'
#' @docType data
#' @keywords datasets
#' @name STAD_met_data_matrix
#' @usage data(TCGA_STAD_data)
#' @format A matrix of numeric values
#' @return The methylation beta matrix for gastric cancer samples
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network. "Comprehensive molecular characterization of gastric adenocarcinoma." Nature 513.7517 (2014): 202-209.
"STAD_met_data_matrix"



#' expression values of genes in STAD samples
#'
#' A data matrix containing the log2 transformed expression levels of genes in
#' gastric cancer samples from TCGA
#'
#' @docType data
#' @keywords datasets
#' @name STAD_exp_data_matrix
#' @usage data(TCGA_STAD_data)
#' @format A matrix of numeric values
#' @return The log2 expression matrix for gastric cancer samples
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network. "Comprehensive molecular characterization of gastric adenocarcinoma." Nature 513.7517 (2014): 202-209.
"STAD_exp_data_matrix"



#' the sample ids of normal samples in TCGA STAD dataset
#'
#' A vector containing sample ids
#'
#' @docType data
#' @keywords datasets
#' @name STAD_control_id
#' @usage data(TCGA_STAD_data)
#' @format A vector of character
#' @return A vector of sample ids
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network. "Comprehensive molecular characterization of gastric adenocarcinoma." Nature 513.7517 (2014): 202-209.
"STAD_control_id"



#' the sample ids of tumor samples in TCGA STAD dataset
#'
#' A vector containing sample ids
#'
#' @docType data
#' @keywords datasets
#' @name STAD_tumor_id
#' @usage data(TCGA_STAD_data)
#' @format A vector of character
#' @return A vector of sample ids
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network. "Comprehensive molecular characterization of gastric adenocarcinoma." Nature 513.7517 (2014): 202-209.
"STAD_tumor_id"






#' the gene coordinates in bed format
#'
#' A data.frame containing the chromosome coordinate information of all the genes
#' nearby the target genes
#'
#' @docType data
#' @keywords datasets
#' @name STAD_ref_gene_bed
#' @usage data(TCGA_STAD_data)
#' @format A data.frame
#' @return A data.frame containing genome coordinates of genes
#' \describe{
#'   \item{name}{name of genes}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of genes}
#'   \item{end}{the ending coordinate of genes}
#'   \item{strand}{the strand of genes}
#' }
#' @source \url{https://genome.ucsc.edu}
"STAD_ref_gene_bed"



#' the CpGs coordinates in bed format
#'
#' A data.frame containing the chromosome coordinate information of all the CpGs
#' nearby the target genes
#'
#' @docType data
#' @keywords datasets
#' @name STAD_ref_CpGs_bed
#' @usage data(TCGA_STAD_data)
#' @format A data.frame
#' @return A data.frame containing genome coordinates of CpGs
#' \describe{
#'   \item{name}{name of CpGs}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of CpGs}
#'   \item{end}{the ending coordinate of CpGs}
#' }
#' @source \url{https://genome.ucsc.edu}
"STAD_ref_CpGs_bed"




#' the classes of samples
#'
#' A data.frame containing the classes(control,tumor,stages or subtypes) information of samples
#'
#' @docType data
#' @keywords datasets
#' @name STAD_sample_class
#' @usage data(TCGA_STAD_data)
#' @format A data.frame
#' @return sample classification information in two column data.frame
#' \describe{
#'   \item{sample_id}{name of samples, should be exactly the same as the colnames of data matrix}
#'   \item{class}{control, tumor, stages, or tumor subtypes}
#' }
#' @source \url{https://genome.ucsc.edu}
"STAD_sample_class"

















#' methylation values of CpGs in HNSC samples
#'
#' A data matrix containing the methylation beta values of CpGs in head and neck cancer
#' samples from TCGA
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_met_data_matrix
#' @usage data(TCGA_HNSC_data)
#' @format A matrix of numeric values
#' @return The methylation beta matrix for head and neck cancer samples
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Network. "Comprehensive genomic characterization of head and neck squamous cell carcinomas." Nature 517.7536 (2015): 576-582.
"HNSC_met_data_matrix"



#' expression values of genes in HNSC samples
#'
#' A data matrix containing the log2 transformed expression levels of genes in
#' head and neck cancer samples from TCGA
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_exp_data_matrix
#' @usage data(TCGA_HNSC_data)
#' @format A matrix of numeric values
#' @return The log2 expression matrix for head and neck cancer samples
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Network. "Comprehensive genomic characterization of head and neck squamous cell carcinomas." Nature 517.7536 (2015): 576-582.
"HNSC_exp_data_matrix"



#' the sample ids of normal samples in TCGA HNSC dataset
#'
#' A vector containing sample ids
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_control_id
#' @usage data(TCGA_HNSC_data)
#' @format A vector of character
#' @return A vector of sample ids
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Network. "Comprehensive genomic characterization of head and neck squamous cell carcinomas." Nature 517.7536 (2015): 576-582.
"HNSC_control_id"



#' the sample ids of tumor samples in TCGA HNSC dataset
#'
#' A vector containing sample ids
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_tumor_id
#' @usage data(TCGA_HNSC_data)
#' @format A vector of character
#' @return A vector of sample ids
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Network. "Comprehensive genomic characterization of head and neck squamous cell carcinomas." Nature 517.7536 (2015): 576-582.
"HNSC_tumor_id"






#' the gene coordinates in bed format
#'
#' A data.frame containing the chromosome coordinate information of all the genes
#' nearby the target genes
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_ref_gene_bed
#' @usage data(TCGA_HNSC_data)
#' @format A data.frame
#' @return A data.frame containing genome coordinates of genes
#' \describe{
#'   \item{name}{name of genes}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of genes}
#'   \item{end}{the ending coordinate of genes}
#'   \item{strand}{the strand of genes}
#' }
#' @source \url{https://genome.ucsc.edu}
"HNSC_ref_gene_bed"



#' the CpGs coordinates in bed format
#'
#' A data.frame containing the chromosome coordinate information of all the CpGs
#' nearby the target genes
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_ref_CpGs_bed
#' @usage data(TCGA_HNSC_data)
#' @format A data.frame
#' @return A data.frame containing genome coordinates of CpGs
#' \describe{
#'   \item{name}{name of CpGs}
#'   \item{chr}{the chromosome id like chr1, chr2, chr3 ...}
#'   \item{start}{the starting coordinate of CpGs}
#'   \item{end}{the ending coordinate of CpGs}
#' }
#' @source \url{https://genome.ucsc.edu}
"HNSC_ref_CpGs_bed"





#' the classes of samples
#'
#' A data.frame containing the classes(control,tumor,stages or subtypes) information of samples
#'
#' @docType data
#' @keywords datasets
#' @name HNSC_sample_class
#' @usage data(TCGA_HNSC_data)
#' @format A data.frame
#' @return sample classification information in two column data.frame
#' \describe{
#'   \item{sample_id}{name of samples, should be exactly the same as the colnames of data matrix}
#'   \item{class}{control, tumor, stages, or tumor subtypes}
#' }
#' @source \url{https://genome.ucsc.edu}
"HNSC_sample_class"



