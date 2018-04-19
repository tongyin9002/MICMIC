

#' example methylation and expression data from TCGA LIHC dataset
#' containing elements neighbouring HDAC11 gene
#' @docType data
#' @keywords datasets
#' @name LIHC
#' @usage data(LIHC)
#' @format matrix, data.frame
#' @return The methylation beta matrix, expression log2 matrix and genome coordinates of genes
#' and CpG sites.
#' @source \url{http://cancergenome.nih.gov/}
#' @references
#' Cancer Genome Atlas Research Network.

#' The annotation data frame of CpG sites used in LIHC example
#' @docType data
#' @keywords datasets
#' @name annotation
#' @usage data(annotation)
#' @format data.frame
#' @return The annotation data frame of CpG sites used in LIHC example.



#' The histone mark score of CpG sites used in LIHC example
#' @docType data
#' @keywords datasets
#' @name histone_mark
#' @usage data(histone_mark)
#' @format data.frame
#' @return CpGs_score_H3K4me1, CpGs_score_p300. The histone score of enhancer marks
#' H3K4me1 and p300 surrounding CpGs used in LIHC example.
#' @source \url{https://www.encodeproject.org/}



#' The histone mark score of CpG sites used in LIHC example
#' @docType data
#' @keywords datasets
#' @name TF_data
#' @usage data(TF_data)
#' @format data.frame
#' @return motif_data, CpG IDs neighbouring (up/down 250bp) TF binding motifs.
#' CpG_target, 10,000 example CpG-target pairs identified by MICMIC from LIHC.
#' control_list, 10,000 randomly selected CpG sites in genome as control list.
#' CpGs_in_SE, two column data.frame annotating CpGs in super enhancer.
#' function_gmt, a list of kegg pathway gene sets for functional study for TF targets.
#' @source \url{https://www.encodeproject.org/}



#' example methylation and expression data from TCGA STAD dataset
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







