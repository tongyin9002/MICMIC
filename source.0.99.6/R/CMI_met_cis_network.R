################MICMIC#############################################
###Copyright (C) 2006-2018 Tong Yin <tongyin9002@gmail.com>########



#' Conditional mutual information learning the methylation cis-acting regulation network
#'
#' This function is to infer the cis-acting regulatory network between DNA methylation and gene expression
#'
#' @param met_data_matrix a SummarizedExperiment object or a numeric matrix containing CpGs methylation data where columns contain samples and rows contain variables(probe site)
#' @param exp_data_matrix a SummarizedExperiment object or a numeric matrix containing gene expression data where columns contain samples and rows contain variables(gene site)
#' @param gene_list a vector containing the names of target genes
#' @param distance integer specifying the upstream/downstream genome range to be analyzed
#' @param ref_gene_bed a GRanges object or a data.frame containing reference gene coorinate with five columns named "name", "chr", "start", "end" and "strand". The coordinates of genes in exp_data_matrix are required to be included in this data.frame.
#' @param ref_CpGs_bed a GRanges object or a data.frame containing reference CpGS coorinate with four columns names "name", "chr", "start" and "end". The coordinates of CpGs/probes in met_data_matrix are required to be included in this data.frame.
#' @param outfiledir a string of file directory to store the result files. If the parameter is not specified, the log file directory will be get by \code{getwd()}.
#' @param pvalue_cut the cutoff of pvalue. The default is 0.01.
#' @param core_num the cpu number using for parallel computation in PC_para
#' @param permutation_times the number of times of permutation to calculate the pvalue
#' @param MI_cut the mutual information cut off to filter the direct regulator
#' @param CpGrange CpG sites within this given range are supposed to be connected before MI/CMI testing
#' @usage CMI_met_cis_network(met_data_matrix,exp_data_matrix,gene_list,distance=300000,
#' ref_gene_bed,ref_CpGs_bed,outfiledir=NA,pvalue_cut=0.001,core_num=1,permutation_times=100)
#' @return the adjacency matrix of the network with value of 0 and 1. 1 means that there is an edge between the rowname and colname of the element.
#' And 0 means there is no edge.
#' @author Tong Yin
#' @export
#' @import parallel
#' @importFrom utils read.table
#' @importFrom SummarizedExperiment assay
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges strand
#' @examples
#' #library(MICMIC)
#' library(SummarizedExperiment)
#' library(GenomicRanges)
#' data("TCGA_STAD_data")
#' gene_name<-"MLH1"
#'
#' ####prepare the data matrix#####
#' colData<-DataFrame(sample_stages=STAD_sample_class[colnames(STAD_met_data_matrix),"stage"],
#' row.names=colnames(STAD_met_data_matrix))
#' STAD_met_data_SE <- SummarizedExperiment(assays=SimpleList(counts=STAD_met_data_matrix),
#'                                          colData=colData)
#' STAD_exp_data_SE <- SummarizedExperiment(assays=SimpleList(counts=STAD_exp_data_matrix),
#'                                          colData=colData)
#' STAD_ref_gene_bed_GR <-  GRanges(seqnames = Rle(STAD_ref_gene_bed[,"chr"]),
#'                                  ranges = IRanges(STAD_ref_gene_bed[,"start"],
#'                                  end=STAD_ref_gene_bed[,"end"],
#'                                  names = STAD_ref_gene_bed[,"name"]),
#'                                  strand = Rle(STAD_ref_gene_bed[,"strand"])
#' )
#' STAD_ref_CpGs_bed_GR <-  GRanges(seqnames = Rle(STAD_ref_CpGs_bed[,"chr"]),
#'                                  ranges = IRanges(STAD_ref_CpGs_bed[,"start"],
#'                                  end=STAD_ref_CpGs_bed[,"end"],
#'                                  names = STAD_ref_CpGs_bed[,"name"])
#' )
#' STAD_met_data_SE
#' STAD_exp_data_SE
#' STAD_ref_gene_bed_GR
#' STAD_ref_CpGs_bed_GR
#'
#' ####infer the cis-acting methylation regulatory network
#' network<-CMI_met_cis_network(met_data_matrix=STAD_met_data_SE,
#' exp_data_matrix=STAD_exp_data_SE,
#' gene_list=gene_name,distance=300000,ref_gene_bed=STAD_ref_gene_bed_GR,
#' ref_CpGs_bed=STAD_ref_CpGs_bed_GR,outfiledir=tempdir(),
#' pvalue_cut=0.0001,core_num=1)
#'
#' network
#'
#' ####generate the regulator information table####
#' result<-generate_regulator_info(met_data_matrix=STAD_met_data_SE,
#' exp_data_matrix=STAD_exp_data_SE,gene_list=gene_name,outfiledir=tempdir(),
#' ref_gene_bed=STAD_ref_gene_bed_GR,ref_CpGs_bed=STAD_ref_CpGs_bed_GR)
#'
#' result
#'
#' ####plotting the cis-acting regulatory network####
#'
#' MICMIC_plotting(gene_name=gene_name,met_data_matrix=STAD_met_data_SE,
#' exp_data_matrix=STAD_exp_data_SE,control_id=STAD_control_id,
#' distance=300000,ref_gene_bed=STAD_ref_gene_bed_GR,
#' ref_CpGs_bed=STAD_ref_CpGs_bed_GR,sample_class=STAD_sample_class,outfiledir=tempdir())
#'
#' ####done####




CMI_met_cis_network<-function(met_data_matrix,exp_data_matrix,gene_list,distance=300000,ref_gene_bed,ref_CpGs_bed,outfiledir=NA,pvalue_cut=0.001,core_num=1,permutation_times=100,MI_cut=0.1,CpGrange=NA,edgemode=c("pvalue","MI"))
{
    options(stringsAsFactors = FALSE)

    if(is.na(outfiledir))
    {
        outfiledir<-getwd()
    }

    if(class(met_data_matrix)=="SummarizedExperiment")
    {
        met_data_matrix<-assay(met_data_matrix)
    }

    if(class(exp_data_matrix)=="SummarizedExperiment")
    {
        exp_data_matrix<-assay(exp_data_matrix)
    }

    if(class(ref_gene_bed)=="GRanges")
    {
        ref_gene_bed<-data.frame(as.vector(seqnames(ref_gene_bed)),
                          start(ranges(ref_gene_bed)),
                          end(ranges(ref_gene_bed)),
                          as.vector(strand(ref_gene_bed)),
                          names(ref_gene_bed),
                          width(ranges(ref_gene_bed)))
        colnames(ref_gene_bed)<-c("chr","start","end","strand","name","length")
        rownames(ref_gene_bed)<-ref_gene_bed[,"name"]
    }
    if(class(ref_CpGs_bed)=="GRanges")
    {
        ref_CpGs_bed<-data.frame(names(ref_CpGs_bed),
                                 as.vector(seqnames(ref_CpGs_bed)),
                                 start(ranges(ref_CpGs_bed)),
                                 end(ranges(ref_CpGs_bed))
                                )
        colnames(ref_CpGs_bed)<-c("name","chr","start","end")
        rownames(ref_CpGs_bed)<-ref_CpGs_bed[,"name"]
    }

    ref_gene_bed<-ref_gene_bed[order(ref_gene_bed[,"start"]),]
    ref_gene_bed<-ref_gene_bed[order(ref_gene_bed[,"chr"]),]
    ref_CpGs_bed<-ref_CpGs_bed[order(ref_CpGs_bed[,"start"]),]
    ref_CpGs_bed<-ref_CpGs_bed[order(ref_CpGs_bed[,"chr"]),]


 ####input inconsistency detection
    for(i in seq_len(length(gene_list)))
    {   if(!gene_list[i]%in%rownames(exp_data_matrix))
        {stop("Error! Your gene is not in your exp data matrix");return("Error! Your gene is not in your exp data matrix")
        }
        if(!gene_list[i]%in%rownames(ref_gene_bed))
        {stop("Error! Cannot find the position of your gene in the reference gene bed");return("Error! Cannot find the position of your gene in the reference gene bed")
        }
    }
    if(!identical(colnames(met_data_matrix),colnames(exp_data_matrix)))
    {stop("Error! The individuals in methylation data matrix is not same as in the expression data matrix!");return("Error! The individuals in methylation data matrix is not same as in the expression data matrix!")
    }

    data_matrix<-rbind(met_data_matrix,exp_data_matrix);data_node_list<-rownames(data_matrix)

 ####find neighbour elements on the genome for target gnes
    cat("\n");cat("######");cat("\n\n");cat("Start to discover neighbour elements for your genes within ");cat(distance);cat(" bp in the genome");cat("\n\n");cat("######");cat("\n\n")

    cat("\n##getting elements....\n")
    elements<-get_nearest_elements(gene_list=gene_list,distance=distance,gene_number=100,ref_gene_bed=ref_gene_bed,ref_CpGs_bed=ref_CpGs_bed,data_node_list=data_node_list)

    #cat("\n##selecting associated elements...\n")
    elements<-select_associated_sites(elements,ref_CpGs_bed,data_matrix,outfiledir,core_num=core_num)
    if(nrow(elements)==0){return();}

    cat("\n##making adjacent matrix...\n\n")
    adj_matrix<-make_cis_matrix(elements=elements,ref_CpGs_bed=ref_CpGs_bed,
                                CpGrange=CpGrange)

 ####slimming the data_matrix according to the neighbour elements
    elements_list<-rownames(adj_matrix);data_matrix<-data_matrix[elements_list,]
    gc()

    for(i in seq_len(nrow(elements)))
    {
        cis_genes<-strsplit(elements[i,2],",")[[1]];cis_CpGs<-strsplit(elements[i,3],",")[[1]]

        cat("The cis-acting regulatory network of your gene ");cat(elements[i,1]);cat(" will be constructed by ");
        cat(" expression level of ");cat(length(cis_genes));cat(" nearby genes and methylation level of ");cat(length(cis_CpGs));cat(" nearby CpG sites.");
        cat("\n\n")
    }

    cat("\n");cat("############");cat("\n\n");
    cat("The data matrix is ready, including ");cat(nrow(data_matrix));cat(" data rows for ");cat(ncol(data_matrix));cat(" individual columns.");
    cat("\n\n");cat("Importing data matrix into CMI-GRP...");cat("\n\n");cat("#######Starting network construction...########");cat("\n\n")

    gc()

 ####run the CMI-GRP: General regulatory patterns inference based on PC-algorithm by conditional mutual information
    CMI_PC_net<-PC_para(data_matrix,max_L=1,method="CMII",pre_adj=adj_matrix,log_file_dir=outfiledir,edgemode=edgemode,pvalue_cut=pvalue_cut,core_num=core_num,permutation_times=permutation_times)

    MI_net<-read.table(paste(outfiledir,"/adj_log0.txt",sep=""),sep="\t")
    colnames(MI_net)<-chartr(".", "-", colnames(MI_net))
    final_edges<-read.table(paste(outfiledir,"/final_edges.txt",sep=""),sep="\t")
    if(file.info(paste(outfiledir,"/del_log1.txt",sep=""))$size>1)
    {
     CMI_del_edges<-read.table(paste(outfiledir,"/del_log1.txt",sep=""),sep="\t")
    }else
    {
     CMI_del_edges<-data.frame("NodeA","NodeB");
    }
    cat("\n\n")

    MIresult<-read.table(paste(outfiledir,"/mi_log0.txt",sep=""),sep="\t")
    CMIresult<-read.table(paste(outfiledir,"/mi_log1.txt",sep=""),sep="\t")

    write.table(MIresult,paste(outfiledir,"/MI_result.txt",sep=""),sep="\t",quote=FALSE)
    write.table(CMIresult,paste(outfiledir,"/CMI_result.txt",sep=""),sep="\t",quote=FALSE)

    cat("####### Result #########")
    cat("\n")

    for(i in seq_len(nrow(elements)))
    {
        cat("\n")

        gene_name<-elements[i,1]
        cis_genes<-strsplit(elements[i,2],",")[[1]]
        cis_CpGs<-strsplit(elements[i,3],",")[[1]]
        elements_list<-c(cis_genes,cis_CpGs)

        regulatory_net<-CMI_PC_net[elements_list,elements_list]

        #### write regulatory_net ####

        write.table(MI_net[elements_list,elements_list],paste(outfiledir,"/",gene_name,"_MI_network.txt",sep=""),sep="\t",quote=FALSE)
        write.table(regulatory_net,paste(outfiledir,"/",gene_name,"_CMI_network.txt",sep=""),sep="\t",quote=FALSE)

        #### write direct and indirect regulation #####

        direct_edges<-final_edges[which((final_edges[,1]%in%elements_list)&(final_edges[,2]%in%elements_list)),]
        direct_edges<-direct_edges[which((direct_edges[,1]==gene_name)|(direct_edges[,2]==gene_name)),]

        if(nrow(direct_edges)>0)
        {
         mi<-unlist(mclapply(1:nrow(direct_edges),function(x){MI(data_matrix[direct_edges[x,1],],data_matrix[direct_edges[x,2],])},mc.cores=core_num))

         type<-rep("direct",length(mi))
         if(!is.na(MI_cut))
         {type[which(mi<MI_cut)]<-"weak_direct"
         }
         direct_edges<-data.frame(direct_edges,mi,type)

         write.table(direct_edges,paste(outfiledir,"/",gene_name,"_direct_regulation.txt",sep=""),sep="\t",quote=FALSE)
        }

        indirect_edges<-CMI_del_edges[which((CMI_del_edges[,1]%in%elements_list)&(CMI_del_edges[,2]%in%elements_list)),]

        if(nrow(indirect_edges)>0)
        {
         indirect_edges<-indirect_edges[which(indirect_edges[,1]==gene_name),]
         if(nrow(indirect_edges)>0)
         {
          write.table(indirect_edges,paste(outfiledir,"/",gene_name,"_indirect_regulation.txt",sep=""),sep="\t",quote=FALSE)
         }
        }

        cat("\n###")
        cat("Discovered ");cat(nrow(direct_edges));cat(" direct regulatory relationship in ");
        cat(gene_name);cat(" cis-acting network.");
        cat("\n")
        cat("Discovered ");cat(nrow(indirect_edges));cat(" indirect regulatory relationship in ");
        cat(gene_name);cat(" cis-acting network.");
        cat("\n")

        #### write cis regulator #####

        regulatory_net<-regulatory_net+t(regulatory_net)

        regulator_list<-colnames(regulatory_net)[which(regulatory_net[gene_name,]==1)]
        regulator_genes<-intersect(regulator_list,cis_genes)
        regulator_CpGs<-intersect(regulator_list,cis_CpGs)

        regulator_genes_info<-data.frame(ref_gene_bed[regulator_genes,],(ref_gene_bed[regulator_genes,"start"]-ref_gene_bed[gene_name,"start"]))
        regulator_CpGs_info<-data.frame(ref_CpGs_bed[regulator_CpGs,],(ref_CpGs_bed[regulator_CpGs,"start"]-ref_gene_bed[gene_name,"start"]))

        colnames(regulator_genes_info)<-c("name","chr","start","end","strand","distance")
        colnames(regulator_CpGs_info)<-c("name","chr","start","end","distance")

        #write.table(regulator_genes_info,paste(outfiledir,"/",gene_name,"_regulator_genes_info.txt",sep=""),sep="\t",quote=FALSE)
        #write.table(regulator_CpGs_info,paste(outfiledir,"/",gene_name,"_regulator_CpGs_info.txt",sep=""),sep="\t",quote=FALSE)

        cat("\n###")
        cat("\n")
        cat("Discovered ");cat(sum(regulatory_net[gene_name,]));cat(" regulators for ");cat(gene_name)
        cat(", including ");cat(nrow(regulator_genes_info));cat(" genes and ");cat(nrow(regulator_CpGs_info));cat(" CpGs.")
        cat("\n\n")
        cat("###")

    }
    return(CMI_PC_net)
}






#' get_nearest_elements
#' This function is to get the neighbour elements on genome for target gene
#'
#' @param gene_list a character vector containing the list of target genes which to be tested
#' @param distance a numeric value determining the genome range of potential cis-acting network
#' @param gene_number a numeric value determining the max number of neighbour genes in this range
#' @param ref_gene_bed a data.frame containing reference gene coorinate with five columns named "name", "chr", "start", "end" and "strand". The coordinates of genes in exp_data_matrix are required to be included in this data.frame.
#' @param ref_CpGs_bed a data.frame containing reference CpGS coorinate with four columns names "name", "chr", "start" and "end". The coordinates of CpGs/probes in met_data_matrix are required to be included in this data.frame.
#' @param data_node_list a list of node names in data matrix. The return elements which are not in the node list will be excluded.
#' @keywords internal
#' @return a data frame containing the neighbour genes and neighbour CpGs
#' @author Tong Yin

get_nearest_elements<-function(gene_list,distance=300000,gene_number=100,ref_gene_bed,ref_CpGs_bed,data_node_list)
{
    elements<-data.frame()

    for(i in seq_len(length(gene_list)))
    {
        gene_name <- gene_list[i]
        nearest_genes_list<-get_neareast_genes_for_one_gene(gene_name=gene_name,distance=distance,gene_number=gene_number,ref_gene_bed=ref_gene_bed)
        nearest_genes_list<-intersect(nearest_genes_list,data_node_list)
        nearest_genes<-paste0(nearest_genes_list,collapse = ",")

        gene1<-nearest_genes_list[1]
        gene2<-nearest_genes_list[length(nearest_genes_list)]
        nearest_CpGs_list<-get_CpGs_between_genes(gene1=gene1,gene2=gene2,CpGsrange=2000,ref_gene_bed=ref_gene_bed,ref_CpGs_bed=ref_CpGs_bed)
        nearest_CpGs_list<-intersect(nearest_CpGs_list,data_node_list)
        nearest_CpGs<-paste0(nearest_CpGs_list,collapse = ",")

        elements<-rbind(elements,c(gene_name,nearest_genes,nearest_CpGs))

    }

    colnames(elements)=c("gene_name","near_genes","near_CpGs")
    return(elements)
}








get_neareast_genes_for_one_gene<-function(gene_name,distance=300000,gene_number=100,ref_gene_bed)
{
    chr<-ref_gene_bed[gene_name,"chr"]
    coordinate<-as.numeric(ref_gene_bed[gene_name,"start"])
    gene_names<-ref_gene_bed[,"name"]

    pos<-which(gene_names==gene_name);pos1<-pos-gene_number;pos2<-pos+gene_number

    if(pos1<=0){pos1<-1}

    gene_names<-gene_names[pos1:pos2]

    gene_chrs<-ref_gene_bed[gene_names,"chr"]
    gene_names<-gene_names[which(gene_chrs==chr)]

    gene_coordinates<-as.numeric(ref_gene_bed[gene_names,"start"])
    gene_names<-gene_names[which(abs(gene_coordinates-coordinate)<distance)]

    return(unique(gene_names))
}








get_CpGs_between_genes<-function(gene1,gene2,CpGsrange=2000,ref_gene_bed,ref_CpGs_bed)
{
    chr1<-ref_gene_bed[gene1,"chr"];chr2<-ref_gene_bed[gene2,"chr"]

    if(chr1!=chr2)
    {stop("Chromosomes of this two genes are not the same")
    }

    pos1<-min(as.numeric(ref_gene_bed[gene1,"start"]),as.numeric(ref_gene_bed[gene1,"end"]))
    pos2<-max(as.numeric(ref_gene_bed[gene2,"start"]),as.numeric(ref_gene_bed[gene2,"end"]))

    if(pos1>pos2)
    {stop("position of geneA is larger than position of geneB")
    }

    pos1<-pos1-CpGsrange;pos2<-pos2+CpGsrange;chr<-chr1

    return(get_CpGs_in_range(chr,pos1,pos2,ref_CpGs_bed))
}








get_CpGs_in_range<-function(chr,range_start,range_end,ref_CpGs_bed)
{
    pos<-which((ref_CpGs_bed[,"chr"]==chr)&(as.numeric(ref_CpGs_bed[,"start"])>range_start)&(as.numeric(ref_CpGs_bed[,"end"])<range_end))
    CpGs_id<-ref_CpGs_bed[pos,1]
    return(unique(CpGs_id))
}




select_associated_sites<-function(elements,ref_CpGs_bed,data_matrix,outfiledir,core_num=1)
{
    cis_genes<-vector();cis_CpGs<-vector()
    ass_num<-vector();
    for(i in seq_len(nrow(elements)))
    {   #cat("select associated sites ");cat(i);cat(" in ");cat(nrow(elements));cat(" ...\n")
        gene_name<-elements[i,1];cis_genes<-strsplit(elements[i,2],",")[[1]];cis_CpGs<-strsplit(elements[i,3],",")[[1]]

        if(length(cis_CpGs)>0)
        {
         cortest<-mclapply(seq_len(length(cis_CpGs)),function(x){result<-cor.test(data_matrix[gene_name,],data_matrix[cis_CpGs[x],]);return(c(result$estimate,result$p.value))},mc.cores = core_num)
         cortestresult<-do.call(rbind,cortest)
         pvalues<-unlist(cortestresult[,2])
         cut_p_lines<-which(pvalues<0.01)
         cut_pvalues<-pvalues[cut_p_lines]
         associate_CpGs<-cis_CpGs[cut_p_lines]

         #cat("find ");cat(length(associate_CpGs));cat(" CpGs are associated with ");cat(gene_name);cat(" with pvalue < 0.01\n");

         regulation<-data.frame(rep(gene_name,length(associate_CpGs)),associate_CpGs,cut_pvalues)
         colnames(regulation)<-c("NodeA","NodeB","pvalue")
         write.table(regulation,paste(outfiledir,"/",gene_name,"_associate_regulation.txt",sep=""),sep="\t",quote=FALSE)
         ass_num<-c(ass_num,length(associate_CpGs))
         associate_CpGs<-paste0(associate_CpGs,collapse = ",")
         elements[i,3]<-associate_CpGs

         cortestresult<-data.frame(gene_name,cis_CpGs,cortestresult)
         colnames(cortestresult)<-c("gene","CpGs","cor","pvalue")
         write.table(cortestresult,paste(outfiledir,"/",gene_name,"_cor_test.txt",sep=""),sep="\t",quote=FALSE)

        }else
        {
         ass_num<-c(ass_num,0)
        }
    }
    elements<-data.frame(elements,ass_num)
    elements<-elements[which(elements[,"ass_num"]>1),]
    #elements<-elements[which(elements[,"ass_num"]<250),]

    return(elements)
}





make_cis_matrix<-function(elements,ref_CpGs_bed,CpGrange=NA)
{   if(is.na(CpGrange))
    {
     CpGrange=15000
    }
    near_genes<-vector();near_CpGs<-vector()
    for(i in seq_len(nrow(elements)))
    {
     near_genes<-c(near_genes,strsplit(elements[i,2],",")[[1]])
     near_CpGs<-c(near_CpGs,strsplit(elements[i,3],",")[[1]])
    }
    near_genes<-unique(near_genes);near_CpGs<-unique(near_CpGs)

    elements_number<-length(near_genes)+length(near_CpGs)
    adj_matrix<-matrix(0,nrow=elements_number,ncol=elements_number)
    colnames(adj_matrix)<-c(near_genes,near_CpGs)
    rownames(adj_matrix)<-c(near_genes,near_CpGs)

    cis_genes<-vector();cis_CpGs<-vector()
    myCpGbed<-ref_CpGs_bed[near_CpGs,]

    for(i in seq_len(nrow(elements)))
    {   cat("making adjacent matrix ");cat(i);cat(" in ");cat(nrow(elements));cat(" ...\n")
        cis_genes<-strsplit(elements[i,2],",")[[1]];cis_CpGs<-strsplit(elements[i,3],",")[[1]]

        #interactions between genes and genes, CpGs and genes are supposed to be true
        adj_matrix[cis_genes,cis_genes]=1
        adj_matrix[cis_genes,cis_CpGs]=1
        adj_matrix[cis_CpGs,cis_genes]=1

        #interactions between CpGs within given range are supposed to be true
        for(j in seq_len(length(cis_CpGs)))
        {
            CpG<-cis_CpGs[j]

            near_CpGs<-get_CpGs_in_range(myCpGbed[CpG,"chr"],myCpGbed[CpG,"start"]-CpGrange,myCpGbed[CpG,"start"]+CpGrange,myCpGbed)
            near_CpGs<-intersect(cis_CpGs,near_CpGs)
            if(length(near_CpGs)>0)
            {adj_matrix[CpG,near_CpGs]<-1
             adj_matrix[near_CpGs,CpG]<-1
            }
        }
    }
    return(adj_matrix)
}





















