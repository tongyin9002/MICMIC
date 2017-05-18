

#' generate_regulator_info
#'
#' This function is to integrate the regulator information in the gene_regulator_info.txt file
#' @param met_data_matrix a SummarizedExperiment object or a numeric matrix containing CpGs methylation data where columns contain samples and rows contain variables(probe site)
#' @param exp_data_matrix a SummarizedExperiment object or a numeric matrix containing gene expression data where columns contain samples and rows contain variables(gene site)
#' @param gene_list a vector containing the names of target genes
#' @param outfiledir a string of file directory to store the result files. If the parameter is not specified, the log file directory will be get by \code{getwd()}.
#' @param ref_gene_bed a GRanges object or a data.frame containing reference gene coorinate with five columns named "name", "chr", "start", "end" and "strand". The coordinates of genes in exp_data_matrix are required to be included in this data.frame.
#' @param ref_CpGs_bed a GRanges object or a data.frame containing reference CpGS coorinate with four columns names "name", "chr", "start" and "end". The coordinates of CpGs/probes in met_data_matrix are required to be included in this data.frame.
#' @usage generate_regulator_info(met_data_matrix,exp_data_matrix,gene_list,
#' outfiledir=NA,ref_gene_bed,ref_CpGs_bed)
#' @return data.frame containing information of direct and indirect regulators
#' @author Tong Yin
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom stats cor.test
#' @importFrom stats shapiro.test
#' @export
#' @import parallel
#' @examples
#'
#' data("TCGA_STAD_data")
#' gene_name<-"MLH1"
#' \dontrun{
#' generate_regulator_info(met_data_matrix=STAD_met_data_matrix,
#' exp_data_matrix=STAD_exp_data_matrix,gene_list=gene_name,
#' ref_gene_bed=STAD_ref_gene_bed,ref_CpGs_bed=STAD_ref_CpGs_bed)
#' }
#' @seealso the usage of CMI_met_cis_network


generate_regulator_info<-function(met_data_matrix,exp_data_matrix,gene_list,outfiledir=NA,ref_gene_bed,ref_CpGs_bed)
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


    ref_gene_list<-ref_gene_bed[,"name"]
    ref_CpGs_list<-ref_CpGs_bed[,"name"]

    data_matrix<-rbind(met_data_matrix,exp_data_matrix)

    cat("##############\n\n")
    for(i in seq_len(length(gene_list)))
    {
        cat("\n###Generation regulator information for ");cat(gene_list[i]);cat(" ");cat(i);cat(" in ");cat(length(gene_list));cat("\n")

        if(file.exists(paste(outfiledir,"/",gene_list[i],"_direct_regulation.txt",sep="")))
        {
         cat("read in regulation information...\n")
         association<-read.table(paste(outfiledir,"/",gene_list[i],"_associate_regulation.txt",sep=""),sep="\t")
         direct_regulation<-read.table(paste(outfiledir,"/",gene_list[i],"_direct_regulation.txt",sep=""),sep="\t")
         indirect_regulation<-read.table(paste(outfiledir,"/",gene_list[i],"_indirect_regulation.txt",sep=""),sep="\t")

         if(association[1,1]=="NodeA"){association<-association[-1,]}
         if(nrow(direct_regulation)==0){direct_regulation<-data.frame("","");direct_regulation<-direct_regulation[-1,]}
         if((indirect_regulation[1,1]=="NodeA")||(indirect_regulation[1,1]=="V1")||(indirect_regulation[1,1]=="X.NodeA.")){indirect_regulation<-indirect_regulation[-1,]}

         swap_lines<-which(association[,2]==gene_list[i])
         if(length(swap_lines)>0)
         {association[swap_lines,2]<-association[swap_lines,1]
          association[swap_lines,1]<-gene_list[i]
         }

         swap_lines<-which(direct_regulation[,2]==gene_list[i])
         if(length(swap_lines)>0)
         {direct_regulation[swap_lines,2]<-direct_regulation[swap_lines,1]
          direct_regulation[swap_lines,1]<-gene_list[i]
         }

         swap_lines<-which(indirect_regulation[,2]==gene_list[i])
         if(length(swap_lines)>0)
         {indirect_regulation[swap_lines,2]<-indirect_regulation[swap_lines,1]
          indirect_regulation[swap_lines,1]<-gene_list[i]
         }

         cat("generate target gene information\n")

         association_sites<-association[,2];
         association_sites<-association_sites[-which(association_sites%in%direct_regulation[,2])]
         association_sites<-association_sites[-which(association_sites%in%indirect_regulation[,2])]

         regulator_name<-c(association_sites,direct_regulation[,2],indirect_regulation[,2])

         target_name<-rep(gene_list[i],length(regulator_name))
         target_chr<-rep(ref_gene_bed[gene_list[i],"chr"],length(regulator_name))
         target_strand<-rep(ref_gene_bed[gene_list[i],"strand"],length(regulator_name))
         if(ref_gene_bed[gene_list[i],"strand"]=="+")
         {
            target_TSS<-ref_gene_bed[gene_list[i],"start"]
         }
         if(ref_gene_bed[gene_list[i],"strand"]=="-")
         {
            target_TSS<-ref_gene_bed[gene_list[i],"end"]
         }

         target_gene_length<-rep(abs(ref_gene_bed[gene_list[i],"start"]-ref_gene_bed[gene_list[i],"end"]),length(regulator_name))

         cat("generate regulator information\n")

         regulator_type<-c(rep("association",length(association_sites)),rep("direct",nrow(direct_regulation)),rep("indirect",nrow(indirect_regulation)))
         elements_type<-rep("CpGs_met",length(regulator_name))
         elements_type[which(regulator_name%in%ref_gene_list)]="gene_exp"

         regulator_chr<-rep("chr0",length(regulator_name))
         regulator_chr[which(elements_type=="gene_exp")]<-ref_gene_bed[regulator_name[which(elements_type=="gene_exp")],"chr"]
         regulator_chr[which(elements_type=="CpGs_met")]<-ref_CpGs_bed[regulator_name[which(elements_type=="CpGs_met")],"chr"]

         regulator_start<-rep(0,length(regulator_name))
         regulator_start[which(elements_type=="gene_exp")]<-ref_gene_bed[regulator_name[which(elements_type=="gene_exp")],"start"]
         regulator_start[which(elements_type=="CpGs_met")]<-ref_CpGs_bed[regulator_name[which(elements_type=="CpGs_met")],"start"]

         regulator_end<-rep(0,length(regulator_name))
         regulator_end[which(elements_type=="gene_exp")]<-ref_gene_bed[regulator_name[which(elements_type=="gene_exp")],"end"]
         regulator_end[which(elements_type=="CpGs_met")]<-ref_CpGs_bed[regulator_name[which(elements_type=="CpGs_met")],"end"]

         indirect_passenger<-rep("NA",length(regulator_name))
         indirect_passenger[which(regulator_type=="indirect")]<-unlist(lapply(rownames(indirect_regulation),function(x){strsplit(x," ")[[1]][3]}))

         cat("locate regulator...\n")

         regulation_dis<-rep("promoter",length(regulator_name))
         if(ref_gene_bed[gene_list[i],"strand"]=="+")
         {
            proximal_start<-ref_gene_bed[gene_list[i],"start"]-2000
            proximal_end<-ref_gene_bed[gene_list[i],"end"]
         }
         if(ref_gene_bed[gene_list[i],"strand"]=="-")
         {
            proximal_start<-ref_gene_bed[gene_list[i],"start"]
            proximal_end<-ref_gene_bed[gene_list[i],"end"]+2000
         }
         regulation_dis[which((regulator_start<proximal_start)|(regulator_start>proximal_end))]="distal"

         regulation_distance<-regulator_start-target_TSS

         if(ref_gene_bed[gene_list[i],"strand"]=="+")
         {
            regulation_dis[which((regulation_dis=="promoter")&(regulation_distance>500))]="gene_body"
         }
         if(ref_gene_bed[gene_list[i],"strand"]=="-")
         {
            regulation_dis[which((regulation_dis=="promoter")&(regulation_distance<(-500)))]="gene_body"
         }

         cat("calculate pearson correlation, standard deviation...\n")


         pearson_cor<-vector()
         pvalue<-vector()
         met_sd<-vector()
         gene_exp<-data_matrix[gene_list[i],]
         for(j in seq_len(length(regulator_name)))
         {cor_result<-cor.test(data_matrix[regulator_name[j],],gene_exp,method = "pearson")
          pearson_cor<-c(pearson_cor,round(cor_result$estimate,3))
          pvalue<-c(pvalue,cor_result$p.value)
         }

         met_sd<-apply(data_matrix[regulator_name,],1,sd)

         regulator_information<-data.frame(target_name,target_chr,target_strand,target_TSS,target_gene_length,regulator_name,regulator_type,elements_type,regulator_chr,regulator_start,regulator_end,regulation_dis,regulation_distance,met_sd,pearson_cor,pvalue,indirect_passenger)
         regulator_information<-regulator_information[order(regulator_start),]

         cat("delete nearby genes from regulators...\n")

         if(length(which(regulator_information[,"elements_type"]=="gene_exp"))>0)
         {
            regulator_information<-regulator_information[-which(regulator_information[,"elements_type"]=="gene_exp"),]
         }

         cat("write results... \n")
         write.table(regulator_information,paste(outfiledir,"/",gene_list[i],"_regulator_info.txt",sep=""),sep="\t",quote=FALSE)

        }else
        {
         cat("file not exist\n")

        }
    }
    cat("\n\n")

    return(regulator_information)

}













#' merge_regulator_info
#'
#' This function is to merge regulation information for multiple genes
#' @param gene_list a vector containing the names of target genes
#' @param outfiledir a string of file directory to store the result files. If the parameter is not specified, the log file directory will be get by \code{getwd()}.
#' @param statisticfiledir summary directory to store merged result. If the parameter is not specified, the file directory will be get by \code{getwd()}.
#' @param ref_gene_bed a data.frame containing reference gene coorinate with five columns named "name", "chr", "start", "end" and "strand". The coordinates of genes in exp_data_matrix are required to be included in this data.frame.
#' @usage merge_regulator_info(gene_list,outfiledir=NA,statisticfiledir=NA,ref_gene_bed)
#' @return numbers of direct and indirect regulators
#' @author Tong Yin
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export
#' @examples
#' \dontrun{
#' data("TCGA_LUAD_data")
#' gene_list<-rownames(LUAD_exp_data_matrix)[41:50]
#'
#' network<-CMI_met_cis_network(met_data_matrix=LUAD_met_data_matrix,
#' exp_data_matrix=LUAD_exp_data_matrix,gene_list=gene_list,distance=300000,
#' ref_gene_bed=LUAD_ref_gene_bed,ref_CpGs_bed=LUAD_ref_CpGs_bed,
#' core_num=1,permutation_times=20)
#'
#' generate_regulator_info(met_data_matrix=LUAD_met_data_matrix,
#' exp_data_matrix=LUAD_exp_data_matrix,gene_list=gene_list,
#' ref_gene_bed=LUAD_ref_gene_bed,ref_CpGs_bed=LUAD_ref_CpGs_bed)
#'
#' merge_regulator_info(gene_list=gene_list,ref_gene_bed=LUAD_ref_gene_bed)
#' }

merge_regulator_info<-function(gene_list,outfiledir=NA,statisticfiledir=NA,ref_gene_bed)
{
    options(stringsAsFactors = FALSE)

    if(is.na(outfiledir))
    {
        outfiledir<-getwd()
    }

    if(is.na(statisticfiledir))
    {
        statisticfiledir<-getwd()
    }

    direct_regulator_information<-data.frame()
    indirect_regulator_information<-data.frame()

    for(i in seq_len(length(gene_list)))
    {
        cat(i);cat("\t");cat(gene_list[i]);cat("\n");
        target_gene<-gene_list[i]
        if(file.exists(paste(outfiledir,"/",target_gene,"_regulator_info.txt",sep="")))
        {
         regulator_information<-read.table(paste(outfiledir,"/",target_gene,"_regulator_info.txt",sep=""),sep="\t")

         if(regulator_information[1,1]!="target_name")
         {
          direct_regulator<-regulator_information[which(regulator_information[,"regulator_type"]=="direct"),]
          indirect_regulator<-regulator_information[which(!(regulator_information[,"regulator_type"]=="direct")),]

          direct_regulator_information<-rbind(direct_regulator_information,direct_regulator)
          indirect_regulator_information<-rbind(indirect_regulator_information,indirect_regulator)
         }
        }
    }

    direct_regulator_bed<-direct_regulator_information[,c("regulator_chr","regulator_start","regulator_end","regulator_name")]
    indirect_regulator_bed<-indirect_regulator_information[,c("regulator_chr","regulator_start","regulator_end","regulator_name")]
    colnames(direct_regulator_bed)<-c("chr","start","end","name")
    colnames(indirect_regulator_bed)<-c("chr","start","end","name")

    direct_regulator_bed<-direct_regulator_bed[which(!duplicated(direct_regulator_bed[,"name"])),]
    indirect_regulator_bed<-indirect_regulator_bed[which(!duplicated(indirect_regulator_bed[,"name"])),]

    TSSpos<-ref_gene_bed[,"start"]
    minus_strand_lines<-which(ref_gene_bed[,"strand"]=="-")
    TSSpos[minus_strand_lines]<-ref_gene_bed[minus_strand_lines,"end"]
    promoter_start<-TSSpos-2000
    promoter_end<-TSSpos+500
    promoter_start[minus_strand_lines]<-TSSpos[minus_strand_lines]-500
    promoter_end[minus_strand_lines]<-TSSpos[minus_strand_lines]+2000

    ref_gene_promoter_bed<-data.frame(ref_gene_bed[,"chr"],promoter_start,promoter_end,ref_gene_bed[,"name"],ref_gene_bed[,"strand"])
    colnames(ref_gene_promoter_bed)<-c("chr","start","end","name","strand")
    ref_gene_promoter_bed<-ref_gene_promoter_bed[order(ref_gene_promoter_bed[,"start"]),]
    ref_gene_promoter_bed<-ref_gene_promoter_bed[order(ref_gene_promoter_bed[,"chr"]),]
    #write.table(ref_gene_promoter_bed,paste(statisticfiledir,"/ref_promoter_bed.txt",sep=""),quote=FALSE,sep="\t")

    cat("Annotating promoter");cat(".......\n")

    names<-c(colnames(direct_regulator_information),"in_promoter_region")
    overlapping_list<-bedoverlapping(direct_regulator_bed,ref_gene_promoter_bed)
    direct_regulator_information<-cbind(direct_regulator_information,direct_regulator_information[,"regulator_name"]%in%overlapping_list)
    colnames(direct_regulator_information)<-names

    names<-c(colnames(indirect_regulator_information),"in_promoter_region")
    overlapping_list<-bedoverlapping(indirect_regulator_bed,ref_gene_promoter_bed)
    indirect_regulator_information<-cbind(indirect_regulator_information,indirect_regulator_information[,"regulator_name"]%in%overlapping_list)
    colnames(indirect_regulator_information)<-names

    direct_regulator_information<-direct_regulator_information[-which((direct_regulator_information[,"regulation_dis"]=="distal")&(direct_regulator_information[,"in_promoter_region"])),]
    indirect_regulator_information<-indirect_regulator_information[-which((indirect_regulator_information[,"regulation_dis"]=="distal")&(indirect_regulator_information[,"in_promoter_region"])),]

    write.table(direct_regulator_information,paste(statisticfiledir,"/direct_information.txt",sep=""),quote=FALSE,sep="\t")
    write.table(indirect_regulator_information,paste(statisticfiledir,"/indirect_information.txt",sep=""),quote=FALSE,sep="\t")

    direct_target_gene<-unique(direct_regulator_information[,"target_name"])
    direct_promoter_regulator<-direct_regulator_information[which(direct_regulator_information[,"regulation_dis"]=="promoter"),"regulator_name"]
    direct_genebody_regulator<-unique(direct_regulator_information[which(direct_regulator_information[,"regulation_dis"]=="gene_body"),"regulator_name"])
    direct_distal_regulator<-unique(direct_regulator_information[which((direct_regulator_information[,"regulation_dis"]=="distal")&(!direct_regulator_information[,"in_promoter_region"])),"regulator_name"])

    indirect_target_gene<-unique(indirect_regulator_information[,"target_name"])
    indirect_promoter_regulator<-indirect_regulator_information[which(indirect_regulator_information[,"regulation_dis"]=="promoter"),"regulator_name"]
    indirect_genebody_regulator<-unique(indirect_regulator_information[which(indirect_regulator_information[,"regulation_dis"]=="gene_body"),"regulator_name"])
    indirect_distal_regulator<-unique(indirect_regulator_information[which((indirect_regulator_information[,"regulation_dis"]=="distal")&(!indirect_regulator_information[,"in_promoter_region"])),"regulator_name"])

    direct_num<-c(length(direct_target_gene),length(direct_promoter_regulator),length(direct_genebody_regulator),length(direct_distal_regulator))
    indirect_num<-c(length(indirect_target_gene),length(indirect_promoter_regulator),length(indirect_genebody_regulator),length(indirect_distal_regulator))
    regulation_num<-cbind(direct_num,indirect_num)
    rownames(regulation_num)<-c("target_gene_num","promoter_regulator_num","genebody_regulator_num","distal_regulator_num")

    write.table(regulation_num,paste(statisticfiledir,"/regulation_number.txt",sep=""),quote=FALSE,sep="\t")

    return(regulation_num)
}









bedoverlapping<-function(mybed,target_bed)
{

    mybed<-mybed[order(mybed[,"start"]),]
    mybed<-mybed[order(mybed[,"chr"]),]
    target_bed<-target_bed[order(target_bed[,"start"]),]
    target_bed<-target_bed[order(target_bed[,"chr"]),]

    bed1<-mybed[,c("chr","start","name")]
    bed2<-cbind(target_bed[,c("chr","start")],rep("targetstart",nrow(target_bed)))
    bed3<-cbind(target_bed[,c("chr","end")],rep("targetend",nrow(target_bed)))

    colnames(bed1)<-c("chr","pos","type")
    colnames(bed2)<-c("chr","pos","type")
    colnames(bed3)<-c("chr","pos","type")

    merge_bed<-rbind(bed1,bed2,bed3)

    merge_bed<-merge_bed[order(merge_bed[,"pos"]),]
    merge_bed<-merge_bed[order(merge_bed[,"chr"]),]

    mylines<-which((merge_bed[,"type"]!="targetstart")&(merge_bed[,"type"]!="targetend"))

    overlapping_list<-vector()
    start_num<-0
    for(i in seq_len(nrow(merge_bed)))
    {
     if(merge_bed[i,"type"]=="targetstart")
     {
         start_num<-start_num+1
     }
     if(merge_bed[i,"type"]=="targetend")
     {
         start_num<-start_num-1
     }
     if(i %in% mylines)
     {
      if(start_num>0)
      {
          overlapping_list<-c(overlapping_list,merge_bed[i,"type"])
      }
     }

    }

    return(overlapping_list)
}








functional_annotation_for_distal_regulators<-function(gene_list,outfiledir,annotation_bed_list,annotation_bed_names,statisticfiledir)
{
    ##read in regulator information
    direct_regulator_information<-read.table(paste(statisticfiledir,"/direct_information.txt",sep=""),sep="\t")
    indirect_regulator_information<-read.table(paste(statisticfiledir,"/indirect_information.txt",sep=""),sep="\t")

    result1<-direct_regulator_information[which(direct_regulator_information[,"regulation_dis"]!="promoter"),]

    result2<-indirect_regulator_information[which(indirect_regulator_information[,"regulation_dis"]!="promoter"),]

    result2<-result2[which(!result2[,"regulator_name"]%in%result1[,"regulator_name"]),]

    direct_regulator_bed<-data.frame(result1[,c("regulator_chr","regulator_start","regulator_end","regulator_name")])
    indirect_regulator_bed<-data.frame(result2[,c("regulator_chr","regulator_start","regulator_end","regulator_name")])

    colnames(direct_regulator_bed)<-c("chr","start","end","name")
    colnames(indirect_regulator_bed)<-c("chr","start","end","name")

    direct_regulator_bed<-direct_regulator_bed[which(!duplicated(direct_regulator_bed[,"name"])),]
    indirect_regulator_bed<-indirect_regulator_bed[which(!duplicated(indirect_regulator_bed[,"name"])),]

    direct_regulator_bed<-direct_regulator_bed[order(direct_regulator_bed[,"start"]),]
    direct_regulator_bed<-direct_regulator_bed[order(direct_regulator_bed[,"chr"]),]

    indirect_regulator_bed<-indirect_regulator_bed[order(indirect_regulator_bed[,"start"]),]
    indirect_regulator_bed<-indirect_regulator_bed[order(indirect_regulator_bed[,"chr"]),]

    direct_regulator_annotation<-direct_regulator_bed
    indirect_regulator_annotation<-indirect_regulator_bed
    for(i in seq_len(length(annotation_bed_list)))
    {
        cat("Annotating ");cat(annotation_bed_names[i]);cat(".......\n")
        overlapping_list<-bedoverlapping(direct_regulator_bed,annotation_bed_list[[i]])
        direct_regulator_annotation<-data.frame(direct_regulator_annotation,direct_regulator_bed[,"name"]%in%overlapping_list)
        overlapping_list<-bedoverlapping(indirect_regulator_bed,annotation_bed_list[[i]])
        indirect_regulator_annotation<-data.frame(indirect_regulator_annotation,indirect_regulator_bed[,"name"]%in%overlapping_list)
    }

    colnames(direct_regulator_annotation)<-c("chr","start","end","name",annotation_bed_names)
    colnames(indirect_regulator_annotation)<-c("chr","start","end","name",annotation_bed_names)
    rownames(direct_regulator_annotation)<-direct_regulator_annotation[,"name"]
    rownames(indirect_regulator_annotation)<-indirect_regulator_annotation[,"name"]

    write.table(direct_regulator_annotation,paste(statisticfiledir,"/direct_regulator_bed_annotation.txt",sep=""),sep="\t",quote=FALSE)
    write.table(indirect_regulator_annotation,paste(statisticfiledir,"/indirect_regulator_bed_annotation.txt",sep=""),sep="\t",quote=FALSE)

}


















