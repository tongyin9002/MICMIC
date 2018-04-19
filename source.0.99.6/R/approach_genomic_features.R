################MICMIC#############################################
###Copyright (C) 2006-2018 Tong Yin <tongyin9002@gmail.com>########



#' Plot of average histone mark score around CpGs detected by MICMIC
#'
#' @param direct_info Regulation information of CpG gene pairs detected by MICMIC.
#' Additional type column should be added to designate the type of regulation.
#' @param CpGs_annotation Annotation data frame of each CpG site.
#' @param name The name of histone mark.
#' @param CpGs_score The data matrix of histone score around each CpG site.
#' @param colors The plotting colors of types of DRE
#' @param range The plotting range (bp) around the CpG site.
#' @param figure_dir The output directory.
#' @usage plot_histone_mark_score(direct_info,name,CpGs_score,colors,range,figure_dir)
#' @return A plotting result would be saved in the output directory to a file named by the histone mark.
#' @export

plot_histone_mark_score<-function(direct_info,name,CpGs_score,colors,range,figure_dir,CpGs_annotation=NA)
{   options(stringsAsFactors = FALSE)
    if(!is.na(CpGs_annotation))
    {
     if("promoter"%in%colnames(CpGs_annotation))
     {direct_info<-direct_info[which(!CpGs_annotation[direct_info[,"regulator_name"],"promoter"]),]
      control_CpGs<-(CpGs_annotation[which((!CpGs_annotation[,"promoter"])),"name"])
     }else
     {
      control_CpGs<-(CpGs_annotation[,"name"])
     }
    }

    score_result<-data.frame()

    random_score<-apply(CpGs_score[control_CpGs,],2,mean)
    y<-data.frame(names(random_score),random_score,"control")
    colnames(y)<-c("coord","score","type")
    score_result<-rbind(score_result,y)

    types<-unique(direct_info$type)
    for(i in 1:length(types))
    {lines<-which(direct_info$type==types[i])
    direct_score<-apply(CpGs_score[direct_info[lines,"regulator_name"],],2,mean)
    x<-data.frame(names(direct_score),direct_score,types[i])
    colnames(x)<-c("coord","score","type")
    score_result<-rbind(score_result,x)
    }

    score_result[,"coord"]<-as.numeric(score_result[,"coord"])

    myscore_result<-score_result

    myplot <- ggplot()+
        geom_point(data=myscore_result,aes(x=coord,y=score,color = factor(type)),size=0.5,alpha=0.5)+
        geom_smooth(data=myscore_result,aes(x=coord,y=score,color = factor(type)),se=FALSE,span = 0.05)+
        scale_x_continuous(breaks=seq(-range,range,by=500))+
        scale_y_continuous()+
        scale_color_manual(values=mycolors)+
        xlab(paste0("distance to CpGs in (bp)"))+
        ylab("score")+
        theme_bw()+
        theme(axis.text.x=element_text(size=10,angle=90,colour="black"),
              axis.text.y=element_text(size=10,colour="black"),
              axis.title.x=element_text(size=10,face="bold"),
              axis.title.y=element_text(size=10,face="bold"),
              plot.title=element_text(size=10,face="bold"),
              legend.title=element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank()
        )

    ggsave(myplot,file=paste0(figure_dir,"/",name,".pdf"),width=4,height=3)
}





#' Count the number of DRE overlapping with chromatin states
#'
#' @param direct_info Regulation information of CpG gene pairs detected by MICMIC.
#' Additional type column should be added to designate the type of regulation.
#' @param CpGs_annotation Annotation data frame of each CpG site.
#' @param chromcol The column number of chromatin state data in CpGs_annotation data.frame.
#' @usage hmm_overlapping(direct_info,CpGs_annotation,chromcol=9)
#' @return The counting of CpGs located in each chromatin state
#' @export

hmm_overlapping<-function(direct_info,CpGs_annotation,chromcol=9)
{   options(stringsAsFactors = FALSE)

    if("promoter"%in%colnames(CpGs_annotation))
    {direct_info<-direct_info[which(!CpGs_annotation[direct_info[,"regulator_name"],"promoter"]),]
    control_CpGs<-(CpGs_annotation[which((!CpGs_annotation[,"promoter"])),"name"])
    }else
    {
        control_CpGs<-(CpGs_annotation[,"name"])
    }

    chrom18names<-c("EnhG1","EnhG2","EnhA1","EnhA2",
                    "EnhWk","ZNF_repeat","Het","TssBiv","EnhBiv",
                    "ReprPC","ReprPCWk","Quies")
    types<-unique(direct_info$type)
    hmm18CpGs<-data.frame()
    for(i in 1:length(chrom18names))
    {
        myname<-chrom18names[i]
        olp<-vector()
        for(j in 1:length(types))
        {olp<-c(olp,length(which(CpGs_annotation[
            unique(direct_info[direct_info$type==types[j],"regulator_name"]),chromcol]==myname)))
        }
        total_num<-length(which(CpGs_annotation[control_CpGs,chromcol]==myname))
        hmm18CpGs<-rbind(hmm18CpGs,
                         c(olp,total_num))
    }

    rownames(hmm18CpGs)<-chrom18names
    colnames(hmm18CpGs)<-c(types,"total_control")

    total_CpGs<-colSums(hmm18CpGs)

    hmm18CpGs<-rbind(total_CpGs,hmm18CpGs)
    rownames(hmm18CpGs)[1]<-"total_num"

    return(hmm18CpGs)
}



