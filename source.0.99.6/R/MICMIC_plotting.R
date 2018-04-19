################MICMIC#############################################
###Copyright (C) 2006-2018 Tong Yin <tongyin9002@gmail.com>########




#' MICMIC_plotting
#'
#' This function is to map genome coordinates to plotting coordinates, and enlarge target gene promoter and gene body
#' @param gene_name The name of target gene to be plotted
#' @param met_data_matrix a SummarizedExperiment object or a numeric matrix containing CpGs methylation data where columns contain samples and rows contain variables(probe site)
#' @param exp_data_matrix a SummarizedExperiment object or a numeric matrix containing gene expression data where columns contain samples and rows contain variables(gene site)
#' @param control_id a vector containing the ids of control/normal samples
#' @param distance Integer specifying the upstream/downstream genome range to be plotted. By default distance will cover all CpGs in analysis result
#' @param ref_gene_bed a GRanges object or a data.frame containing reference gene coorinate with five columns named "name", "chr", "start", "end" and "strand". The coordinates of genes in exp_data_matrix are required to be included in this data.frame.
#' @param ref_CpGs_bed a GRanges object or a data.frame containing reference CpGS coorinate with four columns names "name", "chr", "start" and "end". The coordinates of CpGs/probes in met_data_matrix are required to be included in this data.frame.
#' @param sample_class a data.frame containing the class information for samples
#' @param outfiledir a string of file directory to store the result files. If the parameter is not specified, the log file directory will be get by \code{getwd()}.
#' @param network plot the direct regulator network or the whole network
#' @usage MICMIC_plotting(gene_name,met_data_matrix,exp_data_matrix,control_id,
#' distance=NA,ref_gene_bed,ref_CpGs_bed,sample_class,outfiledir=NA)
#' @return a ggplot object
#' @author Tong Yin
#' @export
#' @import ggplot2
#' @import gridExtra
#' @examples
#' data("TCGA_STAD_data")
#' gene_name<-"MLH1"
#' \dontrun{
#' MICMIC_plotting(gene_name=gene_name,met_data_matrix=STAD_met_data_matrix,
#' exp_data_matrix=STAD_exp_data_matrix,control_id=STAD_control_id,
#' distance=350000,ref_gene_bed=STAD_ref_gene_bed,
#' ref_CpGs_bed=STAD_ref_CpGs_bed,sample_class=sample_class)
#' }
#' @seealso the usage of CMI_met_cis_network



MICMIC_plotting<-function(gene_name,met_data_matrix,exp_data_matrix,control_id=NA,distance=NA,ref_gene_bed,ref_CpGs_bed,sample_class,outfiledir=NA,network="all")
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


    MI_net<-read.table(paste(outfiledir,"/",gene_name,"_MI_network.txt",sep=""),sep="\t")
    CMI_net<-read.table(paste(outfiledir,"/",gene_name,"_CMI_network.txt",sep=""),sep="\t")

    gene_list<-intersect(ref_gene_bed[,"name"],rownames(MI_net))

    regulator_info<-read.table(paste(outfiledir,"/",gene_name,"_regulator_info.txt",sep=""),sep="\t")
    CpGs_type<-regulator_info[,"regulator_type"]
    CpGs_list<-regulator_info[,"regulator_name"]
    CpGs_type[which(CpGs_type=="association")]<-"indirect"


    CpGs_bed<-ref_CpGs_bed[CpGs_list,]
    gene_bed<-ref_gene_bed[gene_list,]

    range<-c(gene_bed[gene_name,"start"]-distance,gene_bed[gene_name,"start"]+distance)
    target_gene_chr<-gene_bed[gene_name,"chr"]
    target_gene_strand<-gene_bed[gene_name,"strand"]
    if(target_gene_strand=="+")
    {
        target_gene_TSS<-gene_bed[gene_name,"start"]
    }else
    {
        target_gene_TSS<-gene_bed[gene_name,"end"]
    }
    target_gene_length<-abs(gene_bed[gene_name,"end"]-gene_bed[gene_name,"start"])

    #################START Plotting

    ###plot genome track

    cat("plot genome track\n")
    genome_annotate_plot<-plot_genome_track(gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,gene_bed,outfiledir)

    ##plot CpGs
    cat("plot CpGs\n")
    CpGs_plot<-plot_CpGs_track(CpGs_list,CpGs_type,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,CpGs_bed,outfiledir)

    ##plot coord conversion
    cat("plot coord conversion\n")
    coord_conversion_plot<-plot_coord_conversion_track(gene_list,CpGs_list,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,gene_bed,CpGs_bed,outfiledir)

    ##plot genes and CpGs names
    cat("plot genes and CpGs names\n")
    names_annotate_plot<-plot_annotation_track(gene_list,CpGs_list,CpGs_type,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,gene_bed,CpGs_bed,outfiledir)

    ##plot CpGs correlation pvalue
    cat("plot CpGs correlation pvalue\n")
    CpGs_pvalue_plot<-plot_CpGs_pvalue_track(CpGs_list,regulator_info,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,outfiledir)

    ##plot MI and CMI network
    cat("plot MI and CMI network\n")
    network_plot<-plot_network_track(gene_list,CpGs_list,MI_net,CMI_net,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,gene_bed,CpGs_bed,outfiledir,CpGs_type,network=network)

    ##plot CpGs correlation scatter
    cat("plot CpGs correlation scatter\n")
    data_matrix<-rbind(exp_data_matrix[gene_list,],met_data_matrix[CpGs_list,])

    CpGs_met_plot<-plot_diff_met_exp_track_new1(CpGs_list,regulator_info,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,outfiledir,sample_class,data_matrix)

    cor_direct_plot_list<-list()
    cor_indirect_plot_list<-list()
    diff_direct_plot_list<-list()

    if(length(which(CpGs_type=="direct"))>0)
    {
     cor_direct_plot_list<-plot_direct_correlation_track(CpGs_list,CpGs_type,gene_name,data_matrix,control_id,outfiledir)
     diff_direct_plot_list<-plot_CpG_diff_met_density_track(CpGs_list,CpGs_type,gene_name,data_matrix,control_id,outfiledir,sample_class)
     diff_box_direct_plot_list<-plot_CpG_diff_met_box_track(CpGs_list,CpGs_type,gene_name,data_matrix,control_id,outfiledir,sample_class)
    }

    ## merge plot

    plot_list<-list()
    plot_list[[1]]<-genome_annotate_plot
    plot_list[[2]]<-CpGs_plot
    plot_list[[3]]<-coord_conversion_plot
    plot_list[[4]]<-names_annotate_plot
    plot_list[[5]]<-CpGs_pvalue_plot
    plot_list[[6]]<-network_plot
    plot_list[[7]]<-CpGs_met_plot

    if(length(cor_direct_plot_list)>8)
    {
      layout<-rbind(rep(6,8),rep(4,8),rep(3,8),rep(2,8),rep(1,8),rep(5,8),rep(7,8))
      sample_number<-0:7
      for(i in 1:floor(length(cor_direct_plot_list)/8))
      {
       plot_list<-c(plot_list,diff_direct_plot_list[((i-1)*8+1):(i*8)])
       sample_number<-sample_number+8
       layout<-rbind(layout,sample_number)
       plot_list<-c(plot_list,cor_direct_plot_list[((i-1)*8+1):(i*8)])
       sample_number<-sample_number+8
       layout<-rbind(layout,sample_number)
       plot_list<-c(plot_list,diff_box_direct_plot_list[((i-1)*8+1):(i*8)])
       sample_number<-sample_number+8
       layout<-rbind(layout,sample_number)
      }

      plot_result<-grid.arrange(grobs=plot_list,
                              layout_matrix=layout,
                              heights=c(2,1.3,0.8,1,2.5,3,1.5,rep(c(1,3,2),(nrow(layout)-7)/3)),
                              widths=c(rep(1,ncol(layout))))

    }else
    {
     layout<-rbind(6,4,3,2,1,5,7)

     plot_result<-grid.arrange(grobs=plot_list[1:7],
                               layout_matrix=layout,
                               heights=c(2,1.3,1,1,2.5,3,2)
                              )
    }

    outfilename<-paste(outfiledir,"/",gene_name,"_met_regulation.pdf",sep="")
    ggsave(plot_result,filename=outfilename,width=45,height=18+2*nrow(layout),dpi=450,limitsize = FALSE)

    return(plot_result)
}







genome_coord_to_plot_coord<-function(mycoords,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)
{
    if(target_gene_strand=="+")
    { promoter_region<-c(target_gene_TSS-2000,target_gene_TSS+500)
    genebody_region<-c(target_gene_TSS+500,target_gene_TSS+target_gene_length)
    }
    if(target_gene_strand=="-")
    {promoter_region<-c(target_gene_TSS-500,target_gene_TSS+2000)
    genebody_region<-c(target_gene_TSS-target_gene_length,target_gene_TSS-500)
    }

    genebody_ratio<-(target_gene_length-500)/(2500)
    TTS_plot_coord<-(target_gene_length-500)/genebody_ratio+500
    distal_ratio<-(range[2]-range[1])/(2500*8)

    plot_coords<-vector()
    for(i in 1:length(mycoords))
    {coord_direction<-(mycoords[i]-target_gene_TSS)/(abs(mycoords[i]-target_gene_TSS))
    if((mycoords[i]>=promoter_region[1])&&(mycoords[i]<promoter_region[2]))
    {plot_coords<-c(plot_coords,mycoords[i]-target_gene_TSS)
    }else
    {if((mycoords[i]>=genebody_region[1])&&(mycoords[i]<genebody_region[2]))
    {plot_coords<-c(plot_coords,((abs(mycoords[i]-target_gene_TSS)-500)/genebody_ratio+500)*coord_direction)
    }else
    {if(((mycoords[i]-target_gene_TSS)*(genebody_region[2]-target_gene_TSS))>0)
    {plot_coords<-c(plot_coords,((abs(mycoords[i]-target_gene_TSS)-target_gene_length)/distal_ratio+TTS_plot_coord)*coord_direction)
    }else
    {plot_coords<-c(plot_coords,((abs(mycoords[i]-target_gene_TSS)-2000)/distal_ratio+2000)*coord_direction)
    }
    }
    }
    }
    return(plot_coords)
}




plot_genome_track<-function(gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,gene_bed,outfiledir)
{
    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    if(target_gene_strand=="+")
    { promoter_start<-(target_gene_TSS-2000);promoter_end<-(target_gene_TSS+500);TTS<-target_gene_TSS+target_gene_length
    genome_coords<-c(promoter_start,promoter_end,target_gene_TSS,TTS)
    }else
    { promoter_start<-(target_gene_TSS-500);promoter_end<-(target_gene_TSS+2000);TTS<-target_gene_TSS-target_gene_length
    genome_coords<-c(promoter_start,promoter_end,target_gene_TSS,TTS)
    }
    coord_annotates<-c("|","||","TSS","TTS")
    coord_labels<-c("_","_","TSS","TTS")

    genome_coords<-c(target_gene_TSS-(1:floor((target_gene_TSS-range[1])/100000))*100000,genome_coords)
    coord_annotates<-c(paste("-",1:floor((target_gene_TSS-range[1])/100000),"00k",sep=""),coord_annotates)
    coord_labels<-c(paste("-",1:floor((target_gene_TSS-range[1])/100000),"00k",sep=""),coord_labels)
    genome_coords<-c(genome_coords,target_gene_TSS+(1:floor((range[2]-target_gene_TSS)/100000))*100000)
    coord_annotates<-c(coord_annotates,paste("+",1:floor((range[2]-target_gene_TSS)/100000),"00k",sep=""))
    coord_labels<-c(coord_labels,paste("+",1:floor((range[2]-target_gene_TSS)/100000),"00k",sep=""))

    plot_coords<-genome_coord_to_plot_coord(genome_coords,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    plot_data<-data.frame(genome_coords,plot_coords,coord_annotates)
    rownames(plot_data)<-coord_annotates
    my_genome_plot<-ggplot(plot_data,aes(x=plot_coords, y=0))+
        geom_line()+
        scale_x_continuous(limits=plot_range,breaks=plot_coords,labels=coord_labels)+
        scale_y_continuous(limits=c(-0.5,0))+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(axis.text.x=element_text(size=28),axis.text.y=element_blank(),legend.position="none")

    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["TTS","plot_coords"],xend=(plot_data["TSS","plot_coords"]+plot_data["TTS","plot_coords"])/2,y=-0.2,yend=-0.05,colour="black",size=12)
    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["TTS","plot_coords"],xend=(plot_data["TSS","plot_coords"]+plot_data["TTS","plot_coords"])/2,y=-0.2,yend=-0.35,colour="black",size=12)
    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["TTS","plot_coords"],xend=plot_data["TSS","plot_coords"],y=-0.2,yend=-0.2,colour="black",size=12)
    my_genome_plot<-my_genome_plot+
        annotate("rect",xmin=plot_data["TSS","plot_coords"],xmax=plot_data["TTS","plot_coords"],ymin=-0.4,ymax=0,colour="lightskyblue",fill="lightskyblue",alpha=0.8)
    my_genome_plot<-my_genome_plot+
        annotate("text",x=(plot_data["TSS","plot_coords"]+plot_data["TTS","plot_coords"])/2,y=-0.25,label=gene_name,size=24,fontface="bold")
    my_genome_plot<-my_genome_plot+
        annotate("rect",xmin=min(plot_data["||","plot_coords"],plot_data["|","plot_coords"]),xmax=max(plot_data["||","plot_coords"],plot_data["|","plot_coords"]),ymin=-0.15,ymax=0,colour="greenyellow",fill="greenyellow",alpha=0.8)


    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["TSS","plot_coords"],xend=plot_data["TSS","plot_coords"],y=-0.5,yend=0,colour="grey",size=10,alpha=0.5)
    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["TTS","plot_coords"],xend=plot_data["TTS","plot_coords"],y=-0.5,yend=0,colour="grey",size=10,alpha=0.5)
    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["|","plot_coords"],xend=plot_data["|","plot_coords"],y=-0.5,yend=0,colour="grey",size=10,alpha=0.5)
    my_genome_plot<-my_genome_plot+
        annotate("segment",x=plot_data["||","plot_coords"],xend=plot_data["||","plot_coords"],y=-0.5,yend=0,colour="grey",size=10,alpha=0.5)




    gene_bed<-gene_bed[which(gene_bed[,"name"]!=gene_name),]

    if(nrow(gene_bed)>0)
    {
    starts<-gene_bed[which(gene_bed[,"strand"]=="-"),"start"]
    gene_bed[which(gene_bed[,"strand"]=="-"),"start"]<-gene_bed[which(gene_bed[,"strand"]=="-"),"end"]
    gene_bed[which(gene_bed[,"strand"]=="-"),"end"]<-starts

    gene_coords1<-genome_coord_to_plot_coord(gene_bed[,c("start")],target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)
    gene_coords2<-genome_coord_to_plot_coord(gene_bed[,c("end")],target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    start<-gene_coords1
    end<-gene_coords2
    y<-(-0.1*((1:nrow(gene_bed))%%5))
    mydata<-data.frame(start,end,y)
    colnames(mydata)<-c("start","end","y")
    my_genome_plot<-my_genome_plot+
        geom_segment(aes(x=start,xend=end,y=y,yend=y),data=mydata,arrow=arrow(length=unit(0.2, "npc")),colour="blue",size=3)
    }

    genome_coords2<-seq(from=(round(range/100000)*100000)[1],to=(round(range/100000)*100000)[2],by=100000)
    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    my_genome_plot2<-ggplot(plot_data,aes(x=plot_coords, y=0))+
        geom_line()+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2,labels=genome_coords2)+
        scale_y_continuous(limits=c(-0.5,0))+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(axis.text.x=element_text(size=30),axis.text.y=element_blank(),legend.position="none")

    my_genome_plot2<-my_genome_plot2+
        annotate("text",x=plot_range[1]+(plot_range[2]-plot_range[1])/2*1,y=-0.25,label=paste(target_gene_chr,":",round(range[1]/1000)/1000,"M -",round(range[2]/1000)/1000,"M"),size=15,fontface="bold")

    layout<-rbind(1,2)

    plot_result<-grid.arrange(my_genome_plot,my_genome_plot2,
                              layout_matrix=layout,
                              heights=c(1,1))

    #outfilename<-paste(outfiledir,"/",gene_name,"_genomeplottest.pdf",sep="")
    #ggsave(plot_result,filename=outfilename,width=48,height=6,dpi=900)

    return(plot_result)

}




plot_CpGs_track<-function(regulator_CpGs_names,CpGs_type,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,ref_CpGs_bed,outfiledir)
{
    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    CpGs_coords<-ref_CpGs_bed[regulator_CpGs_names,"start"]

    plot_coords<-genome_coord_to_plot_coord(CpGs_coords,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    coord_annotates<-regulator_CpGs_names

    h<-rep(0.4,length(CpGs_type))
    h[which(CpGs_type!="direct")]<-0.15
    plot_data<-data.frame(CpGs_coords,plot_coords,coord_annotates,CpGs_type,h)

    rownames(plot_data)<-coord_annotates

    genome_coords2<-seq(from=(round(range/400000)*400000)[1],to=(round(range/400000)*400000)[2],by=400000)
    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    my_CpGs_track<-ggplot(plot_data,aes(x=plot_coords,y=h,colour=CpGs_type))+
        geom_point(size=16,alpha=0.8)+
        scale_colour_manual(values = c("direct"="forestgreen","indirect"="gray60","weak_direct"="gray60","association"="gray60"))+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2)+
        scale_y_continuous(limits=c(0,0.5))+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position="none")
    for(i in 1:nrow(plot_data))
    {
        my_CpGs_track<-my_CpGs_track+
            annotate("segment",x=plot_data[i,"plot_coords"],xend=plot_data[i,"plot_coords"],y=0,yend=(h[i]-0.05),colour="black",size=2,alpha=0.6)
    }

    #outfilename<-paste(outfiledir,"/",gene_name,"_CpGsplot.pdf",sep="")
    #ggsave(my_CpGs_track,filename=outfilename,width=48,height=4,dpi=900)

    return(my_CpGs_track)
}




plot_CpGs_pvalue_track<-function(CpGs_names,regulator_info,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,outfiledir)
{
    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    CpGs_coords<-regulator_info[,"regulator_start"]

    plot_coords<-genome_coord_to_plot_coord(CpGs_coords,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    coord_annotates<-regulator_info[,"regulator_name"]

    pvalue<-regulator_info[,"pvalue"]
    pvalue[which(pvalue==0)]<-min(pvalue[which(pvalue>0)])
    minuslog10p<-(-log(pvalue)/log(10))
    cors<-regulator_info[,"pearson_cor"]
    direction<-rep("positive",nrow(regulator_info))
    direction[which(regulator_info[,"pearson_cor"]<0)]<-"negative"

    plot_data<-data.frame(CpGs_coords,plot_coords,coord_annotates,minuslog10p,cors,direction)

    rownames(plot_data)<-coord_annotates

    genome_coords2<-seq(from=(round(range/100000)*100000)[1],to=(round(range/100000)*100000)[2],by=100000)
    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    my_CpGs_pvalue_track<-ggplot(plot_data)+
        geom_rect(data=plot_data,mapping=aes(xmin=plot_coords-15,xmax=plot_coords+15,ymin=0,ymax=cors,fill=cors),alpha=1)+
        scale_fill_gradient2(low="red", mid = "aliceblue", high="darkseagreen")+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2)+
        scale_y_continuous(limits=c(-1,1))+
        ggtitle("correlation between CpGs methlytion and gene expression (PCC)")+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(title=element_text(size=35),
              axis.text.x=element_blank(),axis.text.y=element_blank(),
              legend.position=c(0.05,0.5),legend.text=element_text(size=30),legend.title=element_text(size=30),
              legend.key.size = unit(1,"cm"))

    my_CpGs_pvalue_track<-my_CpGs_pvalue_track+
        annotate("segment",x=plot_range[1],xend=plot_range[2],y=0,yend=0,size=5,colour="black",alpha=0.3)

    return(my_CpGs_pvalue_track)
}







plot_coord_conversion_track<-function(gene_names,CpGs_names,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,ref_gene_bed,ref_CpGs_bed,outfiledir)
{
    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    TSS<-ref_gene_bed[gene_names,"start"]
    TSS[which(ref_gene_bed[gene_names,"strand"]=="-")]<-ref_gene_bed[gene_names,"end"][which(ref_gene_bed[gene_names,"strand"]=="-")]

    genome_coords<-c(TSS,ref_CpGs_bed[CpGs_names,"start"])
    genome_coords<-genome_coords[order(genome_coords)]

    genome_coords2<-seq(from=(round(range/100000)*100000)[1],to=(round(range/100000)*100000)[2],by=100000)

    plot_coords<-genome_coord_to_plot_coord(genome_coords,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    conversion_coords<-seq(plot_range[1],plot_range[2],by=(plot_range[2]-plot_range[1])/(length(genome_coords)+1))
    conversion_coords<-conversion_coords[2:(length(conversion_coords)-1)]

    plot_data<-data.frame(genome_coords,plot_coords,conversion_coords)

    conversion_track<-ggplot(plot_data,aes(x=plot_coords,y=0))+
        geom_point(size=0.5,alpha=0.1,color="gray")+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2)+
        scale_y_continuous(limits=c(0,0.5))+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position="none")
    for(i in 1:length(genome_coords))
    {
        conversion_track<-conversion_track+
            annotate("segment",x=plot_data[i,"plot_coords"],xend=plot_data[i,"conversion_coords"],y=0,yend=0.35,colour="black",size=2,alpha=0.8)
        conversion_track<-conversion_track+
            annotate("segment",x=plot_data[i,"conversion_coords"],xend=plot_data[i,"conversion_coords"],y=0.35,yend=0.5,colour="black",size=2,alpha=0.8)

    }

    #outfilename<-paste(outfiledir,"/",gene_name,"_CpGsconversion.pdf",sep="")
    #ggsave(conversion_track,filename=outfilename,width=48,height=4,dpi=900)

    return(conversion_track)

}




plot_annotation_track<-function(gene_names,CpGs_names,CpGs_type,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,ref_gene_bed,ref_CpGs_bed,outfiledir)
{
    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    genome_annotation<-c(gene_names,CpGs_names)
    annotation_type<-c(rep("gene",length(gene_names)),CpGs_type)

    TSS<-ref_gene_bed[gene_names,"start"]
    TSS[which(ref_gene_bed[gene_names,"strand"]=="-")]<-ref_gene_bed[gene_names,"end"][which(ref_gene_bed[gene_names,"strand"]=="-")]

    genome_coords<-c(TSS,ref_CpGs_bed[CpGs_names,"start"])
    genome_coords2<-seq(from=(round(range/100000)*100000)[1],to=(round(range/100000)*100000)[2],by=100000)

    genome_annotation<-genome_annotation[order(genome_coords)]
    annotation_type<-annotation_type[order(genome_coords)]
    genome_coords<-genome_coords[order(genome_coords)]

    conversion_coords<-seq(plot_range[1],plot_range[2],by=(plot_range[2]-plot_range[1])/(length(genome_coords)+1))
    conversion_coords<-conversion_coords[2:(length(conversion_coords)-1)]

    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    plot_data<-data.frame(genome_coords,conversion_coords,genome_annotation,annotation_type)

    annotation_track<-ggplot(plot_data,aes(x=conversion_coords,y=0))+
        geom_point(size=0.5,alpha=0.1,color="gray")+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2)+
        scale_y_continuous(limits=c(0,0.4))+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position="none")
    if(length(which(annotation_type=="direct"))>0)
    {
     annotation_track<-annotation_track+
        annotate("text",x=plot_data[which(annotation_type=="direct"),"conversion_coords"],y=0.2,colour="forestgreen",label=plot_data[which(annotation_type=="direct"),"genome_annotation"],size=8,angle=90,fontface="bold")
    }
    if(length(which(annotation_type=="indirect"))>0)
    {
     annotation_track<-annotation_track+
        annotate("text",x=plot_data[which(annotation_type=="indirect"),"conversion_coords"],y=0.2,colour="grey",label=plot_data[which(annotation_type=="indirect"),"genome_annotation"],size=8,angle=90,fontface="bold")
    }
    annotation_track<-annotation_track+
        annotate("text",x=plot_data[which(annotation_type=="gene"),"conversion_coords"],y=0.2,colour="blue",label=plot_data[which(annotation_type=="gene"),"genome_annotation"],size=10,angle=90,fontface="bold")


    #outfilename<-paste(outfiledir,"/",gene_name,"_annotation.pdf",sep="")
    #ggsave(annotation_track,filename=outfilename,width=48,height=4,dpi=900)

    return(annotation_track)


}







plot_network_track<-function(gene_names,CpGs_names,MI_net,CMI_net,gene_name,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,ref_gene_bed,ref_CpGs_bed,outfiledir,CpGs_type,network=c("direct","all"))
{
    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    genome_annotation<-c(gene_names,CpGs_names)
    annotation_type<-c(rep("gene",length(gene_names)),rep("CpG",length(gene_names)))

    TSS<-ref_gene_bed[gene_names,"start"]
    TSS[which(ref_gene_bed[gene_names,"strand"]=="-")]<-ref_gene_bed[gene_names,"end"][which(ref_gene_bed[gene_names,"strand"]=="-")]

    genome_coords<-c(TSS,ref_CpGs_bed[CpGs_names,"start"])
    genome_coords2<-seq(from=(round(range/100000)*100000)[1],to=(round(range/100000)*100000)[2],by=100000)

    genome_annotation<-genome_annotation[order(genome_coords)]
    annotation_type<-annotation_type[order(genome_coords)]
    genome_coords<-genome_coords[order(genome_coords)]

    conversion_coords<-seq(plot_range[1],plot_range[2],by=(plot_range[2]-plot_range[1])/(length(genome_coords)+1))
    conversion_coords<-conversion_coords[2:(length(conversion_coords)-1)]

    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    plot_coord_data<-data.frame(genome_coords,conversion_coords,genome_annotation,annotation_type)
    rownames(plot_coord_data)<-plot_coord_data[,"genome_annotation"]

    network_track<-ggplot(plot_coord_data,aes(x=conversion_coords,y=0))+
        geom_line(alpha=0)+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2)+
        scale_y_continuous(limits=c(0,3))+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position="none")

    network_plot_data<-data.frame()
    for(i in 1:nrow(MI_net))
    {
        partners<-colnames(MI_net)[which(MI_net[i,]==1)]
        if(length(partners)>0)
        {   x1<-rep(plot_coord_data[colnames(MI_net)[i],"conversion_coords"],length(partners))
            x2<-plot_coord_data[partners,"conversion_coords"]
            lines<-which(x1>x2)
            x<-x1[lines];x1[lines]<-x2[lines];x2[lines]<-x
            network_plot_data<-rbind(network_plot_data,data.frame(x1,x2))
        }
    }
    colnames(network_plot_data)<-c("x1","x2")

    if(network!="direct")
    {
     network_track<-network_track+
        geom_curve(aes(x=x1,xend=x2,y=0,yend=0),data=network_plot_data,curvature=-0.2,color="gray60",alpha=0.3)
    }

    network_plot_data<-data.frame()
    for(i in 1:nrow(CMI_net))
    {
        partners<-colnames(CMI_net)[which(CMI_net[i,]==1)]
        if(length(partners)>0)
        {   x1<-rep(plot_coord_data[colnames(CMI_net)[i],"conversion_coords"],length(partners))
        x2<-plot_coord_data[partners,"conversion_coords"]
        lines<-which(x1>x2)
        x<-x1[lines];x1[lines]<-x2[lines];x2[lines]<-x
        network_plot_data<-rbind(network_plot_data,data.frame(x1,x2))
        }
    }
    colnames(network_plot_data)<-c("x1","x2")

    if(network!="direct")
    {
     network_track<-network_track+
        geom_curve(aes(x=x1,xend=x2,y=0,yend=0),data=network_plot_data,curvature=-0.2,color="lightgreen",alpha=0.4)

    }

    network_plot_data<-data.frame()
    i<-which(colnames(CMI_net)==gene_name)
    partners<-colnames(CMI_net)[which(CMI_net[i,]==1)]
    partners<-c(partners,colnames(CMI_net)[which(CMI_net[,i]==1)])
    partners<-intersect(partners,CpGs_names[which(CpGs_type=="direct")])

    if(length(partners)>0)
    {x1<-rep(plot_coord_data[colnames(CMI_net)[i],"conversion_coords"],length(partners))
     x2<-plot_coord_data[partners,"conversion_coords"]
     lines<-which(x1>x2)
     x<-x1[lines];x1[lines]<-x2[lines];x2[lines]<-x
     network_plot_data<-rbind(network_plot_data,data.frame(x1,x2))
    }
    colnames(network_plot_data)<-c("x1","x2")
    network_track<-network_track+
        geom_curve(aes(x=x1,xend=x2,y=0,yend=0),data=network_plot_data,curvature=-0.2,color="forestgreen",size=2)


    #outfilename<-paste(outfiledir,"/",gene_name,"_network.pdf",sep="")
    #ggsave(network_track,filename=outfilename,width=48,height=4,dpi=900)

    return(network_track)


}





plot_direct_correlation_track<-function(CpGs_names,CpGs_type,gene_name,data_matrix,control_id,outfiledir)
{
    n<-0;plots<-list()
    for(i in 1:length(CpGs_names))
    {
     if(CpGs_type[i]=="direct")
     {n<-n+1
      CpGs_name<-CpGs_names[i]
      CpGs_met<-data_matrix[CpGs_name,]
      gene_exp<-data_matrix[gene_name,]
      mydata<-data.frame(data_matrix[CpGs_name,],data_matrix[gene_name,],colnames(data_matrix),0,0)
      colnames(mydata)<-c("CpGs_met","gene_exp","sample_id")

      plot_result<-ggplot(mydata, aes(x=CpGs_met,y=gene_exp))+
         geom_point(colour="black",alpha=1,size=9)+
         geom_point(aes(colour=CpGs_met),size=7,alpha=1)+
         scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","0.2","0.4","0.6","0.8","1"))+
         scale_colour_gradient2(low ="aliceblue",mid="aliceblue",high ="gold")+
         ylab(paste(gene_name))+xlab(paste(CpGs_name))+
         labs(title=paste("PCC =",round(cor(mydata[,"CpGs_met"],mydata[,"gene_exp"]),digits=2)))+
         theme_bw()+
         theme(axis.text=element_text(size=24,face="bold"),axis.title=element_text(size=24,face="bold"),title=element_text(size=24,face="bold"),legend.position="none",
               panel.border=element_rect(linetype="solid",colour="darkslategrey",size=3),
               panel.grid.major=element_line(colour="darkslategrey"))

      mydata<-mydata[which(mydata[,"sample_id"]%in%control_id),]
      plot_result<-plot_result+geom_point(data=mydata,alpha=1,color="springgreen2",size=3,shape=18)

      plots[[n]]<-plot_result
     }
    }

    ml <- marrangeGrob(plots, nrow=1,ncol=n)

    #outfilename<-paste(outfiledir,"/",gene_name,"_direct_cor.pdf",sep="")
    #ggsave(ml,filename=outfilename,width=48,height=6,dpi=900)

    return(plots)

}





plot_indirect_correlation_track<-function(CpGs_names,CpGs_type,gene_name,data_matrix,control_id,outfiledir)
{
    n<-0;plots<-list()
    for(i in 1:length(CpGs_names))
    {
        if(CpGs_type[i]=="indirect")
        {n<-n+1
         CpGs_name<-CpGs_names[i]
         CpGs_met<-data_matrix[CpGs_name,]
         gene_exp<-data_matrix[gene_name,]
         mydata<-data.frame(CpGs_met,gene_exp,colnames(data_matrix),0,0)
         colnames(mydata)<-c("CpGs_met","gene_exp","sample_id")
         plot_result<-ggplot(mydata, aes(x=CpGs_met,y=gene_exp))+
            geom_point(alpha=1,color="black",size=5)+
            geom_point(aes(colour=CpGs_met),size=4)+
            scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","0.2","0.4","0.6","0.8","1"))+
            scale_colour_gradient2(low ="white",mid="white",high ="yellow")+
            ylab(paste(gene_name))+xlab(paste(CpGs_name))+
            labs(title=paste("PCC =",round(cor(mydata[,"CpGs_met"],mydata[,"gene_exp"]),digits=2)))+
            theme_bw()+
            theme(axis.text=element_text(size=24,face="bold"),axis.title=element_text(size=24,face="bold"),title=element_text(size=24,face="bold"),legend.position="none",
                  panel.background=element_rect(fill="snow2"),
                  panel.border=element_rect(linetype="solid",colour="black"))
         mydata<-mydata[which(mydata[,"sample_id"]%in%control_id),]
         plot_result<-plot_result+geom_point(data=mydata,alpha=1,color="springgreen2",size=3,shape=18)

         plots[[n]]<-plot_result
        }
    }

    ml <- marrangeGrob(plots, nrow=1,ncol=n)

    #outfilename<-paste(outfiledir,"/",gene_name,"_indirect_cor.pdf",sep="")
    #ggsave(ml,filename=outfilename,width=48,height=6,dpi=900)

    return(plots)

}






plot_diff_met_exp_track_new1<-function(CpGs_names,regulator_info,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range,outfiledir,sample_class,data_matrix)
{
    CpGs_number<-nrow(regulator_info)

    plot_range<-genome_coord_to_plot_coord(c(range[1],range[2]),target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    CpGs_coords<-regulator_info[,"regulator_start"]

    plot_coords<-genome_coord_to_plot_coord(CpGs_coords,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    genome_coords2<-seq(from=(round(range/100000)*100000)[1],to=(round(range/100000)*100000)[2],by=100000)
    plot_coords2<-genome_coord_to_plot_coord(genome_coords2,target_gene_chr,target_gene_TSS,target_gene_length,target_gene_strand,range)

    coord_annotates<-regulator_info[,"regulator_name"]

    class_names<-unique(sample_class[colnames(data_matrix),2])

    plot_data<-data.frame()
    for(classi in 1:length(class_names))
    {class_name<-class_names[classi]
     class_indivs<-intersect(colnames(data_matrix),sample_class[which(sample_class[,2]==class_name),1])
     sds<-apply(data_matrix[regulator_info[,"regulator_name"],class_indivs],1,sd)
     met_mean<-apply(data_matrix[regulator_info[,"regulator_name"],class_indivs],1,mean)
     myplot_data<-data.frame(CpGs_coords,plot_coords,coord_annotates,met_mean,sds,(met_mean+sds),(met_mean-sds),class_name)
     plot_data<-rbind(plot_data,myplot_data)
    }

    colnames(plot_data)<-c("CpGs_coords","plot_coords","coord_annotates","met_mean","met_sd","met_up","met_low","class")

    colors<-c("chartreuse3","deepskyblue4","gold","lightcoral","darkorchid3","firebrick4","gray8")
    colors<-c("chartreuse","blue","goldenrod1","darkorchid1","firebrick1","black")

    my_CpGs_met_track<-ggplot(data=plot_data,aes(x=plot_coords,y=met_mean,colour=factor(class),group=factor(class)))+
        geom_point(size=5,alpha=0.8)+
        geom_line(size=1.2,alpha=0.95)+
        scale_colour_manual(values=colors[1:length(class_names)])+
        scale_x_continuous(limits=plot_range,breaks=plot_coords2)+
        scale_y_continuous(limits=c(0,1))+
        ggtitle("differential methylation in samples")+
        ylab(NULL)+xlab(NULL)+
        theme_bw()+
        theme(title=element_text(size=35),
              axis.text.x=element_blank(),axis.text.y=element_blank(),
              legend.position=c(0.05,0.5),legend.text=element_text(size=30),legend.title=element_text(size=30),
              legend.key.size = unit(1,"cm"))

    my_CpGs_met_track<-my_CpGs_met_track+
        annotate("segment",x=plot_range[1],xend=plot_range[2],y=0,yend=0,size=5,colour="black",alpha=0.3)

    return(my_CpGs_met_track)
}





plot_CpG_diff_met_density_track<-function(CpGs_names,CpGs_type,gene_name,data_matrix,control_id,outfiledir,sample_class)
{
    n<-0
    plots<-list()

    for(i in 1:length(CpGs_names))
    {
        if(CpGs_type[i]=="direct")
        {n<-n+1
         CpGs_name<-CpGs_names[i]
         CpGs_met<-data_matrix[CpGs_name,]
         gene_exp<-data_matrix[gene_name,]
         mydata<-data.frame(CpGs_met,gene_exp,colnames(data_matrix),sample_class[colnames(data_matrix),2])
         colnames(mydata)<-c("CpGs_met","gene_exp","sample_id","sample_class")
         mydata<-mydata[which(!is.na(mydata[,"sample_class"])),]

         class_names<-unique(mydata[,"sample_class"])


         colors<-c("chartreuse3","deepskyblue4","gold","lightcoral","darkorchid3","firebrick4","gray8")
         colors<-c("chartreuse","blue","goldenrod1","darkorchid1","firebrick1","black")

         plot_result<-ggplot(mydata, aes(x=CpGs_met,group=sample_class,colour=sample_class,fill=sample_class))+
            geom_density(alpha=0.1)+
            scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","0.2","0.4","0.6","0.8","1"))+
            scale_colour_manual(values=colors[1:length(class_names)])+
            scale_fill_manual(values=colors[1:length(class_names)])+
            ylab("density")+xlab("beta value")+
            theme_bw()+
             theme(axis.text.x=element_text(size=24,face="bold"),axis.text.y=element_blank(),
                   axis.title.x=element_text(size=24,face="bold"),axis.title.y=element_text(size=36,face="bold"),
                   title=element_text(size=24,face="bold"),legend.position="none",
                   panel.background=element_rect(fill="gray"),
                   panel.border=element_rect(linetype="solid",colour="darkslategrey",size=3),
                   panel.grid.major=element_line(colour="darkslategrey"))

         plots[[n]]<-plot_result
        }
    }

    ml <- marrangeGrob(plots, nrow=1,ncol=n)

    #outfilename<-paste(outfiledir,"/",gene_name,"_direct_diff_density.pdf",sep="")
    #ggsave(ml,filename=outfilename,width=48,height=5,dpi=900)

    return(plots)

}








plot_CpG_diff_met_box_track<-function(CpGs_names,CpGs_type,gene_name,data_matrix,control_id,outfiledir,sample_class)
{
    n<-0
    plots<-list()

    for(i in 1:length(CpGs_names))
    {
        if(CpGs_type[i]=="direct")
        {n<-n+1
        CpGs_name<-CpGs_names[i]
        CpGs_met<-data_matrix[CpGs_name,]
        gene_exp<-data_matrix[gene_name,]
        mydata<-data.frame(CpGs_met,gene_exp,colnames(data_matrix),sample_class[colnames(data_matrix),2])
        colnames(mydata)<-c("CpGs_met","gene_exp","sample_id","sample_class")
        mydata<-mydata[which(!is.na(mydata[,"sample_class"])),]

        class_names<-unique(mydata[,"sample_class"])


        colors<-c("chartreuse3","deepskyblue4","gold","lightcoral","darkorchid3","firebrick4","gray8")
        colors<-c("chartreuse","blue","goldenrod1","darkorchid1","firebrick1","black")

        plot_result<-ggplot(mydata, aes(factor(sample_class),CpGs_met,fill=sample_class))+
            geom_boxplot()+coord_flip()+
            scale_y_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","0.2","0.4","0.6","0.8","1"))+
            scale_fill_manual(values=colors[1:length(class_names)])+
            ylab("beta value")+xlab("class")+
            theme_bw()+
            theme(axis.text.x=element_text(size=24,face="bold"),axis.text.y=element_blank(),
                  axis.title.x=element_text(size=24,face="bold"),axis.title.y=element_text(size=36,face="bold"),
                  title=element_text(size=24,face="bold"),legend.position="none",
                  panel.background=element_rect(fill="gray"),
                  panel.border=element_rect(linetype="solid",colour="darkslategrey",size=3),
                  panel.grid.major=element_line(colour="darkslategrey"))

        plots[[n]]<-plot_result
        }
    }

    ml <- marrangeGrob(plots, nrow=1,ncol=n)

    #outfilename<-paste(outfiledir,"/",gene_name,"_direct_diff_box.pdf",sep="")
    #ggsave(ml,filename=outfilename,width=48,height=5,dpi=900)

    return(plots)

}

