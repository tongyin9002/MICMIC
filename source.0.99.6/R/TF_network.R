################MICMIC#############################################
###Copyright (C) 2006-2018 Tong Yin <tongyin9002@gmail.com>########
###This program is released under the [GPL], version 3 or later.###



#' Calculate the TF binding enrichment on CpG sites
#'
#' @param motif_data Two-column data.frame containing the TF motif name in the first,
#' and CpG IDs surrounding the TF binding sites in the second column. CpG IDs should be
#' separated by ";".
#' @param CpGs_list A vector of CpG IDs to analyze.
#' @param control_list A vector of randomly selected CpG IDs as background control.
#' @usage TF_binding_enrichment(motif_data,CpGs_list,control_list)
#' @return The odds ratio and confidence interval of each TF on the set the CpG sites.
#' @export

TF_binding_enrichment<-function(motif_data,CpGs_list,control_list)
{   options(stringsAsFactors = FALSE)
    myresult<-list()
    for(i in 1:nrow(motif_data))
    {motif_name<-motif_data[i,1]
     motif_CpGs_list<-unlist(strsplit(motif_data[i,2],";"))

     myCpGs<-intersect(CpGs_list,motif_CpGs_list)
     EE<-length(myCpGs)
     EN<-length(CpGs_list)-EE
     CE<-length(intersect(control_list,motif_CpGs_list))
     CN<-length(control_list)-CE
     OR<-(EE/EN)/(CE/CN)
     ci95l<-exp(log(OR) - 1.96*(sqrt(1/EE+1/EN+1/CE+1/CN)))
     ci95h<-exp(log(OR) + 1.96*(sqrt(1/EE+1/EN+1/CE+1/CN)))
     myresult[[i]]<-c(motif_name,motif_name,EE,OR,ci95l,ci95h,paste(myCpGs,collapse=";"))
    }
    myresult<-data.frame(do.call(rbind,myresult))

    colnames(myresult)<-c("motif_name","TF","CpGs_number","OR","CI95_l","CI95_h","CpGs")
    myresult[,"CpGs_number"]<-as.numeric(myresult[,"CpGs_number"])
    myresult[,"OR"]<-as.numeric(myresult[,"OR"])
    myresult[,"CI95_l"]<-as.numeric(myresult[,"CI95_l"])
    myresult[,"CI95_h"]<-as.numeric(myresult[,"CI95_h"])
    myresult<-myresult[order(-as.numeric(myresult[,"OR"])),]

    return(myresult)
}

#' Generate the TF-target network
#'
#' @param TF_binding TF binding results from \code{TF_binding_enrichment} function
#' @param CpG_target Two-column data.frame containing CpG-target regulation pairs.
#' @usage TF_target_network(TF_binding,CpG_target)
#' @return A two-column data.frame containing the TF-target pairs.
#' @export

TF_target_network<-function(TF_binding,CpG_target)
{options(stringsAsFactors = FALSE)
 TF_network<-list()
 for(i in 1:nrow(TF_binding))
 {
    CpGs<-unlist(strsplit(TF_binding[i,"CpGs"],";"))
    TF<-TF_binding[i,"motif_name"]

    lines<-which(CpG_target[,1]%in%CpGs)
    target<-CpG_target[lines,2]

    if(length(target)>0)
    {
        TF_network[[i]]<-data.frame(TF,target)
    }
 }
 TF_network<-do.call(rbind,TF_network)
 return(TF_network)
}

#' Identify the self-connected TF circuit
#'
#' @param TF_binding TF binding results from \code{TF_binding_enrichment} function
#' @param CpG_target Two-column data.frame containing CpG-target regulation pairs.
#' @param CpGs_in_SE Two-column data.frame annotating CpGs in super enhancer
#' @usage selfconnected_TF_circuit(TF_binding,CpG_target,CpGs_in_SE)
#' @return A data.frame containing the self-regulation relationship within TFs.
#' @export

selfconnected_TF_circuit<-function(TF_binding,CpG_target,CpGs_in_SE)
{options(stringsAsFactors = FALSE)
 core_circuit_result<-data.frame()
 for(i in 1:nrow(TF_binding))
 {
    CpGs<-unlist(strsplit(TF_binding[i,"CpGs"],";"))
    motif_TF<-TF_binding[i,"motif_name"]
    targets_info<-CpG_target[which(CpG_target[,2]==TF_binding[i,"motif_name"]),]
    circuit_CpGs_lines<-which(targets_info[,1]%in%CpGs)
    if(length(circuit_CpGs_lines)>0)
    {
        core_circuit_result<-rbind(core_circuit_result,data.frame(targets_info[circuit_CpGs_lines,],motif_TF,TF_binding[i,c("OR","CI95_l","CI95_h")]))
    }
 }

 circuit_TFs<-unique(core_circuit_result[,"motif_TF"])

 for(i in 1:nrow(TF_binding))
 {
    CpGs<-unlist(strsplit(TF_binding[i,"CpGs"],";"))
    motif_TF<-TF_binding[i,"motif_name"]

    if(motif_TF%in%circuit_TFs)
    {
        targets_info<-CpG_target[which(CpG_target[,2]%in%circuit_TFs),]
        circuit_CpGs_lines<-which((targets_info[,1]%in%CpGs)&(targets_info[,2]!=motif_TF))

        if(length(circuit_CpGs_lines)>0)
        {
            core_circuit_result<-rbind(core_circuit_result,data.frame(targets_info[circuit_CpGs_lines,],motif_TF,TF_binding[i,c("OR","CI95_l","CI95_h")]))
        }
    }
 }

 core_circuit_result<-data.frame(core_circuit_result,CpGs_in_SE[core_circuit_result[,1],"super_enhancer"])

 colnames(core_circuit_result)<-c("CpGs","target","motif_TF",
                                 "OR","CI95_l","CI95_h","super_enhancer")

 return(core_circuit_result)
}


#' Filtering the TF circuit by odds ratio cutoff and super enhancer annotation
#'
#' @param TF_circuit TF circuit results from \code{selfconnected_TF_circuit} function
#' @param OR Odds ratio cutoff
#' @param super_enhancer A logical value to decide whether the TF binding CpGs in the circuit
#' should be all located in super enhancer according to previous studies.
#' @usage select_TF_circuit(TF_circuit,OR=1.05,super_enhancer=TRUE)
#' @return A data.frame containing the filtered self-regulation relationship within TFs.
#' @export

select_TF_circuit<-function(TF_circuit,OR=1.05,super_enhancer=TRUE)
{options(stringsAsFactors = FALSE)
    TF_circuit_core<-TF_circuit[which((TF_circuit[,"OR"]>OR)&(TF_circuit[,"CI95_l"]>1)),]
    if(super_enhancer)
    {
     TF_circuit_core<-TF_circuit_core[which(TF_circuit_core[,"super_enhancer"]),]
    }

    circuit_TFs<-unique(TF_circuit_core[which(TF_circuit_core[,"target"]==TF_circuit_core[,"motif_TF"]),"target"])
    TF_circuit_core<-TF_circuit_core[which((TF_circuit_core[,"target"]%in%circuit_TFs)&(TF_circuit_core[,"motif_TF"]%in%circuit_TFs)),]

    return(TF_circuit_core)
}



#' Functional analysis of TF circuits targets
#'
#' @param TF_network TF-target result from \code{TF_target_network} function
#' @param core_TF_clusters Two-column data.frame define the cluster of each TF
#' @param function_gmt A list object containing the vectors of gene list in different functional gene sets
#' @usage functional_enrichment_of_core_circuit(TF_network,core_TF_clusters,function_gmt)
#' @return A data.frame containing function name, overlapping number, enrichment p-value between gene sets and TF targets.
#' @export

functional_enrichment_of_core_circuit<-function(TF_network,core_TF_clusters,function_gmt)
{options(stringsAsFactors = FALSE)
    myresult<-data.frame()
    for(i in unique(core_TF_clusters[,2]))
    { core_TFs<-core_TF_clusters[core_TF_clusters[,2]==i,1]
      target_genes<-unique(TF_network[which((TF_network[,1]%in%core_TFs)),"target"])
      result<-cancerpathenrich(target_genes,function_gmt)
      result<-data.frame(i,rownames(result),result)
      myresult<-rbind(myresult,result)
    }

    myresult<-myresult[order(myresult[,"adj.pv"]),]
    colnames(myresult)<-c("cluster","function","olp","regSize","pv","adj.pv","olpgene")
    return(myresult)
}


cancerpathenrich<-function(genes,gmt)
{options(stringsAsFactors = FALSE)
    res <- lapply(gmt,
                  function(x,y){
                      pv <- phyper(length(intersect(x,y))-1, length(x), 20000-length(x), length(y), lower.tail = FALSE)
                      return(c(length(intersect(x,y)),length(x),pv,paste(intersect(x,y),collapse=",")))
                  },genes)

    res <- t(as.data.frame(res))

    res <- data.frame(olp=as.numeric(res[,1]),
                      regSize=as.numeric(res[,2]),
                      pv=as.numeric(res[,3]),
                      olpgene=res[,4],stringsAsFactors=FALSE)

    adj.pv <- p.adjust(res[,3],"BH")
    res$adj.pv <- adj.pv
    res <- res[,c(1:3,5,4)]

    res[[3]] <- signif(res[[3]],4)
    res[[4]] <- signif(res[[4]],4)

    res <- res[order(res[[3]]),]

    res<-res[which(res[,"pv"]<0.05),]

    rownames(res)<-sub("hsa\\d+_","",rownames(res))

    return(res)
}



#' Plotting the self-regulation network of the TF circuit
#'
#' @param TF_circuit_core The data.frame of self-regulation relationship from \code{select_TF_circuit} function
#' @param TF_class_data Two-column data.frame define the cluster of each TF
#' @param TF_pathway_data The result from functional analysis \code{functional_enrichment_of_core_circuit}
#' @param outpath The output directory for plotting file.
#' @usage plot_core_circuit(TF_circuit_core,TF_class_data,TF_pathway_data,outpath)
#' @return A network figure representing the self-regulation relationship between TFs and their functions.
#' @import ggplot2
#' @import gridExtra
#' @export

plot_core_circuit<-function(TF_circuit_core,TF_class_data,TF_pathway_data,outpath)
{options(stringsAsFactors = FALSE)
    rownames(TF_class_data)<-TF_class_data[,1]
    TF_class_data<-TF_class_data[order(TF_class_data[,2]),]

    edges<-TF_circuit_core[,c("motif_TF","target")]
    edges<-edges[which(!duplicated(paste(edges[,1],edges[,2]))),]

    core_TFs<-TF_class_data[,1]
    ###
    myplot <- ggplot()+
        theme_bw()+
        theme(axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )
    xmin<-0.8;xmax_m<-1.2;xmin_m<-1.8;xmax<-2.2
    colors<-rainbow(length(unique(TF_class_data[,2]))*1.5)
    for(i in 1:length(core_TFs))
    {myclass<-TF_class_data[core_TFs[i],2]
    myplot<-myplot+
        annotate("rect",xmin=xmin,xmax=xmax_m,ymin=i-0.2,ymax=i+0.2,alpha=1,color="black",fill=colors[myclass])
    myplot<-myplot+
        annotate("rect",xmin=xmin_m,xmax=xmax,ymin=i-0.2,ymax=i+0.2,alpha=1,color="black",fill=colors[myclass])
    myplot<-myplot+
        annotate("text",x=1,y=i,label=paste0("bold(",core_TFs[i],")"),size=6,parse = TRUE)
    myplot<-myplot+
        annotate("text",x=2,y=i,label=paste0("bold(",core_TFs[i],")"),size=6,parse = TRUE)
    myplot<-myplot+
        annotate("segment",x=xmax_m+0.05,y=i,xend=xmin_m-0.05,yend=i,colour="goldenrod",size=2)
    }

    miny<-1-0.5
    maxy<-length(core_TFs)+0.5

    TF_pathway_data<-TF_pathway_data[order(TF_pathway_data[,"cluster"]),]
    TF_pathway_data[,"function"]<-gsub("_"," ",TF_pathway_data[,"function"])
    TF_pathway_data[,"function"]<-gsub("\\.","",TF_pathway_data[,"function"])
    pathways<-unique(TF_pathway_data[,"function"])
    for(i in 1:length(pathways))
    {
        y<-(length(core_TFs)/length(pathways))*(i-1)+1
        myplot<-myplot+
            annotate("rect",xmin=xmax+2,xmax=xmax+2.6,ymin=y-0.3,ymax=y+0.3,alpha=1,color="black",fill="orange")
        lines<-which(TF_pathway_data[,"function"]==pathways[i])
        for(j in lines)
        {
            class=TF_pathway_data[j,"cluster"]
            classes_y<-which(TF_class_data[,2]==class)
            for(k in classes_y)
            {
                myplot<-myplot+
                    annotate("segment",x=xmax+0.05,y=k,xend=xmax+2,yend=y,colour=colors[class],size=2,alpha=0.5)
            }
        }

        myplot<-myplot+
            annotate("text",x=xmax+2.1,y=y,label=pathways[i],size=6)

    }

    miny<-1-0.5
    maxy<-length(core_TFs)+0.5

    myplot<-plot_motif_binding(myplot,core_TFs,edges,xmin,xmax,miny,maxy)

    myplot<-myplot+scale_x_continuous(limits = c(-0.5, 6))

    ggsave(myplot,file=paste0(outpath,"/TF_core_circuit.pdf"),width=10,height=12)
}




plot_motif_binding<-function(myplot,circuit_TFs,edges,xmin,xmax,miny,maxy)
{options(stringsAsFactors = FALSE)
    targets<-unique(edges[,"target"])
    n<-length(targets)
    cols<-rainbow(n*2)[(n/2+1):(n*2)]
    # cols<-heat.colors(n)
    names(cols)<-targets

    for(i in 1:nrow(edges))
    {
        motif_TF<-edges[i,"motif_TF"]
        target<-edges[i,"target"]
        color<-cols[target]

        motif_y<-which(circuit_TFs==motif_TF)
        target_y<-which(circuit_TFs==target)

        bottom_dis<-(abs(motif_y-miny)+abs(target_y-miny))
        top_dis<-(abs(motif_y-maxy)+abs(target_y-maxy))

        if(bottom_dis<top_dis)
        {   margin<-bottom_dis/20
        #margin<-margin*(1+abs(motif_y-target_y)/maxy)
        y<-miny-margin
        target_y<-target_y
        motif_y<-motif_y+target_y/50
        }else
        {   margin<-top_dis/20
        #margin<-margin*(1+abs(motif_y-target_y)/maxy)
        y<-maxy+margin
        target_y<-target_y
        motif_y<-motif_y+target_y/50
        }

        myplot<-myplot+
            annotate("segment",x=xmin-margin-target_y/50,y=y,xend=xmax+margin,yend=y,colour=color,size=1,alpha=0.5)
        myplot<-myplot+
            annotate("segment",x=xmin-margin-target_y/50,y=motif_y,xend=xmin-margin-target_y/50,yend=y,colour=color,size=1,alpha=0.5)
        myplot<-myplot+
            annotate("segment",x=xmax+margin,y=target_y,xend=xmax+margin,yend=y,colour=color,size=1,alpha=0.5)
        myplot<-myplot+
            annotate("segment",x=xmin-0.03,y=motif_y,xend=xmin-margin-target_y/50,yend=motif_y,colour=color,size=1,alpha=0.5)
        myplot<-myplot+
            annotate("segment",x=xmax+0.03,y=target_y,xend=xmax+margin,yend=target_y,colour=color,size=1,alpha=0.5)
    }

    return(myplot)
}




