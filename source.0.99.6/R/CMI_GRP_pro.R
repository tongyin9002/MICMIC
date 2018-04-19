
################MICMIC#############################################
###Copyright (C) 2006-2018 Tong Yin <tongyin9002@gmail.com>########
###This program is released under the [GPL], version 3 or later.###


#' parallel PC network construction based on MI/CMI testing
#'
#' \code{PC_para} is a parallel computation method to infer direct correlation network from data matrix. This method is based on PC-algorithm by conditional mutual information
#' It will generate an adjacent matrix of the infered network.
#' @param data_matrix a numeric data matrix containing data from observation where columns contain samples(observing) and rows contain variables
#' @param max_L The max L of PC. The default value is 1, and that means the network will be infered by CMI testing. If the value is 0, the network will be infered by MI testing.
#' @param method choose a  to test interaction between nodes based on conditional mutual information (CMI), or conditional mutual inclusive information.
#' @param pre_adj the pre-defined adjacent matrix, representing the hypothetical network. The default value is NULL, and that means all nodes are considered to have association between each other in original hypothesis.
#' @param log_file_dir a string of file directory to store the log files. If the parameter is not specified, the log file directory will be get by \code{getwd()}.
#' @param edgemode a string value to select the mode in edge decision
#' @param pvalue_cut the cutoff of pvalue. The default is 0.01.
#' @param core_num the number of CPUs using in the computation.
#' @param permutation_times the number of times of permutation to calculate the pvalue
#' @usage PC_para(data_matrix,max_L=1,method=c("CMII","CMI"),pre_adj=NULL
#' ,log_file_dir=NA,edgemode=c("pvalue","MI"),pvalue_cut=0.01,core_num=1,permutation_times=100)
#' @return the adjacency matrix of the network with value of 0 and 1. 1 means that there is an edge between the rowname and colname of the element.
#' And 0 means there is no edge.
#' @examples
#' x=rnorm(300,mean=20,sd=6)
#' y=x+rnorm(300,mean=0,sd=2)
#' w=y*0.1+rnorm(300,mean=18,sd=1)
#' v=y*0.15+rnorm(300,mean=17,sd=1)
#' z=2*w+v+rnorm(300,mean=0,sd=0.1)
#' a=rnorm(300,mean=20,sd=2)
#' b=0.9*a+rnorm(300,mean=2,sd=1)
#' c=b-rnorm(300,mean=0,sd=2)
#' mydata<-rbind(x,y,w,v,z,a,b,c)
#' MI_PC_net<-PC_para(mydata,max_L=0,log_file_dir=tempdir())
#' CMI_PC_net<-PC_para(mydata,max_L=1,method="CMI"log_file_dir=tempdir())
#' @author Tong Yin
#' @references
#' Zhang, X. (2011). Inferring gene regulatory networkds from gene expression data by path consistency algorithm based on conditional mutual information
#' @references
#' Zhang, X. (2015). Conditional mutual inclusive information enables accurate quatification of associations in gene regulatory networks.
#' @references
#' Kalisch, M. and Buhlmann, P.(2007) Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm.
#' @references
#' Pethel, S.D. and Hahs, D.W. (2014). Exact Test of Independence Using Mutual Information
#' @importFrom utils write.table
#' @export



PC_para<-function(data_matrix,max_L=1,method=c("CMII","CMI"),pre_adj=NULL,log_file_dir=NA,edgemode=c("pvalue","MI"),pvalue_cut=0.01,core_num=1,permutation_times=100,MI_cut=0.1)
{
    options(stringsAsFactors = FALSE)
    method <- match.arg(method)
    edgemode <- match.arg(edgemode)

    cat("edgemode: ");cat(edgemode);cat("\n")
    if(is.na(log_file_dir))
    {
     log_file_dir<-getwd()
    }
 ###define the adjacent matrix

    if(is.null(pre_adj))
    { adj_matrix<-matrix(1,nrow=nrow(data_matrix),ncol=nrow(data_matrix))
      colnames(adj_matrix)<-rownames(data_matrix)
      rownames(adj_matrix)<-rownames(data_matrix)
    }else
    { adj_matrix<-pre_adj
      data_matrix<-data_matrix[rownames(adj_matrix),]
    }

    sample_number<-ncol(data_matrix)
    variables_number<-nrow(data_matrix)
    log_content<-list()
    ##
    cat("The number of samples: ");cat(sample_number);cat("\n")
    cat("The number of variables: ");cat(variables_number);cat("\n")
    ##
    cat("Adjacent matrix is consistent with data names: ");cat(identical(rownames(data_matrix),rownames(adj_matrix)));cat("\n")

 ###generate the original network, using Network Class
    mynet<-Network()
    mynet<-read_adj_matrix(mynet,adj_matrix)

 ###clear the memory
    rm(pre_adj);rm(adj_matrix);gc()

 ###sTART
    cat("Network construction begin, using PC Algorithm!");cat("\n")

   #start L runs of PC
    for(L in 0:max_L)
    {   gc()  #clear memory at each time
        cat("RUN L = ");cat(L);cat("\n")

        if(L==0)
        {
         result<-PC_one_run_para_L0(data_matrix,mynet@edges,L,method=method,log_file_dir=log_file_dir,pvalue_cut=pvalue_cut,core_num=core_num,permutation_times=permutation_times,edgemode=edgemode)
        }

        if(sum(mynet@adj_matrix)>0)
        {
         if(L==1)
         {
          result<-PC_one_run_para_L1(data_matrix,mynet,L,method=method,log_file_dir=log_file_dir,pvalue_cut=pvalue_cut,core_num=core_num,permutation_times=permutation_times,edgemode=edgemode)
         }
        }else
        {
          result<-data.frame()
        }

        cat("This run finished!");cat("\n")
        if(nrow(result)>0)
        {cat("deleting edges... \n")
         mynet<-delete_edge(mynet,result)    #delete the false edges detected by PC
        }
        if(!is.na(log_file_dir))
        {
         write.table(mynet@adj_matrix,paste(log_file_dir,"/adj_log",L,".txt",sep=""),sep="\t",quote=FALSE)
        }
        rm(result);gc()
    }

   #write the result
    nodenames<-rownames(mynet@adj_matrix)
    write.table(cbind(nodenames[mynet@edges[,1]],nodenames[mynet@edges[,2]]),paste(log_file_dir,"/final_edges.txt",sep=""),sep="\t",quote=FALSE)

    cat("Construction of network finished!");cat("\n")

    result<-mynet@adj_matrix
    rm(mynet);gc()

    return(result)
}









#' @importFrom utils write.table

PC_one_run_para_L0<-function(data_matrix,L0_edges,L,method=c("CMII","CMI"),log_file_dir=NA,edgemode=c("pvalue","MI"),pvalue_cut=0.01,core_num=1,permutation_times=100,MI_cut=0.1)
{   options(stringsAsFactors = FALSE)
    method <- match.arg(method)
    edgemode <- match.arg(edgemode)

    delete_edges<-data.frame()

    NodeA<-vector();NodeB<-vector();condition<-vector()

    cat(nrow(L0_edges));cat(paste(" edges having ",L," alternative routes and to be determined",sep=""));cat("\n")

    if(nrow(L0_edges)>0)
    {
        if(!is.na(log_file_dir))
        {write.table(L0_edges,paste(log_file_dir,"/testedges_log",L,".txt",sep=""),sep="\t",quote=FALSE)
        }
 ###Start parallel MI test
        time_consuming<-system.time(result_para<-callmcL0test(L0_edges,data_matrix,core_num=core_num,permutation_times=permutation_times))
        cat(time_consuming);cat("\n");cat("Parallel MI test finished! \n")

 ###build a result data frame to store the MI value and pvalue for each edge
        nodenames<-rownames(data_matrix)
        L0_result<-cbind(nodenames[L0_edges[,1]],nodenames[L0_edges[,2]],rep(NA,nrow(L0_edges)),as.numeric(result_para[,1]),as.numeric(result_para[,2]),as.numeric(result_para[,3]),as.numeric(result_para[,4]))
        colnames(L0_result)<-c("NodeA","NodeB","condition","nor_MI_value","median_of_test","adj.pvalue","MI")
        L0_result2<-cbind(nodenames[L0_edges[,2]],nodenames[L0_edges[,1]],rep(NA,nrow(L0_edges)),as.numeric(result_para[,1]),as.numeric(result_para[,2]),as.numeric(result_para[,3]),as.numeric(result_para[,4]))
        colnames(L0_result2)<-c("NodeA","NodeB","condition","nor_MI_value","median_of_test","adj.pvalue","MI")
        L0_result<-rbind(L0_result,L0_result2)

        rownames(L0_result)<-paste(L0_result[,1],L0_result[,2])

        cat("tested edges number: ");cat(nrow(L0_result)/2);cat("\n")

        if(edgemode=="pvalue")
        {
         cat("pvalue cutoff:");cat(pvalue_cut);cat("\n")
         pvalue<-as.numeric(L0_result[,"adj.pvalue"])
         remain<-((pvalue)<(pvalue_cut))
        }
        if(edgemode=="MI")
        {
         cat("MI cutoff:");cat(MI_cut);cat("\n")
         MI<-as.numeric(L0_result[,"MI"])
         remain<-((MI)>=(MI_cut))
        }

        cat("remain edges: ");cat(length(which(remain))/2)
        cat(" in ");cat(length(remain)/2);cat("\n")

        delete_edges<-L0_result[which(!remain),c("NodeA","NodeB")]

        cat("L0 finished!\n");cat(paste(nrow(delete_edges)/2," edges deleted!\n",sep=""))
    }

    write.table(L0_result,paste(log_file_dir,"/mi_log",L,".txt",sep=""),sep="\t",quote=FALSE)
    write.table(delete_edges,paste(log_file_dir,"/del_log",L,".txt",sep=""),sep="\t",quote=FALSE)

    rm(data_matrix);rm(pvalue);rm(L0_result);rm(L0_edges);gc()

    return(delete_edges)
}















#' @importFrom utils write.table

PC_one_run_para_L1<-function(data_matrix,mynet,L,method=c("CMII","CMI"),log_file_dir=NA,edgemode=c("pvalue","MI"),pvalue_cut=0.01,core_num=1,permutation_times=100,MI_cut=0.1)
{ options(stringsAsFactors = FALSE)
  method <- match.arg(method)
  edgemode <- match.arg(edgemode)

  delete_edges<-data.frame()
  totaledgenum<-get_edge_number(mynet)

  L1_edges_num<-nrow(mynet@edges)
  NodeA<-vector()
  NodeB<-vector()
  condition<-vector()

  L1edges<-mynet@edges

  sharing_partners<-mclapply(seq_len(nrow(L1edges)),function(x){get_sharing_partners(mynet,L1edges[x,1],L1edges[x,2])},mc.cores = core_num)

  sharing_if<-sapply(seq_len(nrow(L1edges)),function(x){!is.na(sharing_partners[[x]][1])})
  L1_edges_num<-sum(as.numeric(sharing_if))

  NodeA<- unlist(lapply(seq_len(nrow(L1edges)),function(x){if(sharing_if[x]){rep(L1edges[x,1],length(sharing_partners[[x]]))}}))
  NodeB<- unlist(lapply(seq_len(nrow(L1edges)),function(x){if(sharing_if[x]){rep(L1edges[x,2],length(sharing_partners[[x]]))}}))
  condition<-unlist(sharing_partners)

  L1_edgesinCondition<-matrix(0,ncol=3,nrow=length(NodeA))
  storage.mode(L1_edgesinCondition)<-"integer"
  L1_edgesinCondition[,1]<-NodeA
  L1_edgesinCondition[,2]<-NodeB
  L1_edgesinCondition[,3]<-condition

  rm(sharing_if);rm(sharing_partners);rm(mynet);rm(L1edges);rm(NodeA);rm(NodeB);rm(condition);gc()

  cat(L1_edges_num);cat(paste(" edges having ",L," alternative routes and to be determined",sep=""))
  cat(paste(" in ",nrow(L1_edgesinCondition)," conditions",sep=""));cat("\n")

  L1_result<-data.frame(NA,NA,NA,NA,NA,NA,NA);colnames(L1_result)<-c("NodeA","NodeB","condition","nor_CMI_value","median_of_test","adj.pvalue","CMI")
  L1_result<-L1_result[-1,]
  delete_edges<-data.frame(NA,NA)
  delete_edges<-delete_edges[-1,]

  if(nrow(L1_edgesinCondition)>0)
  {
    if(!is.na(log_file_dir))
    {
      write.table(L1_edgesinCondition,paste(log_file_dir,"/testedges_log",L,".txt",sep=""),sep="\t",quote=FALSE)
    }

    system.time(result_para<-callmcL1test(L1_edgesinCondition,data_matrix,core_num=core_num,permutation_times=permutation_times))

    cat("Parallel CMI test finished! \n")

    nodeNames<-rownames(data_matrix)
    L1_result<-data.frame(nodeNames[L1_edgesinCondition[,1]],nodeNames[L1_edgesinCondition[,2]],nodeNames[L1_edgesinCondition[,3]],result_para)

    colnames(L1_result)<-c("NodeA","NodeB","condition","nor_CMI_value","median_of_test","adj.pvalue","CMI")
    L1_result2<-cbind(L1_result[,2],L1_result[,1],L1_result[,c(3:7)])
    colnames(L1_result2)<-c("NodeA","NodeB","condition","nor_CMI_value","median_of_test","adj.pvalue","CMI")
    L1_result<-rbind(L1_result,L1_result2)

    rm(L1_edgesinCondition)
    rm(L1_result2)
    gc()

    rownames(L1_result)<-paste(L1_result[,1],L1_result[,2],L1_result[,3])

    #remain<-((L1_result[,"nor_CMI_value"]>L1_result[,"median_of_test"])&(L1_result[,"adj.pvalue"]<pvalue_cut))
    if(edgemode=="pvalue")
    {
     remain<-((as.numeric(L1_result[,"adj.pvalue"]))<(pvalue_cut))
    }
    if(edgemode=="MI")
    {
     remain<-((as.numeric(L1_result[,"CMI"]))>=(MI_cut))
    }

    if(method=="CMI")
    {
      delete_edges<-L1_result[which(!remain),c("NodeA","NodeB")]
    }
    if(method=="CMII")
    {
      L1_result<-data.frame(L1_result,remain)
      alterpathA<-paste(L1_result[which(!remain),1],L1_result[which(!remain),3],L1_result[which(!remain),2])
      alterpathB<-paste(L1_result[which(!remain),2],L1_result[which(!remain),3],L1_result[which(!remain),1])
      delete_result<-L1_result[which(!remain),][which(L1_result[alterpathA,"remain"]&L1_result[alterpathB,"remain"]),]
      delete_edges<-delete_result[,1:2]
      delete_edges<-delete_edges[!duplicated(paste(delete_edges[,1],delete_edges[,2])),]
    }

    cat("L1 finished!\n");cat(paste(nrow(delete_edges)/2," edges deleted!\n",sep=""))
  }

  write.table(L1_result,paste(log_file_dir,"/mi_log",L,".txt",sep=""),sep="\t",quote=FALSE)
  write.table(delete_edges,paste(log_file_dir,"/del_log",L,".txt",sep=""),sep="\t",quote=FALSE)

  rm(L1_result)
  gc()

  return(delete_edges)
}





#' @import parallel

callmcL0test<-function(L0_edges,data_matrix,core_num=1,permutation_times=30)
{
    system.time(result_para<-unlist(mclapply(seq_len(nrow(L0_edges)),function(x){MI(X=data_matrix[L0_edges[x,1],],Y=data_matrix[L0_edges[x,2],],unit="normalized",pvalue=TRUE,permutation_times=permutation_times)}, mc.cores = core_num)))
    result<-matrix(result_para,ncol=7,nrow=length(result_para)/7,byrow=TRUE)
    colnames(result)<-names(result_para[1:7])
    myMI<-(unlist(mclapply(seq_len(nrow(L0_edges)),function(x){MI(X=data_matrix[L0_edges[x,1],],Y=data_matrix[L0_edges[x,2],])}, mc.cores = core_num)))
    result<-data.frame(result,myMI)
    #rownames(result)<-rownames(L0_edges)
    return(result[,c("nor_MI_value","median_of_test","adj.pvalue","myMI")])
}





#' @import parallel

callmcL1test<-function(L1_edges,data_matrix,core_num=1,permutation_times=30)
{
    system.time(result_para<-unlist(mclapply(seq_len(nrow(L1_edges)),function(x){CMI(X=data_matrix[L1_edges[x,1],],Y=data_matrix[L1_edges[x,2],],Z=data_matrix[L1_edges[x,3],],unit="normalized",pvalue=TRUE,permutation_times=permutation_times)}, mc.cores = core_num)))
    result<-matrix(result_para,ncol=8,nrow=length(result_para)/8,byrow=TRUE)
    colnames(result)<-names(result_para[1:8])
    myCMI<-(unlist(mclapply(seq_len(nrow(L1_edges)),function(x){CMI(X=data_matrix[L1_edges[x,1],],Y=data_matrix[L1_edges[x,2],],Z=data_matrix[L1_edges[x,3],])}, mc.cores = core_num)))
    result<-data.frame(result,myCMI)
    return(result[,c("nor_CMI_value","median_of_test","adj.pvalue","myCMI")])
}
















