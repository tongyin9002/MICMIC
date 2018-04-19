################MICMIC#############################################
###Copyright (C) 2006-2018 Tong Yin <tongyin9002@gmail.com>########





######################################################################
# Create a network class
#
# This is used to store the data of simple network
######################################################################

#' Network
#' An S4 class to store the network in adjacent matrix
#' @slot vertex, a string vector to store the node names
#' @slot edges, a numeric matrix with two columns to store the edges
#' @slot adj_matrix, a adjacent matrix to store the unidirectional network
#' @slot bi_adj_matrix, a adjacent matrix to sotre the bidirectional network
#' @name Network
#' @author Tong Yin

Network<-setClass("Network",slots=

                                 c(vertexs = "vector"
                                   ,edges = "matrix"
                                   ,adj_matrix = "matrix"
                                   ,bi_adj_matrix = "matrix"
                                 )


)



#' read_adj_matrix
#'
#' This function is to Create a method to read genes annotation data
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @param adj_matrix, an adjacent matrix to load
#' @return the new object loaded with adjacent matrix
#' @import methods
#' @importFrom utils object.size
#' @keywords internal
#' @author Tong Yin
setGeneric("read_adj_matrix",function(object,adj_matrix) standardGeneric("read_adj_matrix"))
#' @aliases read_adj_matrix
setMethod("read_adj_matrix","Network",function(object,adj_matrix){


  vertexs1<-rownames(adj_matrix)
  vertexs2<-colnames(adj_matrix)

  if(!all.equal(vertexs1,vertexs2))
  {
    stop("Error! The rownames and colnames of the adjacent matrix should be the same as vertexs name.")
    return("Error! The rownames and colnames of the adjacent matrix should be the same as vertexs name.")
  }

  if(!identical(adj_matrix,t(adj_matrix)))
  {
    stop("Error! The adjacent matrix should be the same as tranpose of itself.")
    return("Error! The adjacent matrix should be the same as tranpose of itself.")
  }

  object@vertexs<-vertexs1

  storage.mode(object@adj_matrix)<-"integer"
  storage.mode(object@bi_adj_matrix)<-"integer"

  object@bi_adj_matrix<-adj_matrix

  diag(object@bi_adj_matrix)<-0

  object@adj_matrix<-object@bi_adj_matrix

  n<-length(object@vertexs)
  x<-1:(n*n)
  cols<-ceiling(x/n)
  rows<-ceiling(x%%n)
  rows[which(rows==0)]<-n
  object@adj_matrix[which(cols<rows)]<-0

  cat("read in network successful!")
  cat("\n")
  cat("The size of this network is ")
  cat((object.size(object@adj_matrix)/1000000))
  cat("Mb ~ \n")

  if(length(which(object@adj_matrix==1))==0)
  {
    return(object)
  }



  object<-update_edge(object)

  #edge_pos<-which(object@adj_matrix==1)


  #edge_vertex1<-object@vertexs[ceiling(edge_pos/n)]
  #edge_vertex2<-ceiling(edge_pos%%n)
  #edge_vertex2[which(edge_vertex2==0)]<-n
  #edge_vertex2<-object@vertexs[edge_vertex2]

  #object@edges<-data.frame(edge_vertex1,edge_vertex2)

  #rownames(object@edges)=paste(object@edges[,1],object@edges[,2],sep=" ")

  cat("create edges successful!")
  cat("\n")

  return(object)

}
)



#' get_partners
#'
#' This function is to Create a method to get partners for any vertex
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @param vertex, a vector that store the names of nodes
#' @import methods
#' @return a vector of partners
#' @keywords internal
#' @author Tong Yin

setGeneric("get_partners",function(object,vertex) standardGeneric("get_partners"))
#' @aliases get_partners
setMethod("get_partners","Network",function(object,vertex){


    if(!(vertex%in%object@vertexs))
    {
        stop("Error!Your vertex is not in this network.")
    }

    return(object@vertexs[which(object@bi_adj_matrix[vertex,]==1)])

}
)


#' get_sharing_partners
#'
#' This function is to Create a method to get sharing partners for vertex
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @param vertex1, vector that store the name of the first node
#' @param vertex2, vector that store the name of the secode node
#' @return a vector of partners
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("get_sharing_partners",function(object,vertex1,vertex2) standardGeneric("get_sharing_partners"))
#' @aliases get_sharing_partners
setMethod("get_sharing_partners","Network",function(object,vertex1,vertex2){


    if((!(vertex1%in%(1:length(object@vertexs))))||(!(vertex2%in%(1:length(object@vertexs)))))
    {
        stop("Error!Your vertex is not in this network.")
    }

    return(which((object@bi_adj_matrix[vertex1,]+object@bi_adj_matrix[vertex2,])==2))

}
)






#' delete_edge
#'
#' This function is to Create a method to delete edges
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @param vertex_pairs, numeric matrix that store the pairs of vertexs
#' @return a new object that store the network which deleted edges
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("delete_edge",function(object,vertex_pairs) standardGeneric("delete_edge"))
#' @aliases delete_edge
setMethod("delete_edge","Network",function(object,vertex_pairs){


    if(test_vertex_pairs(object,vertex_pairs)=="empty")
    {
        return(object)
    }

    if((is.null(nrow(vertex_pairs)))||(nrow(vertex_pairs)==1))
    {
        object@bi_adj_matrix[as.character(vertex_pairs[1]),as.character(vertex_pairs[2])]<-0
        object@bi_adj_matrix[as.character(vertex_pairs[2]),as.character(vertex_pairs[1])]<-0

    }else
    {
        vertex_pairs<-cbind(vertex_pairs[,1],vertex_pairs[,2])
        object@bi_adj_matrix[vertex_pairs]<-0
        object@bi_adj_matrix[vertex_pairs[,c(2,1)]]<-0
    }
    object@adj_matrix<-object@bi_adj_matrix

    n<-length(object@vertexs)
    x<-1:(n*n)
    cols<-ceiling(x/n)
    rows<-ceiling(x%%n)
    rows[which(rows==0)]<-n
    object@adj_matrix[which(cols<rows)]<-0

    object<-update_edge(object)

    return(object)
}
)



#' add_edge
#'
#' This function is to Create a method to add edges
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @param vertex_pairs, numeric matrix that store the pairs of vertexs
#' @return a new object that store the network which added edges
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("add_edge",function(object,vertex_pairs) standardGeneric("add_edge"))
#' @aliases add_edge
setMethod("add_edge","Network",function(object,vertex_pairs){


    if(test_vertex_pairs(object,vertex_pairs)=="empty")
    {
        return(object)
    }

    if(is.null(nrow(vertex_pairs)))
    {
        object@bi_adj_matrix[vertex_pairs[1],vertex_pairs[2]]<-1
        object@bi_adj_matrix[vertex_pairs[2],vertex_pairs[1]]<-1
        diag(object@bi_adj_matrix)<-0
    }else
    {
        object@bi_adj_matrix[vertex_pairs]<-1
        object@bi_adj_matrix[vertex_pairs[,c(2,1)]]<-1
        diag(object@bi_adj_matrix)<-0
    }

    object@adj_matrix<-object@bi_adj_matrix

    n<-length(object@vertexs)
    x<-1:(n*n)
    cols<-ceiling(x/n)
    rows<-ceiling(x%%n)
    rows[which(rows==0)]<-n
    object@adj_matrix[which(cols<rows)]<-0

    object<-update_edge(object)

    return(object)

}
)




#' update_edge
#'
#' This function is to Create a method to update edge data from object
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @return a new object that store the network which added edges
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("update_edge",function(object) standardGeneric("update_edge"))
#' @aliases update_edge
setMethod("update_edge","Network",function(object){

    n<-nrow(object@adj_matrix)

    edge_pos<-which(object@adj_matrix==1)

    edge_vertex1<-ceiling(edge_pos/n)
    edge_vertex2<-ceiling(edge_pos%%n)
    edge_vertex2[which(edge_vertex2==0)]<-n

    object@edges<-matrix(0,ncol=2,nrow=length(edge_vertex1))
    storage.mode(object@edges)<-"integer"
    object@edges[,1]<-edge_vertex1
    object@edges[,2]<-edge_vertex2

    return(object)

}
)



#' test_vertex_pairs
#'
#' This function is to test whether the nodes in the input edges are in the network or not
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @param vertex_pairs, numeric matrix that store the pairs of vertexs
#' @return testing information
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("test_vertex_pairs",function(object,vertex_pairs) standardGeneric("test_vertex_pairs"))
#' @aliases test_vertex_pairs
setMethod("test_vertex_pairs","Network",function(object,vertex_pairs){


    vertexs<-as.vector(as.matrix(vertex_pairs))

    if(length(vertexs)==0)
    {
        return("empty")
    }

    if(!all(vertexs%in%object@vertexs))
    {
        stop("Error!Your edge containing vertexs which is not in this network!")
    }


    if(is.vector(vertex_pairs))
    {
        if(length(vertexs)<2)
        {
            stop("Error!We need at least one pair of vertexs to add/delete the edge!")
        }
    }else
    {
        if(ncol(vertex_pairs)<2)
        {
            stop("Error!We need at least one pair of vertexs to add/delete the edge!")
        }
        if(ncol(vertex_pairs)>2)
        {
            stop("Error!We need pairs of vertexs to add/delete the edges, and please put them in two columns.")
        }
    }

    return("done")

}
)



#' get_vertex_number
#'
#' This function is to get number of vertexs
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @return number
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("get_vertex_number",function(object) standardGeneric("get_vertex_number"))
#' @aliases get_vertex_number
setMethod("get_vertex_number","Network",function(object){

    return(length(object@vertexs))
}
)






#' get_edge_number
#'
#' This function is to get number of edges
#' @docType methods
#' @rdname Network
#' @param object, Object of class Network.
#' @return number
#' @import methods
#' @keywords internal
#' @author Tong Yin
setGeneric("get_edge_number",function(object) standardGeneric("get_edge_number"))
#' @aliases get_edge_number
setMethod("get_edge_number","Network",function(object){

    return(nrow(object@edges))
}
)
























