
#' Mutual Informaiton
#'
#' \code{MI} takes two continuous variables as input and calculate the mutual information between them in various units. Different from estimator method based on data discretization, this fucntion will use covarians transformation or density estimation to estimate the continuous probabilities distribution of x and y values.
#'
#' @param X a numeric vector to test
#' @param Y a numeric vector to test
#' @param method choose an estimator method to test the mutual information: "covariance" or "KDE" (the default is "covariance").
#' @param unit The unit of the result: "bits", "nats", "hartley" and "normalized" (the default is "bits"). The normalized result will be between 0 and 1.
#' @param pvalue a logical value to determine whether to calculate the pvalue or not
#' @param permutation_times integral value to determin the permutation times in calculating p value.
#' @usage MI(X,Y,method=c("covariance","KDE"),unit=c("bits",
#' "nats","hartley","normalized"),pvalue=FALSE,permutation_times=100)
#' @return a numeric value of mutual information between X and Y
#' @author Tong Yin
#' @export
#' @importFrom MASS kde2d
#' @examples
#' x=rnorm(100);y1=rnorm(100);y2=x+rnorm(100)
#' MI(x,y1)
#' MI(x,y2)
#' MI(x,y2,pvalue=TRUE)

MI<-function(X,Y,method=c("covariance","KDE"),unit=c("bits","nats","hartley","normalized"),pvalue=FALSE,permutation_times=100)
{

  method <- match.arg(method)
  unit <- match.arg(unit)

  if(method == "covariance")
  {
    MI_value<-MI.ZhangXJ(X,Y,unit=unit,pvalue=pvalue,permutation_times=permutation_times)
  }
  if(method == "KDE")
  {
    MI_value<-MI.KDE(X,Y,unit=unit)
  }

  return(MI_value)
}



#' Conditional Mutual Information
#'
#' \code{CMI} takes three continuous variables as input and calculate the conditional mutual information between X and Y based on the condition of Z. Different from estimator method based on data discretization, this fucntion will use covarians transformation to estimate the continuous probabilities distribution of x and y values.
#'
#' @param X a numeric vector to test
#' @param Y a numeric vector to test
#' @param Z a numeric vector as the condition
#' @param method the estimator method to test the CMI: "covariance"
#' @param unit The unit of the result: "bits", "nats", "hartley" and "normalized" (the default is "bits"). The normalized result will be between 0 and 1.
#' @param pvalue a logical value to determine whether to calculate the pvalue or not
#' @param permutation_times integral value to determin the permutation times in calculating p value.
#' @usage CMI(X,Y,Z,method=c("covariance"),unit=c("bits",
#' "nats","hartley","normalized"),pvalue=FALSE,permutation_times=100)
#' @return a numeric value of conditional mutual information between X and Y based on condition of Z
#' @author Tong Yin
#' @export
#' @examples
#' x<-rnorm(100)
#' y<-0.7*x+rnorm(100,sd=0.1)
#' z<-0.8*x+rnorm(100,sd=0.1)
#' cor(x,y);cor(x,y);cor(y,z)  #correlation test cannot identify the direct connection
#' CMI(x,y,z)  # CMI identify the direct connection between x and y is not relying on
#'             #  the condiction of z
#' CMI(y,z,x)  # CMI identify the direct connection between y and z is not relying on
#'             #  the condiction of x
#' CMI(x,z,y)  # CMI identify the connection between x and z is depending on the condition of y
#' CMI(x,y,z,pvalue=TRUE)$adj.pvalue
#' CMI(x,z,y,pvalue=TRUE)$adj.pvalue
#' @references
#' Zhang, X. (2011). Inferring gene regulatory networkds from gene expression data by path consistency algorithm based on conditional mutual information
#' @references
#' Pethel, S.D. and Hahs, D.W. (2014). Exact Test of Independence Using Mutual Information

CMI<-function(X,Y,Z,method=c("covariance"),unit=c("bits","nats","hartley","normalized"),pvalue=FALSE,permutation_times=100)
{
  method <- match.arg(method)
  unit <- match.arg(unit)

  if(method == "covariance")
  {
    CMI_value<-CMI.ZhangXJ(X,Y,Z,unit=unit,pvalue=pvalue,permutation_times=permutation_times)
  }

  return(CMI_value)
}






#' entropy
#'
#' \code{entropy} takes discrete or continuous as input and calculate the entropy of X or joint entropy of X and Y.
#'
#' @param X a numeric vector to test
#' @param Y a numeric vector to test, default is NULL. If Y is given, then the joint entropy of X and Y will be calculated.
#' @param method the method to estimate the probability distribution: "covariance" or "density" method. The covariance method uses equation covariance matrix which was describled by Zhang, X in 2012. And the density method use the \code{density()} and \code{kde2d()} function to estimate the variables' density.
#' @param unit The unit of the result: "bits", "nats", "hartley" (the default is "bits").
#' @param variable variable type: "continuous" or "discrete"
#' @usage entropy(X,Y=NULL,method=c("covariance","density"),
#' unit=c("bits","nats","hartley"),variable=c("continuous","discrete"))
#' @return a numeric value of entropy
#' @author Tong Yin
#' @export
#' @examples
#' x1<-rnorm(100,mean=50,sd=16);x2<-c(1:100);x3<-c(1:100)+rnorm(100)
#' entropy(x1)
#' entropy(x2)
#' entropy(x3)
#' entropy(X=x1,Y=x3)
#' @references
#' Zhang, X., Zhao, X. M., He, K., Lu, L., Cao, Y., Liu, J., ... & Chen, L. (2012). Inferring gene regulatory networks from gene expression data by path consistency algorithm based on conditional mutual information. Bioinformatics, 28(1), 98-104.
#' @references
#' Moon, Y. I., Rajagopalan, B., & Lall, U. (1995). Estimation of mutual information using kernel density estimators. Physical Review E, 52(3), 2318.
#' @references
#' Venables, W. N., & Ripley, B. D. (2013). Modern applied statistics with S-PLUS. Springer Science & Business Media.


entropy<-function(X,Y=NULL,method=c("covariance","density"),unit=c("bits","nats","hartley"),variable=c("continuous","discrete"))
{
    method <- match.arg(method)
    unit <- match.arg(unit)
    variable <- match.arg(variable)

     if(method == "covariance")
     {
        entropy<-entropy.ZhangXJ(X=X,Y=Y,unit=unit)
     }

     if(method == "density")
     {
       if(length(Y)==0)
       {
        oneD_den_func_A<-oneD_kernel_density(X)
        entropy<-continuous_entropy(oneD_den_func_A,min(X),max(X),n=500,unit=unit)
       }else
       {
         twoD_den_func<-twoD_kernel_density(X,Y)
         entropy<-continuous_joint_entropy(twoD_den_func,lowerLimit<-c(min(X),min(Y)),upperLimit<-c(max(X),max(Y)),method="twoD_entropy")
       }
    }

    return(entropy)

}




entropy.ZhangXJ<-function(X,Y=NULL,unit=c("bits","nats","hartley"))
{
  unit <- match.arg(unit)

  XY<-as.matrix(cbind(X,Y))

  n=ncol(XY)

  C_XY<-cov(XY)

  entropy_xy<-0.5*log((2*pi*exp(1))^n*det(C_XY))

  if(unit == "bits")
  {
    entropy_xy=entropy_xy/log(2)
  }

  if(unit == "hartley")
  {
    entropy_xy=entropy_xy/log(10)
  }

  return(entropy_xy)
}


MI.ZhangXJ<-function(X,Y,unit=c("bits","nats","hartley","normalized"),pvalue=FALSE,permutation_times=100)
{
  unit <- match.arg(unit)

  X_mat<-as.matrix(X)
  Y_mat<-as.matrix(Y)
  XY_mat<-as.matrix(cbind(X,Y))

  C_X<-cov(X_mat)
  C_Y<-cov(Y_mat)
  C_XY<-cov(XY_mat)

  MI_value<-0.5*log((det(C_X)*det(C_Y))/(det(C_XY)))

  if(unit == "nats")
  {
    result_MI_value=MI_value
  }

  if(unit == "bits")
  {
    result_MI_value=MI_value/log(2)
  }

  if(unit == "hartley")
  {
    result_MI_value=MI_value/log(10)
  }


  if((unit == "normalized")||(pvalue))
  {
   orderX<-X[order(X)]

   if(cor(X,Y)>0)
   {orderY<-Y[order(Y)]
   }else
   {
    orderY<-Y[order(-Y)]
   }

   X_mat<-as.matrix(orderX)
   Y_mat<-as.matrix(orderY)
   XY_mat<-as.matrix(cbind(orderX,orderY))
   C_X<-cov(X_mat)
   C_Y<-cov(Y_mat)
   C_XY<-cov(XY_mat)
   max_MI_value<-0.5*log((det(C_X)*det(C_Y))/(det(C_XY)))
   nor_MI_value<-MI_value/max_MI_value

   if(unit == "normalized")
   {
    result_MI_value=nor_MI_value
   }

  }


  if(pvalue)
  {

    total<-permutation_times
    mat<-cbind(seq_len(total),seq_len(total))
    MIs<-apply(mat,1,function(x){MI.ZhangXJ(sample(X),sample(Y),unit="normalized")})

    ####test by permutation
    #pass_number<-length(which(MIs>nor_MI_value))
    #p_value<-pass_number/(total/2)


    ###test by Fisher's Z transformation

    sample_size<-length(X)

    hypothesis_r<-mean(MIs)

    hypothesis_z_value<-0.5*log((1+hypothesis_r)/(1-hypothesis_r))

    sd_of_z_value<-sd(0.5*log((1+MIs)/(1-MIs)))

    observed_z_value<-0.5*log((1+nor_MI_value)/(1-nor_MI_value))

    z_score<-(observed_z_value-hypothesis_z_value)/sd_of_z_value

    adj.pvalue<-2*pnorm(-abs(z_score))

    result<-list(result_MI_value,median(MIs),nor_MI_value,observed_z_value,sd_of_z_value,z_score,adj.pvalue)

    names(result)<-c("MI","median_of_test","nor_MI_value","observed_z_value","sd_of_z_value","z_score","adj.pvalue")

    return(result)
  }else
  {
   return(result_MI_value)
  }

}



MI.KDE<-function(series_A,series_B,unit=c("bits","nats","hartley","normalized"),method=c("entropy","RSteuer"),normalization=FALSE)
{
  method = match.arg(method)
  unit = match.arg(unit)

  if(unit=="normalized")
  {
    normalization<-TRUE
  }

  sample_number<-length(series_A)
  if(is.na(match(unit,c("nats","bits","hartley","normalized"))))
  {
    return("Error: unit should be nats, bits, or hartley")
  }

  if(length(series_A)!=length(series_B))
  {
    return("Error: please provide indentical length series of two variables")
  }

  oneD_den_func_A<-oneD_kernel_density(series_A)
  oneD_den_func_B<-oneD_kernel_density(series_B)
  twoD_den_func<-twoD_kernel_density(series_A,series_B)

  entropy_a<-continuous_entropy(oneD_den_func_A,min(series_A),max(series_A),n=500)
  entropy_b<-continuous_entropy(oneD_den_func_B,min(series_B),max(series_B),n=500)

  if(method == "entropy")
  {
    entropy_ab<-continuous_joint_entropy(twoD_den_func,lowerLimit<-c(min(series_A),min(series_B)), upperLimit<-c(max(series_A),max(series_B)),method="twoD_entropy")
    MI<-entropy_a+entropy_b-entropy_ab
  }


  if(method == "RSteuer")
  {

    MI<-0
    for(i in seq_len(sample_number))
    {
      fa<-oneD_den_func_A(series_A[i])
      fb<-oneD_den_func_B(series_B[i])
      fab<-densityofxy(series_A[i],series_B[i],twoD_den_func)
      MI<-MI+log(fab/(fa*fb))
    }
    MI<-MI/sample_number
  }

  if(unit=="nats")
  {

  }

  if(unit=="bits")
  {
    MI<-MI/log(2)
  }
  if(unit=="hartley")
  {
    MI<-MI/log(10)
  }

  if(normalization)
  {
    MI<-MI/(entropy_a+entropy_b)
  }
  return(MI)
}


#' @importFrom stats density
#' @importFrom stats approxfun

oneD_kernel_density<-function(x,kernel = "gaussian")
{
  n=2^round(log(length(x))+5)
  density_result<-density(x,n=n,kernel = kernel)
  OneD_den_func<-approxfun(density_result$x,density_result$y)
  return(OneD_den_func)
}



#' @importFrom MASS kde2d

twoD_kernel_density<-function(x,y)
{
  n=2^round(log(length(x))+5)
  twoD_density<-kde2d(x,y,n=n)
  return(twoD_density)
}


densityofxy<-function(x1,y1,twoD_density)
{
    numX<-which.min(abs((x1)-(twoD_density$x)))
    numY<-which.min(abs((y1)-(twoD_density$y)))
    return(twoD_density$z[numX,numY])
}



#' @importFrom stats integrate

continuous_entropy<-function(fx,xmin,xmax,n=100,unit = c("nats","bits","hartley","normalized"))
{
  integrand <- function(x) {fx(x)*log(fx(x))}
  shannon_entropy<-(-(integrate(integrand,xmin,xmax,subdivisions=n)$value))
  if(unit[1]=="bits")
  {
    shannon_entropy<-shannon_entropy/log(2)
  }
  if(unit[1]=="hartley")
  {
    shannon_entropy<-shannon_entropy/log(10)
  }
  return(shannon_entropy)
}


#' @importFrom cubature adaptIntegrate

continuous_joint_entropy<-function(den_func_xy,lowerlimit,upperlimit,method=c("twoD_entropy","adaptintegration"),unit=c("nats","bits","hartley"))
{
  method = match.arg(method)
  unit = match.arg(unit)
  minx<-lowerlimit[1]
  miny<-lowerlimit[2]
  maxx<-upperlimit[1]
  maxy<-upperlimit[2]
  if(method=="twoD_entropy")
  {
    z<-den_func_xy$z
    z<-z*log(z)
    sumz<-sum(z)
    entropy<-(sumz/length(den_func_xy$z))*((maxx-minx)*(maxy-miny))
  }
  if(method=="adaptintegration")
  {
    integrand <- function(x)
    { x1<-x[1]
      y1<-x[2]
      numX<-which.min(abs((x1)-(den_func_xy[[1]])))
      numY<-which.min(abs((y1)-(den_func_xy[[2]])))
      fxy<-den_func_xy$z[numX,numY]
      fxy*log(fxy)
    }
    entropy<-adaptIntegrate(integrand,lowerLimit<-c(minx,miny), upperLimit<-c(maxx,maxy), tol=1e-3)$integral
  }
  entropy<-(0-entropy)
  if(unit == "units")
  {
    entropy<-entropy/log(2)
  }
  if(unit == "hartley")
  {
    entropy<-entropy/log(10)
  }
  return(entropy)
}



#' @importFrom stats cov
#' @importFrom stats cor
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats median

CMI.ZhangXJ<-function(X,Y,Z,unit=c("bits","nats","hartley","normalized"),pvalue=FALSE,permutation_times=100)
{
  unit <- match.arg(unit)

  XZ<-as.matrix(cbind(X,Z))
  YZ<-as.matrix(cbind(Y,Z))
  XYZ<-as.matrix(cbind(X,Y,Z))
  Z_m<-as.matrix(Z)
  C_XZ<-cov(XZ)
  C_YZ<-cov(YZ)
  C_Z<-cov(Z_m)
  C_XYZ<-cov(XYZ)

  value<-((det(C_XZ)*det(C_YZ))/(det(C_Z)*det(C_XYZ)))
  if(value<=0)
  {
    CMI_value <- 0
  }else
  {
    CMI_value<-0.5*log(value)
  }



  if(unit == "nats")
  {
    result_CMI_value=CMI_value
  }

  if(unit == "bits")
  {
    result_CMI_value=CMI_value/log(2)
  }

  if(unit == "hartley")
  {
    result_CMI_value=CMI_value/log(10)
  }

  if((unit == "normalized")||(pvalue))
  {
    if(cor(X,Y)>0)
    {
     orderX<-X[order(X)[rank(Y)]]
    }else
    {
     orderX<-X[order(-X)[rank(Y)]]
    }

    XZ<-as.matrix(cbind(orderX,Z))
    XYZ<-as.matrix(cbind(orderX,Y,Z))
    C_XZ<-cov(XZ)
    C_XYZ<-cov(XYZ)
    value<-((det(C_XZ)*det(C_YZ))/(det(C_Z)*det(C_XYZ)))
    if(value<=0)
    {
      max_CMI_value <- 0
    }else
    {
      max_CMI_value<-0.5*log(value)
    }

    nor_CMI_value<-CMI_value/max_CMI_value

    if(nor_CMI_value==Inf)
    {
      nor_CMI_value<-0.999999
    }
    if(nor_CMI_value>1)
    {
      nor_CMI_value<-0.999999
    }

    if(unit == "normalized")
    {
      result_CMI_value=nor_CMI_value
    }

  }



  if(pvalue)
  {
    sample_size<-length(X)

    total<-permutation_times
    mat<-cbind(seq_len(total),seq_len(total))
    CMIs<-apply(mat,1,function(x){CMI.ZhangXJ(sample(X),Y,Z,unit="normalized")})

    ####test by permutation
    #pass_number<-length(which(CMIs>CMI_value))
    #p_value<-pass_number/(total/2)

    sample_size<-length(X)

    hypothesis_r<-mean(CMIs)

    hypothesis_z_values<-0.5*log((1+CMIs)/(1-CMIs))

    hypothesis_z_value<-mean(hypothesis_z_values)

    sd_of_z_value<-sd(hypothesis_z_values)

    observed_z_value<-0.5*log((1+nor_CMI_value)/(1-nor_CMI_value))

    z_score<-(observed_z_value-hypothesis_z_value)/sd_of_z_value

    adj.pvalue<-2*pnorm(-abs(z_score))

    result<-list(result_CMI_value,nor_CMI_value,median(CMIs),observed_z_value,hypothesis_z_value,sd_of_z_value,z_score,adj.pvalue)

    names(result)<-c("CMI","nor_CMI_value","median_of_test","observed_z_value","hypothesis_z_value","sd_of_z_value","z_score","adj.pvalue")
    return(result)
  }else
  {
    return(result_CMI_value)
  }


}














