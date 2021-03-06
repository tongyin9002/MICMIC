% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMII_functions.R
\name{entropy}
\alias{entropy}
\title{entropy}
\usage{
entropy(X,Y=NULL,method=c("covariance","density"),
unit=c("bits","nats","hartley"),variable=c("continuous","discrete"))
}
\arguments{
\item{X}{a numeric vector to test}

\item{Y}{a numeric vector to test, default is NULL. If Y is given, then the joint entropy of X and Y will be calculated.}

\item{method}{the method to estimate the probability distribution: "covariance" or "density" method. The covariance method uses equation covariance matrix which was describled by Zhang, X in 2012. And the density method use the \code{density()} and \code{kde2d()} function to estimate the variables' density.}

\item{unit}{The unit of the result: "bits", "nats", "hartley" (the default is "bits").}

\item{variable}{variable type: "continuous" or "discrete"}
}
\value{
a numeric value of entropy
}
\description{
\code{entropy} takes discrete or continuous as input and calculate the entropy of X or joint entropy of X and Y.
}
\examples{
x1<-rnorm(100,mean=50,sd=16);x2<-c(1:100);x3<-c(1:100)+rnorm(100)
entropy(x1)
entropy(x2)
entropy(x3)
entropy(X=x1,Y=x3)
}
\references{
Zhang, X., Zhao, X. M., He, K., Lu, L., Cao, Y., Liu, J., ... & Chen, L. (2012). Inferring gene regulatory networks from gene expression data by path consistency algorithm based on conditional mutual information. Bioinformatics, 28(1), 98-104.

Moon, Y. I., Rajagopalan, B., & Lall, U. (1995). Estimation of mutual information using kernel density estimators. Physical Review E, 52(3), 2318.

Venables, W. N., & Ripley, B. D. (2013). Modern applied statistics with S-PLUS. Springer Science & Business Media.
}
\author{
Tong Yin
}
