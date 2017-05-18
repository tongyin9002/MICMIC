test_MICMIC <- function() {
    x<-rnorm(10000)
    y<-0.5*x+rnorm(10000,sd=0.1)
    z<-0.8*y+rnorm(10000,sd=0.1)

    checkTrue(CMI(x,y,z)>CMI(x,z,y))
    checkTrue(CMI(y,z,x)>CMI(x,z,y))
    checkTrue(CMI(x,y,z,pvalue=TRUE)$adj.pvalue<0.0001)
    checkTrue(CMI(x,z,y,pvalue=TRUE)$adj.pvalue>0.05)
}


