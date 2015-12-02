### A Gibbs Sampler for linear regression with binary outcomes


```R

## Libraries
 library(bayesm)

## A wrapper for sampling truncated normal
rtrun2<-function(sigma,mu,bounds){
	rtrun(sigma=sigma,mu=mu,a=bounds[1],b=bounds[2])
}
rTruncNormal<-function(sigma,mu,a,b){
	apply(rtrun2,sigma=1,mu=0,X=cbind(a,b),MARGIN=1)
}

## Data
  ## Read data here, the code below expects:
      # y (nx1) time to event or time to censoring #*#
      # d (nx1) 1 for event, 0 for censoring       #*#
      # X (nxp) an incidence matrix of effects 
      # groups (px1) grouping of effects (integers from 1 to q, p of them mapping effects into groups)
      # isRandom (qx1) are the effects of the group random?
      # Example if you have two groups of effects, with incidence matrices X1 and X2, the first one random the 2nd one fixed,
      # then:  X=cbind(X1,X2) ; groups=c(rep(1,ncol(X1)),rep(2,ncol(X2))); isRandom=c(TRUE,FALSE)

 if(TRUE){
  library(BGLR)
  data(wheat)
  X=cbind(1,svd(scale(wheat.X))$u[,1:10])
  y0=wheat.Y[,1]
  y=ifelse(y0>0,1,0)
  groups=rep(1,ncol(X))
  isRandom=FALSE
 }


## Parameters
  nIter<-1200 # use more iterations, this is just for illustration
  df0<-1  # df0 and R0 are used to determine
  R20<-.7 # the prior scale and DF of variances.
  verbose<-TRUE

## Renumbering groups from 1:K
   groups<-as.integer(as.factor(groups))
   nGroups<-length(unique(groups))

## Calculating hyper-parameters
  Sb<-R20/ncol(X)*(df0+2)
  p<-as.numeric(table(groups))

# Objects that will store samples
 B<-matrix(nrow=nIter, ncol=sum(p),0)
 varB<-matrix(nrow=nIter,ncol=nGroups)

# some useful computations
 postDFb<-p+df0
 sumSqX=colSums(X^2)

# Initialization

 B[1,1]=qnorm(p=mean(y)) # initialize intercept to match the observed mean(y) and other effects to zero
 varB[1,]=Sb/(df0+2)
 
 isOne<-y==1
 Xb=X%*%B[1,]
 a<-ifelse(isOne,-Xb,-Inf)
 b<-ifelse(isOne,Inf,-Xb)
 error<-rTruncNormal(sigma=1,mu=0,a=a,b=b)
 
 for(i in 1:nGroups){
    varB[,i]<-ifelse(isRandom[i],varB[,i],1e4)
 }


 beta<-B[1,]

# Gibbs Sampler
 for(i in 1:nIter){

   # sampling  paramters by group
   for(j in 1:nGroups){
       # sampling variances
       if(isRandom[j]){              
          S=sum(beta[groups==j]^2)+Sb
          varB[i,j]<-S/rchisq(df=postDFb[j],n=1)
       }
    }

    # sampling effects
    z<-rnorm(ncol(X))
    for(j in 1:ncol(X)){
        xj=X[,j]
        error<-error+xj*beta[j]
        C=sumSqX[j]+1/varB[i,groups[j]]
        rhs<-sum(xj*error)
        sol<-rhs/C
        beta[j]<-sol+z[j]/sqrt(C)
        error<-error-xj*beta[j]
    } 
    B[i,]=beta

    Xb=X%*%beta
     a<-ifelse(isOne,-Xb,-Inf)
 	b<-ifelse(isOne,Inf,-Xb)
 	error<-rTruncNormal(sigma=1,mu=0,a=a,b=b)
 
   if(verbose){ cat(i,'\n') }

 }
 
 
 
 ```
 
