 ### Example
 
 if(FALSE){
  library(survival)
  y0=rnorm(1000,sd=4)
  threshold=2
  d=ifelse(y0<threshold,1,0)
  y=y0
  y[d==0]=threshold
  mean(y0); sd(y0)
  mean(y); sd(y)
  fm=survreg(Surv(time=y,event=d)~1,dist="gaussian")
  summary(fm)

  X=matrix(nrow=length(y),ncol=1,1)
  groups=1
  isRandom=FALSE

 
 }

## Libraries
 library(bayesm)

## Data
  ## Read data here, the code below expects:
      # y (nx1) time to event or time to censoring #*#
      # d (nx1) 1 for event, 0 for censoring       #*#
      # X (nxp) an incidence matrix of effects 
      # groups (px1) grouping of effects (integers from 1 to q, p of them mapping effects into groups)
      # isRandom (qx1) are the effects of the group random?
      # Example if you have two groups of effects, with incidence matrices X1 and X2, the first one random the 2nd one fixed,
      # then:  X=cbind(X1,X2) ; groups=c(rep(1,ncol(X1)),rep(2,ncol(X2))); isRandom=c(TRUE,FALSE)

## Parameters
  nIter<-1200 # use more iterations, this is just for illustration
  df0<-1  # df0 and R0 are used to determine
  R20<-.7 # the prior scale and DF of variances.
  verbose<-TRUE

## Renumbering groups from 1:K
   groups<-as.integer(as.factor(groups))
   nGroups<-length(unique(groups))

## Calculating hyper-parameters
  Se<-var(y ,na.rm=T)*(1-R20)*(df0+2)
  Sb<-var(y,na.rm=T)*R20/ncol(X)*(df0+2)
  p<-as.numeric(table(groups))

## determining # of censored points
  whichCensored=which(d==0) #*#
  nCensored=length(d)  #*#

# Objects that will store samples
 B<-matrix(nrow=nIter, ncol=sum(p),0)
 varE<-rep(NA,nIter)
 varB<-matrix(nrow=nIter,ncol=nGroups)

# some useful computations
 postDFe<-nrow(X)+df0
 postDFb<-p+df0
 sumSqX=colSums(X^2)

# Initialization

 B[1,1]=mean(y,na.rm=TRUE) # initialize intercept to mean(y) and other effects to zero
 varE[1]=var(y,na.rm=TRUE)*(1-R20)
 varB[1,]=Sb/(df0+2)
 error<-y-mean(y,na.rm=TRUE)
 if(nCensored>0){ error[whichCensored]=0 }

 for(i in 1:nGroups){
    varB[,i]<-ifelse(isRandom[i],varB[,i],1e4)
 }


 beta<-B[1,]

# Gibbs Sampler
 for(i in 1:nIter){

   # sampling error variance
   S=sum(error^2)+Se
   varE[i]<-S/rchisq(df=postDFe,n=1)

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
        C=sumSqX[j]/varE[i]+1/varB[i,groups[j]]
        rhs<-sum(xj*error)/varE[i]
        sol<-rhs/C
        beta[j]<-sol+z[j]/sqrt(C)
        error<-error-xj*beta[j]
    } 
    B[i,]=beta

    if(nCensored>0){#*# sampling censored points from truncated normal
      Xb=X%*%beta
      lowerBound=y[whichCensored]-Xb[whichCensored]
      error[whichCensored]<-lapply(FUN=rtrun,X=lowerBound,mu=0,sigma=sqrt(varE[i]),b=Inf)
    }

   if(verbose){ print(c(i, ' ' ,round(varE[i],3))) }

 }
 
 
 
