library(bayesm)

gibbsLM_RC<-function(y,d,X,groups,isRandom,R20=.5,df0=1,verbose=TRUE,nIter=150){
  #renumbering groups from 1:K
  groups <- as.integer(as.factor(groups))
  nGroups <- length(unique(groups))
  
  #calculating hyper-parameter
  Se <- var(y ,na.rm=T)*(1-R20)*(df0+2)
  Sb <- var(y,na.rm=T)*R20/ncol(X)*(df0+2)
  p <- as.numeric(table(groups))
  
  #determing how many censored point
  whichCensored=which(Y$death==0) ## 
  nCensored=length(whichCensored) ##
  
  #store sample
  B <- matrix(nrow=nIter, ncol=sum(p), 0) 
  varE <- rep(NA,nIter) 
  varB <- matrix(nrow=nIter,ncol=nGroups)
  
  #computation
  postDFe <- nrow(X)+df0 
  postDFb <- p+df0 
  sumSqX=colSums(X^2)
  
  #initialization
  B[1,1]=mean(y,na.rm=TRUE) # initialize intercept to mean(y) and other effects to zero 
  varE[1]=var(y,na.rm=TRUE)*(1-R20) 
  varB[1,]=Sb/(df0+2) 
  error <- y-mean(y,na.rm=TRUE) 
  if(nCensored>0){
    error[whichCensored]=0
  }
  for(i in 1:nGroups){ 
    varB[,i]<-ifelse(isRandom[i],varB[,i],1e4) 
  }
  beta<-B[1,]
  
  # Gibbs Sampler
  for(i in 1:nIter){
    # sampling error variance (varE)
    S=sum(error^2)+Se
    varE[i]<-S/rchisq(df=postDFe,n=1)
    
    # sampling  paramters by group
    for(j in 1:nGroups){
      # sampling variances (vaeB)
      if(isRandom[j]){              
        S=sum(beta[groups==j]^2)+Sb
        varB[i,j]<-S/rchisq(df=postDFb[j],n=1)
      }
    }
    
    # sampling effects (B)
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
      error[whichCensored]<-unlist(lapply(FUN=rtrun,X=lowerBound,mu=0,sigma=sqrt(varE[i]),b=Inf))
      #    Ystar[whichCensored] <- Xb[whichCensored]+error[whichCensored]
    }
    
    if(verbose){ cat(i,round(varE[i],3),'\n') }
  }
  OUT=list(varE=varE,varB=varB,B=B)
  return(OUT)
}

