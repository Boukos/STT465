## Two Gibbs Samplers for Linear Models

### (1) Mixed-effects model with fixed and random effects, blocked sampling of effects



### (2) Mixed-effects model with an arbitrary number of groups of effects

```R

## Gibbs sampler linear model 

## Data
  ## Read data here, the code below expects:
      # a vector y (response)
      # an incidence matrix of effects (X)
      # groups: a vector indicating how to group the effects associated to X
      # isRandom: a vector with dimensions equal to the number of groups stating                  # whether an effect is fixed or random

## Parameters
  nIter<-12000
  burnIn<-2000
  df0<-1  # df0 and R0 are used to determine
  R20<-.7 # the prior scale and DF of variances.
  verbose<-TRUE

## Scaling X
   X<-scale(X)
   X<-cbind(1,X)

## Renumbering groups from 1:K
   groups<-as.integer(as.factor(groups))
   nGroups<-length(unique(groups))

## Calculating hyper-parameters
  Se<-var(y)*(1-R20)*(df0+2)
  Sb<-var(y)*R20/ncol(X)*(df0+2)
  p<-as.numeric(table(groups))

# Objects that will store samples
 B<-matrix(nrow=nIter, ncol=sum(p),0)
 varE<-rep(NA,nIter)
 varB<-matrix(nrow=nIter,ncol=nGroups)
 
# some useful computations
 postDFe<-nrow(X)+df0
 postDFb<-p+df0
 sumSqX=colSums(X^2)

# Initialization

 B[1,1]=mean(y) # initialize intercept to mean(y) and other effects to zero
 varE[1]=var(y)*(1-R0)
 varB[1,]=Sb/(df0+2)
 error<-y-mean(y)
 varB[1,]<-ifelse(isRandom,varB[i,],1e4)
 beta<-B[1,]

# Gibbs Sampler
 for(i in 1:nIter){
   # sampling error variance
   S=sum(error^2)+Se  
   varE[i]<-S/rchisq(df=postDFe,n=1)

   # sampling  paramters by group
    for(j in 1:nGroups){
       whichEffects<-which(groups==j)

       # sampling variances
       if(isRandom[i]){
          bj=beta[whichEffects]
          S=sum(bj^2)+Sb
          varB[i,j]<-S/rchisq(df=postDFb,n=1)
       }
     }
       # sampling effects
        xj=X[,j]
        error<-error+xj*beta[j]
          C=sumSqX[j]/varE+1/
        error<-error-xj*beta[j]
    }  
  if(verbose){ print(i) }
 }
```
