### Gibbs sampler for a muliple linear regression with 'fixed' and random effectslinear

  The code below can be used to run a Gibbs Sampler for a linear model of the form

    y=X1b2+X2b2+....+Xqbq+e

where:
   y (nx1) is the response
   X=[X1,X2,...,Xq] is the incidence matrix for effects

Errors are assumed to be iid normal. And the bj (j=1,..,q) are assumed to be normally and independently distributed with null mean and group-specific variance. Variances are assigned scaled-inverse chi-square priors.

```R 
## Data
  ## Read data here, the code below expects:
      # y (nx1) the response
      # X (nxp) an incidence matrix of effects 
      # groups (px1) grouping of effects 
      # isRandom (nGroupsx1) are the effects of the group random?

## Parameters
  nIter<-120 # use more iterations, this is just for illustration
  df0<-1  # df0 and R0 are used to determine
  R20<-.7 # the prior scale and DF of variances.
  verbose<-TRUE

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
 varE[1]=var(y)*(1-R20)
 varB[1,]=Sb/(df0+2)
 error<-y-mean(y)
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
   if(verbose){ print(c(i, ' ' ,round(varE[i],3))) }

 }
 
```
