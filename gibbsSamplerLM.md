### Gibbs sampler for a muliple linear regression with 'fixed' and random effectslinear


```R

gibbsLM<-function(y,X,groups,isRandom,nIter, df0=1,R20=.5,verbose=T){

   	## Inputs
	# y (nx1) the response, can have NAs.
    # X (nxp) an incidence matrix of effects 
    # groups (px1) grouping of effects (integers from 1 to q, p of them mapping effects into groups)
    # isRandom (qx1) are the effects of the group random?
    # Example if you have two groups of effects, with incidence matrices X1 and X2, the first one random the 2nd one fixed,
    # then:  X=cbind(X1,X2) ; groups=c(rep(1,ncol(X1)),rep(2,ncol(X2))); isRandom=c(TRUE,FALSE)

    ## Renumbering groups from 1:K
    groups<-as.integer(as.factor(groups))
    nGroups<-length(unique(groups))

	## Calculating hyper-parameters
  	Se<-var(y ,na.rm=T)*(1-R20)*(df0+2)
  	Sb<-var(y,na.rm=T)*R20/ncol(X)*(df0+2)
  	p<-as.numeric(table(groups))

	## determining NAs 
  	whichNA=which(is.na(y))
  	nNA=length(whichNA)
  
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
 	if(nNA>0){ error[whichNA]=0 }
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
    
    	if(nNA>0){
      		error[whichNA]<-rnorm(n=nNA,mean=0,sd=sqrt(varE[i]))
    	}
   		if(verbose){ print(c(i, ' ' ,round(varE[i],3))) }
 	}

  	OUT=list(varE=varE,varB=varB,B=B)
}


```

###  Example

Fitting a regression to a training data set and evaluating prediction accuracy in a testing data set.

```R
 library(BGLR)
 data(wheat)
 X=cbind(1,wheat.X)
 y=wheat.Y[,1]
 
 tst=sample(1:nrow(X),size=150)
 
 
 yTRN=y[-tst]
 XTRN=X[-tst,]
 yTST=y[tst]
 XTST=X[tst,]
 burnIn=1000
 
 fm=samples=gibbsLM(y=yTRN,X=XTRN,groups=c(1,rep(2,ncol(X)-1)),isRandom=c(FALSE,TRUE),nIter=11000,R20=.5)
 bHat=colMeans(fm$B[-(1:burnIn),])
 yHatTST=XTST%*%bHat
 cor(yTST,yHatTST)
 
```
