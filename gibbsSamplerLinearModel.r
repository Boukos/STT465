

X=as.matrix(model.matrix(~factor(GENDER)+factor(Litter)+factor(cage)))[,-1])
X=scale(X)
X=cbind(1,X)












rm(list=ls())

library(BGLR)
data(mice)
XF=scale(as.matrix(model.matrix(~factor(GENDER)+factor(Litter),data=mice.pheno))[,-1])
XR=scale(as.matrix(model.matrix(~factor(cage)-1,data=mice.pheno))[,-1])
XF=cbind(1,XF)
y=scale(mice.pheno$Obesity.BMI)
XR=XR/sqrt(ncol(XR))




fmOLS=lm(y~XF[,-1]+XR)
bHatOLS=coef(fmOLS)
bHatOLS<-ifelse(is.na(bHatOLS),0,bHatOLS)

fmBGLR=BGLR(y=y,ETA=list( 
                            list(X=XF[,-1],model='FIXED'), 
                            list(X=XR,model='BRR')),
            nIter=12000,burnIn=2000)


bHatBGLR=c(fmBGLR$mu,fmBGLR$ETA[[1]]$b,fmBGLR$ETA[[2]]$b)


## Parameters

 nIter=1200
 burnIn=200

 df0=1
 Se=var(y)/2*(df0+2)
 Sb=c(1e4, var(y)/3/ncol(XR)*(df0+2))
 
 groups<-c(rep(1,each=ncol(XF)),rep(2,ncol(XR)))


 X<-cbind(XF,XR)
 sumX2<-colSums(X^2)
 p=ncol(X)
 n=nrow(X)

 B=matrix(nrow=nIter,ncol=p,0)
 nGroups<-length(unique(groups))
 varB=matrix(nrow=nIter,ncol=nGroups,0)
 varE<-rep(0,nIter)
 varB[1,2]=Sb[2]/(df0+2)
 varB[,1]=1e4

 B[1,1]=mean(y)
 error=y-mean(y)

 beta=B[1,]

 DFe=n+df0
 DFb=p+df0

 for(i in 2:nIter){
	## Variance
  	 S=sum(error^2)+Se
	 varE[i]=S/rchisq(df=DFe,n=1)

     S=sum(beta^2)+Sb[2]
     varB[i,2]=S/rchisq(df=DFb,n=1)
 
 	  
	for(j in 1:p){
		xj=X[,j]
		error=error+xj*beta[j]
			C=sumX2[j]/varE[i]+1/varB[i,groups[j]]
			rhs<-sum(xj*error)/varE[i]
			sol=rhs/C
			beta[j]=rnorm(1,sd=sqrt(1/C))+sol
		error=error-X[,j]*beta[j]
			
	
	}
	B[i,]=beta
    print(i)
    

}





fmBayes=BGLR(y=y,ETA=list( 
                            list(X=XF,model='FIXED')),
            nIter=52000,burnIn=2000)


## Plot of estimated fixed effects
 plot(bHatOLS,bHatBayes)








install.packages('BGLR')
library(BGLR)
library(MASS)   #for mvrnorm;

data(mice)
rm(mice.A); rm(mice.X)
y<-mice.pheno$Obesity.BMI
sex<-factor(mice.pheno$GENDER)
litterSize<-factor(mice.pheno$Litter)
cage<-factor(mice.pheno$cage)

#### 3 ####
XF <- as.matrix(model.matrix(~sex+litterSize))[,-1]

XR <- as.matrix(model.matrix(~cage))[,-1]
XR <- scale(XR, scale=T, center=T)
XF <- scale(XF, scale=T, center=T)
XF=cbind(1,XF)
y <-  mice.pheno$Obesity.BMI  ; y=scale(y) ###
#parameter
nIter <- 2200
burnIn <- 200
dfE=dfBR=4
varBF <- 1e6
pXR <- ncol(XR); pXF <- ncol(XF)
sE <- var(y)/3*(dfE+2) 
sBR <- var(y)/3*(dfBR+2)/pXR



#objects store samples
BF <- matrix(nrow=nIter, ncol=pXF,0) 
BR <- matrix(nrow=nIter, ncol=pXR,0)
varBR = varE = rep(NA, nIter)
#simple computation
post.dfBR <- ncol(XR)+dfBR 
post.dfE <- nrow(XR)+dfE
XFXF <- crossprod(XF)
XRXR <- crossprod(XR)
#initialization
varE[1] <- var(y)/3
varBR[1] <- var(y)/3/pXR
BR[1,] <- rep(0, pXR) 

#Gibbs Sampler
for(i in 2:nIter){
  #smaple BF
  CBF <- XFXF/varE[(i-1)]
  diag(CBF) <- diag(CBF)+(1/varBF)
  yF <- y-(XR%*%BR[(i-1),])
  XFyF <- crossprod(XF,yF)
  rhsBF <- XFyF/varE[(i-1)]
  CBFInv <- chol2inv(chol(CBF))
  solBF <- CBFInv%*%rhsBF
  BF[i,] <- mvrnorm(mu=solBF, Sigma=CBFInv, n=1)
  
  #sample BR
   CBR <- XRXR/varE[(i-1)]
   diag(CBR) <- diag(CBR)+(1/varBR[(i-1)])
   yR <- y-(XF%*%BF[i,])
   XRyR <- crossprod(XR, yR)
   rhsBR <- XRyR/varE[(i-1)]
   CBRInv <- chol2inv(chol(CBR))
   solBR <- CBRInv%*%rhsBR
   BR[i,] <- mvrnorm(mu=solBR, Sigma=CBRInv, n=1)
  #sample varBR
   s <- sBR+sum(BR[(i),]^2)
   varBR[i] <- s/rchisq(df=post.dfBR, n=1)
  #sample varE
  
  error <- y-XF%*%BF[(i),]-XR%*%BR[(i),]
  s <- sE+sum(error^2)
  varE[i] <- s/rchisq(df=post.dfE, n=1)
  
  print(varE[i])
}


####posterior mean and CI
pBF <- BF[-c(1:burnIn),]
pBR <- BR[-c(1:burnIn),]
p.varE <- varE[-c(1:burnIn)]
p.varBR <- varBR[-c(1:burnIn)]
post <- function(n){
  n.post.mean <- mean(n)
  CI <- quantile(n, probs=c(0.025,0.975))
  list(n.post.mean, CI)
}
post(pBF[,1])
post(pBF[,2])
post(pBF[,3])
post(pBF[,4])
post(pBF[,5])
post(pBF[,6])
post(pBF[,7])
post(pBF[,8])
post(p.varE)
post(p.varBR)

#trace plot
plot(varBR, cex=0.5, type='o')
plot(varE, cex=0.5, type='o')
plot(density(varBR))
plot(density(varE))
