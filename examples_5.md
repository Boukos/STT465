### Examples Chapters 5

#### (1) Inference on the mean, conditional on the variance

The following example evaluates the likelihood, prior and posterior density of the mean in a normal model
with known variance.

```R
 mu0=27
 v0=10
 x=scan('~/GitHub/STT465/bmi.txt')
 varX=var(x)
 
 ### Frequentist analysis
 ## Let's assume x is the population and let's sample from it
 n<-5
 y=sample(x,size=n,replace=F)
 meanY=mean(y)

 # a grid of values for evaluating the likelihood function
 myGrid=seq(from=10,to=50,by=.1)

 logLik<-function(x,mu,var){
    tmp=sum(dnorm(x=x,mean=mu,sd=sqrt(var),log=TRUE))
    return(tmp)
 }
 
 likelihood<-rep(NA,length(myGrid))
 for(i in 1:length(myGrid)){
    likelihood[i]<-exp(logLik(x=y,mu=myGrid[i],var=varX))
 }
 plot(x=myGrid,y=likelihood,type='l', main='Likelihood', xlab='Mean',ylab='Likelihood')
 abline(v=meanY)
 

 ### Bahyesian analysis
 priorDensity=dnorm(x=myGrid,mean=mu0,sd=sqrt(v0))

 rhs=sum(y)/varX+mu0/v0
 C=(length(y)/varX+1/v0)
 postMean=rhs/C

 postVar=1/C
 
 
 ## Comparing the likelihood function, prior and posterior density
 postDensity=dnorm(x=myGrid,mean=postMean,sd=sqrt(postVar))
 likelihood<-likelihood/max(likelihood)*max(postDensity)
 
 plot(postDensity~myGrid,type='l',col=2)
 lines(x=myGrid,y=likelihood,col=4)
 lines(x=myGrid,y=priorDensity)
 abline(v=c(mu0,meanY,postMean),col=c(1,4,2),lty=2)
 
```


#### (2) Composition Sampling

```R
 y=scan('~/GitHub/STT465/bmi.txt')
 n=length(y)
 
 DF0=4
 S0=var(y)*(DF0+2) # a bit of Bayesian 'cheating'; this gives E[var]=var(y)
 mu0=28
 k0=3
 meanY=mean(y)
 SS<-sum((y-meanY)^2)+ ((meanY-mu0)^2)*(n*k0/(k0+n))
 S=S0+ SS
 DF= n+DF0
 
 nSamps=100000
 mu=rep(NA,nSamps)
 varE=rep(NA,nSamps)
 
 for(i in 1:nSamps){
  varE[i]=S/rchisq(n=1,df=DF)
  v0=varE[i]/k0
  C=n/varE[i]+1/v0
  rhs<-meanY*n/varE[i] + mu0/v0
  sol=rhs/C
  mu[i]=rnorm(n=1,sd=sqrt(1/C),mean=sol) 
  #print(i)
 }
 mean(varE)
```


#### (3) Gibbs Sampler

```R
 mu_gibbs=rep(NA,nSamps+1)
 varE_gibbs=rep(NA,nSamps+1)
 
 mu_gibbs[1]<-mean(y)

 for(i in 2:(nSamps+1)){
  # sampling error variance conditional on the mean
   SSy=sum((y-mu_gibbs[i-1])^2) 
   S=SSy+S0
   varE_gibbs[i]=S/rchisq(n=1,df=DF)
  
  # sampling the mean conditional on the error variance
   C=n/varE_gibbs[i]+1/v0
   rhs<-meanY*n/varE_gibbs[i] + mu0/v0
   sol=rhs/C
   mu_gibbs[i]=rnorm(n=1,sd=sqrt(1/C),mean=sol) 
   #print(i)
 }
 
 mu_gibbs<-mu_gibbs[-1]
 varE_gibbs<-varE_gibbs[-1]
 
 mean(varE_gibbs)
 
 plot(density(mu),col=2)
 tmp=density(mu_gibbs)
 lines(x=tmp$x,y=tmp$y,col=4)
 
