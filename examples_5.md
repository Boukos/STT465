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
 n<-100
 y=sample(x,size=n,replace=F)
 meanY=mean(y)


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
