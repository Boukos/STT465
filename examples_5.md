
 mu0=27
 v0=10
 x=scan('~/GitHub/STT465/bmi.txt')
 
 n<-10

 y=sample(x,size=n,replace=F)

 meanY=mean(y)
 varY=var(y)
 myGrid=seq(from=min(y),to=max(y),by=.1)

 priorDensity=dnorm(x=z,mean=mu0,sd=sqrt(v0))

 rhs=sum(y)/varY+mu0/v0
 C=(length(y)/varY+1/v0)
 postMean=rhs/C

 postVar=1/C
 postDensity=dnorm(x=z,mu=postMean,sd=sqrt(postVar))
 postDensity=dnorm(x=z,mean=postMean,sd=sqrt(postVar))

