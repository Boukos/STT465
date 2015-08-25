**Sampling from binomial distirbution**

**Estimating the expected value and variance of an estimator & it's MSE**

```R
##
 theta=0.8
 n=10
##

 x=rbinom(n=n,size=1,prob=theta)
 mean(x)
```

** No let's play the sampling excercise**

```R
 nRep=1000
 thetaHat=rep(NA,nRep)

 for(i in 1:nRep){
   x=rbinom(n=n,size=1,prob=theta)
   thetaHat[i]=mean(x)
 }
 plot(hist(x,30)); abline(v=theta,col=4,lwd=2)

 #bias
 mean(thetaHat)-theta

 # variance
 var(thetaHat)

 # MSE
  Ex=mean(x)
  mean((x-Ex)^2)
```

