###Sampling from binomial distirbution

**Estimating the expected value and variance of an estimator & it's MSE**

```R
##
 theta=0.9
 n=35
##

 x=rbinom(n=n,size=1,prob=theta)
 mean(x)
```

**Now let's play the (frequentist) sampling excercise**

```R
 nRep=1000
 thetaHat=rep(NA,nRep)

 for(i in 1:nRep){
   x=rbinom(n=n,size=1,prob=theta)
   thetaHat[i]=mean(x)
 }
 plot(hist(thetaHat)); abline(v=theta,col=4,lwd=2)

 #bias
 mean(thetaHat)-theta

 # variance
 var(thetaHat)

 # MSE
  Ex=mean(thetaHat)
  mean((thetaHat-Ex)^2)
```

Note: 

 - Above we have computed Monte Carlo (MC) estimates of the expected value and varianc of the estimator. The MC error depends on the number of MC samples (nRep, above).

 - In this simple example we can actually compute the mean and variance of the estimator analythically.

