### Examples MC Approximations (Ch. 4)


#### Example 1: Approximating the standard normal distribution

```R
  S=c(20,50,100,500,1000,5000)
  par(mfrow=c(2,3))
  z=seq(from=-3, to=3, by=.001)
  for(i in 1:6){
    x=rnorm(n=S[i])
    plot(density(x),xlim=c(-2.5,2.5),col=2,ylim=c(0,.5),main=paste0('S=',S[i]))
    lines(x=z,y=dnorm(z),col=4)
  }
```


#### Example 2: Computing the areas within a unit-square

```R
 S=10000 # number of MC samples

 # x and y are independent uniform RVs
 x=runif(S)
 y=runif(S)
 
 # area below the diagonal
 z=x>y
 mean(z)
 
 # set satysfying x<y^2
  z= x < y^2
  mean(z)
  
```


#### Example 3: Convergence and MC error

```R
 S=5e4
 x=rgamma(rate=2,shape=4,n=S)
 runMean=cumsum(x)/(1:length(x))
 plot(runMean,type='l',col=4) ;  abline(h=2,col=2)
 
 ## Let's now see the MC error on estimating the pth quantile
 prob=.99
 rate=2
 shape=4
 x=rgamma(rate=rate,shape=shape,n=S)
 estPercentile<-rep(NA,length=S-100)
 for(i in 1:length(estPercentile)){
    estPercentile[i]<-quantile(x=x[1:(100+i)],prob=prob)
    print(i)
 }
  plot(estPercentile,type='l',col=4); abline(h=qgamma(rate=rate,shape=shape,p=prob),col=2)
 
```

#### Example 4: Approximating the distribution of a deterministic function of a RV

```R
  ## Approximating the posterior density of log(lambda) in the gamma-poisson model
  
  ## Gamma Poisson Model
  N=100 # sample size
  lambda=4
  S=1e5 # # of MC samples
  
  y=rpois(n=N,lambda=lambda) # data
  # prior hyperparameters
    bPrior=5 # our 'prior sample size'
    aPrior=10 # our 'prior n*yBar'
    
  # posterior hyper-parameters
    aPost=aPrior+ sum(y)
    bPost=bPrior+N
    
  # The posterior density is rgamma(rate=bPost, shape=aPost,..)
  x<-rgamma(rate=bPost, shape=aPost,n=S) # samples from the posterior distribution
  mean(  x) #MC Estimate of the posterior mean
  aPost/bPost # true posterior mean
  mean(y)       # MLE
  aPiror/bPrior # prior mean
    
  ## Estimating the posterior density of log(lambda)
  y=log(x)
  plot(density(y))
  
```

## Example 5: approximating the marginal likelihood

```R
  N=100 # sample size
  lambda=4
  S=1e5 # # of MC samples
  
  y=rpois(n=N,lambda=lambda) # data
  
  # prior hyperparameters
    bPrior=5 # our 'prior sample size'
    aPrior=10 # our 'prior n*yBar'
  
  # we first sample lambda from the prior and then evaluat the conditional likelihood
  # averaging integrates that function.
    log_condLikelihood<-rep(NA,S)
    for(i in 1:S){
       lambda<-rgamma(rate=bPrior,shape=aPrior,n=1)
       tmp<-dpois(x=y,lambda=lambda,log=T)
       log_condLikelihood[i]=sum(tmp)
    }
```
