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


#### Example 2: Compute the area in the unit square that satisfies y>=x^2

```R

```

#### Example 3: 

```R

```

#### Example 4: 

```R

```

