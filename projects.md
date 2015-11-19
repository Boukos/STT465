## Ideas for Final Project
-----------------------------------------------------------------------------------------------------	
####Project 1	

   **Title**: Multiple-Linear Regression with Binary Outcomes

   **Programming task**: Implement a Gibbs sampler for a multiple-linear regression model and extend it to handle binary outcomes using the probit link.

   **Analysis**:  Multiple linear regession of Gout (Yes/No) on: Uric Acid, Creatinine, BMI, Glucose, HDL, LDL, Triglycerides, Race, Sex and Age. 
    
   **Expected outcomes**: 
       (1) Descriptive statistics for each variable and for the response against each predictor.
       (3) Estimates of effects and 95% CI derived from GLM (i.e., maximum likelihood).
       (4) Estimates of effects and 95 posterior crediblity regions derived froma Bayesian model.
       (5) Estimated posterior means and estimated 95% credibility regions for the change in risk of developing Gout for the following comparisons:
              - Male versus Female
              - Black versus White
              - 2 point increase in Uric Acid
        (6) A  2 paragraph summary of your findings.
        (7) Appendix, including: (a) the code you use, and (b) convergence diagnosis (e.g., trace plots, density plots, MCErrors, etc.) for the Bayesian analyses.
        
Notes: for the Bayesian analysis, treat effects as 'Fixed' and run a sufficiently long chain. For the glm analysis, a sample code is provided. 
   
  
```R
 DATA=read.table('~/Dropbox/STT_465_FALL_2015/gout.txt',header=T,as.is=T)
 y=ifelse(DATA$Gout=='Y',1,0)
 fmGLM=glm( y~Sex+Age+BMI+Race+Sex+UricAcid+Glucose+Creatinine+SBP+LDL+HDL+SBP+Triglycerides,
            data=DATA,family=binomial(link=probit))
 summary(fmGLM)
```
              
   **[Data-link](https://www.dropbox.com/s/ho3p0uwohjnoln3/gout.txt?dl=0)**

-----------------------------------------------------------------------------------------------------	
####Project 2	

**Title**: Multiple-linear regression with censoring.

**Programming task**: Implement a Gibbs sampler for a linear regression model and extend it to handle right, left and interval censoring.

**Analysis**:  using a real data set (this will be provided) analyze it using OLS (lm) and your software, accounting and ignoring censoring.

**Expected outcome**: a comparison of your results with those provided by lm().

**[Data](https://www.dropbox.com/s/1rw7s4z1ta3kehy/DATA_STT465.RData?dl=0)**

-----------------------------------------------------------------------------------------------------	

####Project 3	
									

**Title**: High dimensional regression.

**Programming task**: implement a Gibbs sampler for a linear regression model and extend it to handle two sets of predictors, one will be treated as fixed and the second one as random.

**Analysis**:   using a real data set (this will be provided) analyze it using: BGLR() and your software.

**Expected outcome**: a comparison of your results Compare your results with those obtained with BGLR.
```R


library(BGLR)
data(mice)
nFolds=5
y=scale(mice.pheno$Obesity.BMI)

set.seed(12345)
fold=sample(1:nFolds,size=length(y),replace=T)


ETA=list(   fixed=list(~GENDER+Litter,data=mice.pheno,model='FIXED'),
            cage=list(~cage,data=mice.pheno,model='BRR'),
            markers=list(X=scale(mice.X), model='BRR')

         )


COR=matrix(nrow=nFolds,ncol=2)

for(i in 1:nFolds){
  print(i)
  tst<-which(fold==i)
  yNA=y
  yNA[tst]<-NA
  fm0=BGLR(y=yNA,ETA=ETA[c(1:2)],nIter=12000)
  fmA=BGLR(y=yNA,ETA=ETA,nIter=12000)
  
  COR[i,1]<-cor(y[tst],fm0$yHat[tst])
  COR[i,2]<-cor(y[tst],fmA$yHat[tst])
  
}

```

-----------------------------------------------------------------------------------------------------	
####Project 4	

**Title**: non-parametric regression and automatic knot selection using Bayesian models

**Programming task**: implement a Gibbs Sampler for a linear regression model, using a simulated data set (this will be provided) compare the use of natural splines with different degree of freedom and with the Bayesian approach.

**Analysis**:   analyze the simulated data set using splines with variable degree of freedom using lm() and the Bayesian software that you developed. 

**Expected outcome**: a comparison of your results  with those obtained with lm().


-----------------------------------------------------------------------------------------------------	
