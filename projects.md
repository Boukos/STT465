## Ideas for Final Project
-----------------------------------------------------------------------------------------------------	
####Project 1	

   **Title**: Multiple-Linear Regression with Binary Outcomes

   **Programming task**: Implement a Gibbs sampler for a multiple-linear regression model and extend it to handle binary outcomes using the probit link.

   **Analysis**:  Multiple linear regession of Gout (Yes/No) on: Uric Acid, Creatinine, BMI, Glucose, HDL, LDL, 
                Triglycerides, Race, Sex and Age. 
    
   **Expected outcomes**: 
   
	- Descriptive statistics for each variable and for the response against each predictor.
	- Estimates of effects and 95% CI derived from GLM (i.e., maximum likelihood).
	- Estimates of effects and 95 posterior crediblity regions derived froma Bayesian model.
	- Estimated posterior means and estimated 95% credibility regions for the change in risk of developing Gout for:
        	- Male versus Female
        	- Black versus White
        	- 2 point increase in Uric Acid
	- A  2 paragraph summary of your findings
	- Appendix, including: (a) the code you use, and (b) convergence diagnosis (e.g., trace plots, density plots, MC-SEs, etc.) for the Bayesian analyses.
        
Note: for the Bayesian analysis, treat effects as 'Fixed' and run a sufficiently long chain. For the glm analysis, a sample code is provided below. 
   
  
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

The main goal is to assess wheather gene expression information derived from the tumor is associated to expected years of life after diagnosis of Glioblastoma. The [data set](https://www.dropbox.com/s/1rw7s4z1ta3kehy/DATA_STT465.RData?dl=0) provided contains survival information (days to last folowup and survival status at last follow up), a series of clinical predictors and 30 principal components derived from gene expression profiles assessed at the tumor cell.

**Programming task**: Implement a Gibbs sampler for a linear regression model and extend it to handle right censoring.

**Analysis**: Regression of survival time of Glioblastoma patients on clinical covariates and principal components derived from gene expression.

**Expected outcomes**

	- Estimated effects from a baseline model obtaiend by regressing survival time on clinical and demographic covariates. The report should include both results from g maximum likelihood (sruvreg, see sample code below) and Bayesian analysis. For maximum likelihood report parameter estimates and 95% confidence intervals. For Bayesian report estimated posterior means and estiamted 95% crediblity regions.
	- Extend the baseline model by adding the random effects of the 30 principal components provided in the dataset. Estimate, using a  Bayesian model, the proportion of variance explained by the principal components, and assess which ones seem to have sronger association with survival.
	- An appendix with: (a) the code you use for analysis, and (b) convergence diagnostics for Bayesian analysis (trace plots, estimated MC-SEs, etc.).

```R
  load('~/Dropbox/STT_465_FALL_2015/TCGA_GB/DATA_STT465.RData')
  library(survival)
  y<-Surv(time=log(Y$days_to_last_followup),event=Y$death)
  fmML=survreg(y~race+gender+initial_pathologic_diagnosis_method+age_group10,
               data=Y,dist='gaussian')
  summary(fmML)
  
```

-----------------------------------------------------------------------------------------------------	

####Project 3	
									

**Title**: High dimensional regression.

**Programming task**: implement a Gibbs sampler for a multiple linear regression model with arbitrary sets of effets, and extend it to handle missing values on the response.

**Analysis**:  Compare, using cross-validation, a high-dimensional regression accounting for systematic and environmental effects (M0, regression of BMI on sex, litter size and cage) with one that includes both environmental and genetic effects (M1, regression of BMI on sex, litter size, cage and genetic markers).

**Expected outcomes**: 

	- A report of estimates of variance components (variance of random effects and error variances) and 95% posterior credibility regions for M0 and M1.
	- A comparison of prediction accuracy (cross-validation prediction correlation) of models M1 and M2 derived from a 5 fold cross-validation.
	- A 2 paragraph summary statment sumarizing your findings.
	- An appendix including: (a) convergence diagnostics (trace plots, estimates of MC SEs, etc.), and (b) the code used to carry out analyses.

Note: in your analyses treat sex and litter size as 'fixed effects' and the effects of cages and of DNA markers as random. Assign two different variance components for markers and cage.

The code below illustrate how to obtain the data from the BGLR package.

```R

library(BGLR)
data(mice)
nFolds=5
y=scale(mice.pheno$Obesity.BMI)
set.seed(12345)
fold=sample(1:nFolds,size=length(y),replace=T) # this vector randomly assigns each mice to a fold

XF=as.matrix(model.matrix(~GENDER+Litter,data=mice.pheno))[,-1]
XCage<-as.matrix(model.matrix(~cage-1,data=mice.pheno))
XMarkers<-scale(mice.X)/sqrt(ncol(mice.X))

COR=matrix(nrow=nFolds,ncol=2)
for(i in 1:nFolds){
  print(i)
  tst<-which(fold==i) # this determines which entries of y are used for testing in the ith fold.
  yNA=y
  yNA[tst]<-NA
  fm0=BGLR(y=yNA,ETA=list(list(X=XF,model='FIXED'),list(X=XCage,model='BRR')),nIter=12000)
  fmA=BGLR(y=yNA,ETA=list(list(X=XF,model='FIXED'),list(X=XCage,model='BRR'),list(X=XMarkers,model='BRR')),nIter=12000)
  COR[i,1]<-cor(y[tst],fm0$yHat[tst])
  COR[i,2]<-cor(y[tst],fmA$yHat[tst])
}

```

-----------------------------------------------------------------------------------------------------	
####Project 4	

**Title**: Non-parametric regression and automatic knot selection using Bayesian models

**Programming task**: implement a Gibbs Sampler for a linear regression model, using a simulated data set (this will be provided) compare the use of natural splines with different degree of freedom and with the Bayesian approach.

**Analysis**:   analyze the simulated data set using splines with variable degree of freedom using lm() and the Bayesian software that you developed. 

**Expected outcome**: a comparison of your results  with those obtained with lm().


-----------------------------------------------------------------------------------------------------	
