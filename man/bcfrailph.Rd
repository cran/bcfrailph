% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcfrailph.R
\name{bcfrailph}
\alias{bcfrailph}
\title{Bivariate correlated frailty model with Proportional hazard.}
\usage{
bcfrailph(
  formula,
  data,
  frail_distrn = c("gamma", "lognormal"),
  initfrailp = NULL,
  control,
  ...
)
}
\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the Surv function.}

\item{data}{A dataframe contain survival time, censor, covariate etc with data in columns.}

\item{frail_distrn}{A type of frailty distribution to be used in fit. Either gamma or lognormal. The default is gamma.}

\item{initfrailp}{Initial estimates for the frailty parameters. The default is c(0.5,0.5).}

\item{control}{Arguments to control the fit. The default is \code{\link{bcfrailph.control}}.}

\item{...}{further arguments}
}
\value{
An object of  that contains  the following components.
\itemize{
  \item \code{coefficients} - {A vector of estimated Covariate coefficients.}
  \item \code{frailparest} - {A vector of estimated Frailty parameters i.e. frailty variance and correlation.}
  \item \code{vcov2}- {Variance Covariance matrix of the Estimated Covariate coefficients obtained from the observed information matrix.}
  \item \code{vcovth2}-Variance Covariance matrix of the Estimated Frailty parameters obtained from the observed information matrix of the mariginal likelihood.
  \item \code{loglilk0}- Log likelihood of without frailty model.
  \item \code{loglilk}-Log likelihood of Cox PH model with frailty.
  \item \code{Iloglilk}- Log likelihood of with frailty model after integrating out the frailty term.
  \item \code{cbashaz}- array containing Cummulative baseline hazard.
  \item \code{X}-{Matrix of observed covariates.}
  \item \code{time}-{the observed survival time.}
  \item \code{censor}-{censoring indicator.}
  \item \code{resid}-{the martingale residuals.}
  \item \code{lin.prid}-{the vector of linear predictors.}
  \item \code{frail}-{estimated Frailty values.}
  \item \code{stderr}-{A vector containing the Standard error of the Estimated parameters both covariate coefficients and  frailty parameters.}
  \item \code{iteration}-{Number of outer iterations.}
  \item \code{e.time}-{the vector of unique event times.}
  \item \code{n.event}- {the number of events at each of the unique event times.}
  \item \code{converg}-  {TRUE if converge, FALSE otherwise.}
  }
}
\description{
Fit a semiparametric Bivariate correlated frailty model with Proportional Hazard structure.
}
\note{
Parameters of Bivariate correlated gamma frailty model was estimated basically using the EM-approach proposed by Iachine, I. A. (1995) with modifications.
The main modification that made on the original EM-approach was similar to the modification made on EM approach for univariate gamma frailty model by Duchateau and Janssen (2008).
This means following more or less similar procedure as Duchateau and Janssen (2008), frailty parameters are estimated from the marginal log likelihood function.
The results of both EM- approach and the modified EM- approach are similar. The difference is that the modified one is much faster.

Bivariate correlated gamma frailty model is constructed in a way that the estimated parameters except the correlation are more or less similar with the estimates
obtained from univariate gamma frailty model. Thus, the standard errors of the estimated covariate coefficients and the frailty variance parameter is obtained using
the standard errors estimation approach for univariate gamma frailty model given in Klein and Moeschberger (2003). The standard error of the estimated frailty
correlation parameter is obtained from the observed information matrix of the marginal log likelihood. Simulation study showed the estimation approach is reasonably good.

Parameters of Bivariate correlated lognormal frailty model was using penalized likelihood approach used by Ripatti and Palmgren (2000).
}
\examples{
set.seed(24)
simdata<-simbcfrailph(p.size=300, c.rate= c(0.3),fraildistrn=c("gamma"),frail.par=c(0.5,0.5),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

fitbcfrailph=bcfrailph(Surv(time,censor)~ X1+frailty(PID) ,data=dataa,frail_distrn=c("gamma"))
fitbcfrailph


\dontshow{
set.seed(5)
simdata<-simbcfrailph(p.size=300, c.rate= c(0.2),fraildistrn=c("lognormal"),frail.par=c(0.5,0.6),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

fitbcfrailph=bcfrailph(Surv(time,censor)~ X1++frailty(PID) ,
data=dataa,frail_distrn=c("lognormal"),control=bcfrailph.control(max.iter2=5))
fitbcfrailph

set.seed(24)
simdata<-simbcfrailph(p.size=300, c.rate= c(0.3),fraildistrn=c("gamma"),frail.par=c(0.5,0.5),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data

fitbcfrailph=bcfrailph(Surv(time,censor)~ X1+frailty(PID) ,data=dataa,frail_distrn=c("gamma"))
fitbcfrailph
}

\donttest{
# now for lognormal
#Weibull baseline hazard with parameters shape= 5 and scale=0.1.
#a dataset with 300 pairs. Lognormal frailty distribution and the frailty parameters are taken to
#be variance=0.5 and rho=0.6. One binomial B(1,0.5) with regression coefficient 4.
#Each observed covariate for the two individuals in a
#pair is taken to be independent and 20 percent of the observations are censored.
# simulate the data set

set.seed(5)
simdata<-simbcfrailph(p.size=300, c.rate= c(0.2),fraildistrn=c("lognormal"),frail.par=c(0.5,0.6),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(2),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data ## the simulated data set

#fit
fitbcfrailph=bcfrailph(Surv(time,censor)~ X1++frailty(PID) ,data=dataa,frail_distrn=c("lognormal"))
fitbcfrailph
# the output looks like
#   Call:
#   bcfrailph(formula = Surv(time, censor) ~ X1, data = dataa, frail_distrn = c("lognormal"))
#
#   n=  600 and number of events= 466
#
#   Regression Coefficients:
#   Estimate  StdErr     se2 z.value   p.value
#   4.15678 0.22724 0.18627  18.292 < 2.2e-16 ***
#   ---
#   Frailty Distribution:Bivariate Correlated  lognormal
#   Variance of random effect = 0.7578053 ( 0.06592505 )
#   Correlation Estimate of random effects = 0.5205348 ( 0.06406965 )
#   Log likelihood with frailty = -2216.765
#   Log likelihood without frailty= -2455.313

# gamma fit in uncensored data

# simulate the data set
set.seed(3)
simdata<-simbcfrailph(p.size=300, c.rate= c(0),fraildistrn=c("gamma"),frail.par=c(0.5,0.6),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(1.5),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data ## the simulated data set

#fit
fitbcfrailph=bcfrailph(Surv(time,censor)~ X1+cluster(PID) ,data=dataa,frail_distrn=c("gamma"))
fitbcfrailph

# lognormal fit in uncensored data
set.seed(4)
simdata<-simbcfrailph(p.size=300, c.rate= c(0),fraildistrn=c("lognormal"),frail.par=c(0.5,0.6),
bhaz.arg=list(distrn = c("weibull"),shape =c(5), scale = c(0.1)),
covar.arg=list(coefs=c(4),types = c("B"),size=c(1),prob=c(0.5)))
dataa<-simdata$data ## the simulated data set

#fit
fitbcfrailph=bcfrailph(Surv(time,censor)~ X1++frailty(PID) ,data=dataa,frail_distrn=c("lognormal"))
fitbcfrailph

# let us try with two covariates in detail
####try bcfrailph on the following SIMULATED DATASET
####Gompertz baseline hazard with parameters shape= 0.1 and scale=0.0001.
####a dataset with 600 pairs. The frailty parameters are taken to
####be variance=0.25 and rho=0.5. One continuous (U [0,1]) and one
####categorical (Binomial (1,0.4)) with regression coefficients beta1=3 and beta2=-1
####Each observed covariates for the two individuals in a
####pair is taken to be independent and all are uncensored.

n=600; n1=n*2   ### 600 pairs
indic1=2*array(1:n)-1;indic2=2*array(1:n)
PID=1;e1=array(1:n);PID[indic1]=e1;PID[indic2]=e1  ### PID is cluster indicator

X11<-runif(n,  min=0, max=1) ;X12<-runif(n,  min=0, max=1)
#####covariate 1 for the first and second individuals
X21<-rbinom(n, size=1, prob=0.4) ;X22<-rbinom(n, size=1, prob=0.4)
#####covariate 2 for the first and second individuals
# gamma frailty variabbles  with sigma^2 =0.25 and roh=0.5
lam=1/0.25;k0=0.5/0.25;k1=k2=0.5/0.25
y0=rgamma(n,shape=k0,scale=1/lam) ;y1=rgamma(n,shape=k1,scale=1/lam)
y2=rgamma(n,shape=k2,scale=1/lam)
z1=(y0+y1);z2=(y0+y2)#frailty variables
#  with this set up both z1 and z2 are gamma
#######distributed frailty variables with mean 1 and variance 0.25
# and the correlation between z1 and z2 is 0.5
#
u1<-runif(n,  min=0, max=1)
u2<-runif(n,  min=0, max=1)
# survival times
t1<- 1/0.1*log(1-0.1*log(u1)/(0.0001*exp(3*X11-1*X21)*z1))
t2<- 1/0.1*log(1-0.1*log(u2)/(0.0001*exp(3*X12-1*X22)*z2))
# let us organize the data set
time=X1=X2=1
time[indic1]=t1;time[indic2]=t2
X1[indic1]=X11;X1[indic2]=X12
X2[indic1]=X21;X2[indic2]=X22
censor=rep(1,n1)
#data set
dataa <- data.frame(time=time, X2=X2, X1=X1,censor=censor,ID=ID)
#
#fit
fitbcfrailph=bcfrailph(Surv(time,censor)~ X1+X2+frailty(PID) ,data=dataa)
fitbcfrailph

## one can set the initial parameter for the frailty parameters
## the default is initfrailp = c(0.5,0.5)
fitbcfrailph=bcfrailph(Surv(time,censor)~ X1+X2+frailty(PID),data=dataa,initfrailp = c(0.1,0.5))
fitbcfrailph

# since the estimated covariate coeficients and frailty variance
#parameter are more or less similar with coxph with univariate gamma frailty

IID=array(1:nrow(dataa))# individual id
cphfit <- coxph(Surv(time, censor, type = "right") ~ X1+X1+X2+frailty(IID),data =  dataa)
cphfit

# Not run

#if covariates are not included
fitmoe=bcfrailph(Surv(time,censor)~0,data=dataa,frail_distrn=c("lognormal"))
fitmoe
fitmoe=bcfrailph(Surv(time,censor)~1,data=dataa,frail_distrn=c("lognormal"))
fitmoe

#if fraility id is not specified correctly
#or if it is not specified in a way that it indicates pairs.
ID=array(1:nrow(dataa))# this is not pair id rather it is individual id.
fitmoe=bcfrailph(Surv(time,censor)~ X1+frailty(ID),data=dataa,frail_distrn=c("lognormal"))
fitmoe
fitmoe=bcfrailph(Surv(time,censor)~ X1+cluster(ID),data=dataa,frail_distrn=c("lognormal"))
fitmoe

# if control is not specified correctly.
# if one needs to change only max.iter to be 100,

fitmoe=bcfrailph(Surv(time,censor)~ X1+frailty(PID),data=dataa,control=c(max.iter=100))
fitmoe

#the correct way is
fitmoe=bcfrailph(Surv(time,censor)~ X1+frailty(PID),data=dataa,
control=bcfrailph.control(max.iter=100))
fitmoe

#if initial frailty parameters are in the boundary of parameter space
fitmoe=bcfrailph(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0.2,1))
fitmoe
fitmoe=bcfrailph(Surv(time,censor)~ X1,data=dataa,initfrailp=c(0,0.1))
fitmoe

#if a frailty distribution other than gamma and lognormal are specified

fitmoe=bcfrailph(Surv(time,censor)~ X1,data=dataa,,frail_distrn=c("exp"))
fitmoe
# End Not run
}

}
\references{
Duchateau, L., Janssen, P. (2008) The Frailty Model. Springer, New York.

Iachine, I. A. (1995). Correlated frailty concept in the analysis of bivariate survival data. Bachelor project, Odense University, Department of Mathematics and Computer Science, Denmark.

Klein, J. P., and Moeschberger, M. L. (2003), Survival analysis: techniques for censored and truncated data, New York: Springer.

Rippatti, S. and Palmgren, J (2000). Estimation of multivariate frailty models using penalized partial likelihood. Biometrics, 56: 1016-1022.
}
\seealso{
\code{\link{bcfrailph.control}},\code{\link{simbcfrailph}}
}
