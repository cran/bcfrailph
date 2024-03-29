% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitbccv.lognsp.R
\name{fitbccv.lognsp}
\alias{fitbccv.lognsp}
\title{Bivariate correlated lognormal frailty model fitting function.}
\usage{
fitbccv.lognsp(X, Y, initfrailp, control = bcfrailph.control())
}
\arguments{
\item{X}{Matix of predictors. This should not include an intercept.}

\item{Y}{a Surv object containing 2 columns (coxph.fit).}

\item{initfrailp}{Initial estimates for the frailty parameters. If not specified, initial frailty variance will be obtained from coxph with univariate lognormal frailty model and for correlation c(0.5) will be used.}

\item{control}{Arguments to control the fit. The default is \code{\link{bcfrailph.control}}.}
}
\value{
An object of  that contains  the following components.
\itemize{
  \item \code{coefficients} - {A vector of estimated Covariate coefficients.}
  \item \code{frailparest} - {A vector of estimated Frailty parameters i.e. frailty variance and correlation.}
  \item \code{stderr}-{A vector containing the Standard error of the Estimated parameters both covariate coefficients and  frailty parameters.}
  \item \code{loglilk0}- Log likelihood of without frailty model or loglik of coxph fit.
  \item \code{loglilk}-Log likelihood of Cox PH model with frailty.
  \item \code{Iloglilk}- Log likelihood of with frailty. For gamma fit it is I-likelihood or the likelihood after integrating out the frailty term.For lognormal fit it is the approximate likelihood.
  \item \code{bhaz}- an array containing unique event times and estimated baseline hazard.
  \item \code{X}-{Matrix of observed covariates.}
  \item \code{time}-{the observed survival time.}
  \item \code{censor}-{censoring indicator.}
  \item \code{resid}-{the martingale residuals.}
  \item \code{lin.prid}-{the vector of linear predictors.}
  \item \code{frail}-{estimated Frailty values.}
  \item \code{iteration}-{Number of outer iterations.}
  \item \code{e.time}-{the vector of unique event times.}
  \item \code{n.event}- {the number of events at each of the unique event times.}
  \item \code{convergence}-{an indicator, 1 if converge and 0 otherwise.}
  \item \code{history}-{an array containing records of estimates and other information on each iterations.}
  }
}
\description{
Semi-parametric Bivariate correlated lognormal frailty model fitting function.
}
\note{
This function is important especially for simulation studies as it reduced checking time.
Parameters of Bivariate correlated lognormal frailty model is based on the penalized partial likelihood approach by Rippatti and Palmgren (2000).
}
\examples{
set.seed(18)
simdata<-simbcfrailph(psize=100, cenr= c(0.2),beta=c(1,-0.7,0.5),frailty=c("lognormal"),
frailpar=c(0.5,-0.25),bhaz=c("exponential"),
bhazpar=list(scale = c(0.1)),covartype= c("N","N","B"),
covarpar=list(fargs=c(0,0,1),sargs=c(1,1,0.5)),comncovar=2)
Y<-simdata$Y;X<-simdata$X

bcspfit<-fitbccv.lognsp(X=X,Y=Y,initfrailp=NULL)
bcspfit$coef 
bcspfit$frailpar

}
\references{
Rippatti, S. and Palmgren, J (2000). Estimation of multivariate frailty models using penalized partial likelihood. Biometrics, 56: 1016-1022.
}
\seealso{
\code{\link{bcfrailph}}
}
