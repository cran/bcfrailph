#'
#' @name simbcfrailph
#' @title Simulate data from bivariate correlated frailty models.
#' @description Simulate data from bivariate correlated gamma or lognormalfrailty models with one covariate.
#'
#' @param p.size pair size.
#' @param c.rate censored rate. The default is zero..
#' @param fraildistrn A type of frailty distribution to be used. Either gamma or lognormal.
#' @param frail.par vector of frailty parameters, variance and correlation respectively. The default is c(0.5,0.25) meaning variance 0.5 and correlation 0.25.
#' @param bhaz.arg is a \code{\link{list}} i.e,\code{list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01), rate = c(0.5))}.\code{distrn} is the type of baseline hazard to be used. weibull, gompertz and exponential is possible. \code{shape}, \code{scale} and \code{rate} are the parameters of the coresponding baseline hazard parameters. rate needs to be specified if the baseline is exponential.
#' @param covar.arg is a \code{\link{list}} i.e,\code{coefs} is covariate coefficient, \code{types} is the type of covariate to be used. \code{B} is for binomial and \code{U} is for uniform.\code{size} and \code{prob} needs to be specified if the covariate is binomial and if the covariate is uniform, then \code{min} and \code{max} should be specified.
#'
#' @return
#' An object of class \code{simbcfrailph} that contain the following:
#'
#' \itemize{
#'
#' \item{\code{data}}  {A data frame i.e, the simulated data set. IID is individual Id, PID is pair ID, time is the simulated survival time, censor is censoring indicator and X1 denote the simulated covariate.}
#'
#' \item{\code{numberofpair}}  {The specified number of pairs.}
#'
#' \item{\code{censoredrate} } {The specified censored rate.}
#'
#' \item{\code{fraildist} } {The specified frailty distribution.}
#'
#' \item{\code{frailpar}} {The specified frailty parameters.}
#' }
#'
#'
#' @export simbcfrailph
#'
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rbinom
#' @importFrom stats quantile
#'
#' @seealso \code{\link{bcfrailph}}
#'
#' @examples
#' simdata<-simbcfrailph(p.size=1000, c.rate= c(0.2),fraildistrn=c("gamma"),frail.par=c(0.5,0.5),
#' bhaz.arg=list(distrn = c("gompertz"),shape =c(3), scale = c(0.1)),
#' covar.arg=list(coefs=c(1),types = c("U"),min=0,max=1))
#' simdata
#'
#' \donttest{
#' #Let us simulate a data set with the following parameters
#' #weibull baseline hazard with parameters shape= 2.5 and scale=0.01.
#' #a dataset with 1000 pairs. Frailty distribution is gamma
#' # and the frailty parameters are taken to
#' #be variance=0.4 and correlation =0.6. One binomial covariate i.e,
#' #(Binomial (1,0.5)) with regression coefficient 0.5.
#' #Each observed covariate for the two individuals in a
#' #pair is taken to be independent and 20 percent of the observations are censored.
#'
#' #simulate the data set
#' #first for gamma
#'
#' simdata<-simbcfrailph(p.size=1000, c.rate= c(0.2),fraildistrn=c("gamma"),frail.par=c(0.4,0.6),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01)),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#' simdata  # a list contain the parameters set by user and the simulated data set
#'
#' #to extract the simulated data set
#' dataa<-simdata$data ## the simulated data set
#' dataa[1:4,] # the first four rows looks like
#'
#' #    IID PID     time censor X1
#' #  1   1   1 1.793324      0  0
#' #  2   2   1 5.245163      1  1
#' #  3   3   2 6.730729      1  0
#' #  4   4   2 4.963159      1  0
#'
#' # IID is individual indicator
#' # PID is pair indicator
#' # time is the simulated survival time
#' # censor is the simulated censoring indicator
#' # X1 is the simulated covariate
#'
#' # if lognormal frailty is desired
#'
#' simdata<-simbcfrailph(p.size=1000, c.rate= c(0.2),fraildistrn=c("lognormal"),frail.par=c(0.4,0.6),
#' bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01)),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#' dataa<-simdata$data ## the simulated data set
#'
#' # Not run
#' # if p.size, pair size missed
#' simdata<-simbcfrailph( c.rate= c(0.2),fraildistrn=c("gamma"),frail.par=c(0.4,0.6),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#'
#' # if frailty distribution other than gamma and lognormal specified
#' simdata<-simbcfrailph(p.size=100, c.rate= c(0.2),fraildistrn=c("exp"),frail.par=c(0.4,0.6),
#' covar.arg=list(coefs=c(0.5),types = c("B"),size=1,prob=0.5))
#' # End Not run
#' }
#'
simbcfrailph<-function(p.size, c.rate= c(0),fraildistrn,frail.par=c(0.5,0.25),
bhaz.arg=list(distrn = c("weibull"),shape =c(2.5), scale = c(0.01), rate = c(0.5)),
covar.arg=list(coefs=c(0.5),types = c("B","U"),size=1,prob=0.5,min=0,max=3)){
Call <- match.call()
if (missing(p.size)) {stop("number of pairs must be supplied")}
p.size<-round(p.size,digits=0)
p.size=p.size[1];c.rate=c.rate[1]
if (p.size<=0) {stop(" number of pairs must be positive")}
if ((c.rate<0) |(c.rate>=1)) {stop(" Censoring rate must be between 0 and 1")}
if (missing(fraildistrn)) {fraildistrn<-c("gamma")}
if(length(fraildistrn)>1){stop("simbcfrailph doesn't support more than one frailty distributions")}
if(!(c(fraildistrn)%in%c("gamma","lognormal"))){
stop("simbcfrailph only support gamma or lognormal frailty distributions")}
if (length(frail.par)<2) {stop(" frailty parameters i.e. variance and correlation must be supplied")}
if (length(frail.par)>2) {frail.par<-frail.par[1:2]}
if(any(frail.par<=0)){stop("At least one invalid frailty parameter")}
if (frail.par[2]>=1) {stop(" invalid frailty correlation parameter")}
bhaz.arg$distrn=bhaz.arg$distrn[1];bhaz.arg$shape=bhaz.arg$shape[1]
bhaz.arg$scale=bhaz.arg$scale[1];bhaz.arg$rate=bhaz.arg$rate[1]
covar.arg$coefs=covar.arg$coefs[1];covar.arg$types=covar.arg$types[1]
covar.arg$size=covar.arg$size[1];covar.arg$prob=covar.arg$prob[1]
covar.arg$min=covar.arg$min[1];covar.arg$max=covar.arg$max[1]
if(!(c(bhaz.arg$distrn)%in%c("weibull","gompertz","exponential"))){
stop("simbcfrailph only support weibull gompertz or exponential baseline hazards")}
if((pmatch(c(bhaz.arg$distrn),c("weibull","gompertz","exponential"))<3)){
if(any(c(bhaz.arg$shape,bhaz.arg$scale)<=0)){stop(" invalid baseline hazard parameters")}}
if((pmatch(c(bhaz.arg$distrn),c("weibull","gompertz","exponential"))==3)){
if(bhaz.arg$rate<=0){stop(" invalid rate parameter for exponential baseline hazard")}}
if(!(c(covar.arg$types)%in%c("B","U"))){
stop("incorrect covariates type and types should be either B or U")}
dataa<-.gendata1(p.size,c.rate,fraildistrn,frail.par,bhaz.arg,covar.arg)
gendat <-list(numberofpair=p.size,censoredrate= c.rate,
fraildist=fraildistrn,frailpar=frail.par,basehazdis=bhaz.arg$distrn,
covartype=covar.arg$types,covarcoef=covar.arg$coefs)
gendat$data<-dataa
gendat$call <- match.call()
class(gendat) <- c("simbcfrailph",class(gendat))
gendat}


