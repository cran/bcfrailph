
##########
.expfrailfun <- function(x,k0,k1,k2){
q1<-q2<-q3<-q4<-p1<-p2<-p3<-p4<-su<-NULL
B0<-B1<-B2<-frai<-NULL
if((x[1]==0) & (x[2]==0)){
B0<-(k0/x[5])
B1<-(k1/x[3]);B2<-(k2/x[4])
frai<-c(B0,B1,B2)}
if((x[1]==1) & (x[2]==0)){
q1=(k0/x[5]); q2=(k1/x[3])
p1=(q1/(q1+q2)); p2=(q2/(q1+q2))
B0=((p1*(k0+1)+p2*k0)/x[5])
B2=(k2/x[4])
B1=((p1*k1+p2*(k1+1))/x[3])
frai<-c(B0,B1,B2)}
if((x[1]==0) & (x[2]==1)){
q1=(k0/x[5]); q2=(k2/x[4])
p1=(q1/(q1+q2)); p2=(q2/(q1+q2))
B0=((p1*(k0+1)+p2*k0)/x[5])
B1=(k1/x[3])
B2=((p1*k2+p2*(k2+1))/x[4])
frai<-c(B0,B1,B2)}
if((x[1]==1) & (x[2]==1)){
q1=((k0*(k0+1))/x[5]^2); q2=((k0*k1)/(x[5]*x[3]))
q3=((k0*k2)/(x[5]*x[4])); q4=((k1*k2)/(x[3]*x[4]))
su=(q1+q2+q3+q4)
p1=(q1/su); p2=(q2/su); p3=(q3/su); p4=(q4/su)
B0=((p1*(k0+2)+(p2+p3)*(k0+1)+p4*k0)/(x[5]))
B1=(((p1+p3)*k1+(p2+p4)*(k1+1))/(x[3]))
B2=(((p1+p2)*k2+(p3+p4)*(k2+1))/(x[4]))
frai<-c(B0,B1,B2)}
frai
}


.Expefrail=function(newtht,cen1,cen2,H1,H2){
a<-newtht[1];R<-newtht[2]
k0<-(R/a);k1<-k2<-(1-R)/a
x<-cbind(cen1,cen2,(1/a+H1),(1/a+H2),(1/a+H1+H2))
frailcomp <- apply(x,1,.expfrailfun,k0=k0,k1=k1,k2=k2)
frailcomp =t(as.matrix(frailcomp))
z1=c(frailcomp[,1])+c(frailcomp[,2])
z2=c(frailcomp[,1])+c(frailcomp[,3])
list(frail1=log(z1),frail2=log(z2))
}

.SEgammacv=function(sig,n_eve,etime,nonzero_h0,exp_bet_xo,Ho,ho,censoro,timeo,Xo){
a<-as.vector(sig);n_eve<-as.vector(n_eve);etime<-as.vector(etime)
nonzero_h0<-as.vector(nonzero_h0);EXB<-exp_bet_xo
H=Ho;h0=ho;di=censoro;t0=timeo
Xo <- as.matrix(Xo)
n_cov_coef= ncol(Xo);data.n1= nrow(Xo)
interacmat <- function(Xo,u){Xo*u}
resinteracmat<-apply(as.matrix(Xo),2,interacmat,u=Xo)
IM<-matrix(resinteracmat,data.n1,n_cov_coef^2)
risskseth0 <- function(timeo,etime) ifelse(timeo[1]<etime,0,1)
RIh00 <- apply(as.array(timeo),1,risskseth0,etime=etime)
RIh00 <-t(as.matrix(RIh00))
nev_nozero<-n_eve;nev_nozero[nev_nozero>0]<-1
RIh0=t((c(nev_nozero)*t(RIh00)))
id_ind_zero <- which(apply(RIh0==0, 2, all))
if(length(id_ind_zero)>0){RIh0<- RIh0[,-id_ind_zero]}
Da=sum((1/a^2)*log(1+a*(H*EXB))-(1/a+di)*((H*EXB)/(1+a*(H*EXB))))
D2a= sum(-(2/a^3)*log(1+a*(H*EXB))+(2/a^2)*((H*EXB)/(1+a*(H*EXB)))+(1/a+di)*((H*EXB)^2/(1+a*(H*EXB))^2))
Db0=c(-((1/a+di)/(1+a*(H*EXB))))*(c(a*H*EXB)*Xo)+Xo
Db=colSums(Db0)
nev_nozero<-n_eve;nev_nozero<-nev_nozero[nev_nozero>0]
SW1=c((sqrt((1/a+di))/(1+a*(H*EXB)))*(a*H*EXB))*Xo
D2b=(t(SW1)%*%SW1)-matrix(c(colSums(c(((1/a+di)*(a*H*EXB))/(1+a*(H*EXB)))*IM)),n_cov_coef,n_cov_coef)
RI1=c(c(sqrt((1/a+di)))*c((EXB)^2/(1/a+(H*EXB))^2))*RIh0
D2h=(t(RI1)%*%(c(sqrt((1/a+di)))*RIh0))
diag(D2h)<-diag(D2h)-nev_nozero/nonzero_h0^2
DhDa=colSums(c((1/a^2)*((EXB/(1/a+(H*EXB)))-(1/a+di)*(EXB/(1/a+(H*EXB))^2)))*RIh0)
DbDa=colSums(c((1/a^2)*((H*EXB)/(1/a+(H*EXB)))-(1/a+di)*((1/a^2)*(H*EXB)/(1/a+(H*EXB))^2))*Xo)
DbDh=t(t(c(-(1/a+di)*((a*EXB)/(1+a*(H*EXB))-(a*EXB)*(a*H*EXB)/(1+a*(H*EXB))^2))*Xo)%*%RIh0)
MA=rbind(D2b,DbDh,DbDa);MB=rbind(t(DbDh),D2h,DhDa);MC=c(c(DbDa),c(DhDa),c(D2a))
HES=-cbind(MA,MB,MC)
INVE<-solve(HES)
seofthet=sqrt(diag(INVE))
SEofbetandsig=c(c(seofthet[1:n_cov_coef]),seofthet[length(seofthet)])
SEofbetandsig
}


.Llikgammacv = function (theta, cen1,cen2,H1,H2){
a = abs(theta[1]);R = abs(theta[2])
LSB=log(1+a*(H1+H2));LS1=log(1+a*H1);LS2=log(1+a*H2)
SB1=(1+a*(H1+H2))^(-1);SB2=(1+a*(H1+H2))^(-2)
S1=(1+a*H1)^(-1);S2=(1+a*H2)^(-1)
P11=(R+a)*R*SB2+(1-R)*R*SB1*(S1+S2)+((1-R)^2)*S1*S2
P10=R*SB1+(1-R)*S1
P01=R*SB1+(1-R)*S2
logli=sum((-((1-R)/a)*(LS1+LS2)-(R/a)*LSB)+cen1*cen2*log(P11)+
cen1*(1-cen2)*log(P10)+(1-cen1)*cen2*log(P01))
-logli
    }

.fdLlikgammacv =function (theta, cen1,cen2,H1,H2){
a = abs(theta[1]);R = abs(theta[2])
da=sum(((1 - R)/a^2 * (log(1 + a * H1) + log(1 + a * H2)) - ((1 - R)/a) * (H1/(1 + a * H1) + H2/(1 + a * H2)) - ((R/a) * ((H1 + H2)/(1 +
a * (H1 + H2))) - R/a^2 * log(1 + a * (H1 + H2))) + cen1 * cen2 * ((R * (1 + a * (H1 + H2))^(-2) + (R + a) * R * ((1 +
a * (H1 + H2))^((-2) - 1) * ((-2) * (H1 + H2))) + ((1 - R) * R * ((1 + a * (H1 + H2))^((-1) - 1) * ((-1) * (H1 + H2))) *
((1 + a * H1)^(-1) + (1 + a * H2)^(-1)) + (1 - R) * R * (1 + a * (H1 + H2))^(-1) * ((1 + a * H1)^((-1) - 1) * ((-1) *
H1) + (1 + a * H2)^((-1) - 1) * ((-1) * H2))) + (((1 - R)^2) * ((1 + a * H1)^((-1) - 1) * ((-1) * H1)) * (1 + a * H2)^(-1) +
((1 - R)^2) * (1 + a * H1)^(-1) * ((1 + a * H2)^((-1) - 1) * ((-1) * H2))))/((R + a) * R * (1 + a * (H1 + H2))^(-2) +
(1 - R) * R * (1 + a * (H1 + H2))^(-1) * ((1 + a * H1)^(-1) + (1 + a * H2)^(-1)) + ((1 - R)^2) * (1 + a * H1)^(-1) *
(1 + a * H2)^(-1))) + cen1 * (1 - cen2) * ((R * ((1 + a * (H1 + H2))^((-1) - 1) * ((-1) * (H1 + H2))) + (1 - R) * ((1 +
a * H1)^((-1) - 1) * ((-1) * H1)))/(R * (1 + a * (H1 + H2))^(-1) + (1 - R) * (1 + a * H1)^(-1))) + (1 - cen1) * cen2 * ((R *
((1 + a * (H1 + H2))^((-1) - 1) * ((-1) * (H1 + H2))) + (1 - R) * ((1 + a * H2)^((-1) - 1) * ((-1) * H2)))/(R * (1 + a *
(H1 + H2))^(-1) + (1 - R) * (1 + a * H2)^(-1)))))
dR=sum((1/a * (log(1 + a * H1) + log(1 + a * H2)) - 1/a * log(1 + a * (H1 + H2)) + cen1 * cen2 * (((R + (R + a)) * (1 + a * (H1 +
H2))^(-2) + ((1 - R) - R) * (1 + a * (H1 + H2))^(-1) * ((1 + a * H1)^(-1) + (1 + a * H2)^(-1)) - 2 * (1 - R) * (1 + a *
H1)^(-1) * (1 + a * H2)^(-1))/((R + a) * R * (1 + a * (H1 + H2))^(-2) + (1 - R) * R * (1 + a * (H1 + H2))^(-1) * ((1 +
a * H1)^(-1) + (1 + a * H2)^(-1)) + ((1 - R)^2) * (1 + a * H1)^(-1) * (1 + a * H2)^(-1))) + cen1 * (1 - cen2) * (((1 +
a * (H1 + H2))^(-1) - (1 + a * H1)^(-1))/(R * (1 + a * (H1 + H2))^(-1) + (1 - R) * (1 + a * H1)^(-1))) + (1 - cen1) *
cen2 * (((1 + a * (H1 + H2))^(-1) - (1 + a * H2)^(-1))/(R * (1 + a * (H1 + H2))^(-1) + (1 - R) * (1 + a * H2)^(-1)))))
gd=c(da,dR)
-gd
}


D2THTT=function (theta, cen1,cen2,H1,H2){
a = abs(theta[1]);R = abs(theta[2])
LSB=log(1+a*(H1+H2));dLSB=(H1+H2)/(1+a*(H1+H2));d2LSB=-(H1+H2)^2/(1+a*(H1+H2))^2
LS1=log(1+a*H1);dLS1=H1/(1+a*H1);d2LS1=-H1^2/(1+a*H1)^2
LS2=log(1+a*H2);dLS2=H2/(1+a*H2);d2LS2=-H2^2/(1+a*H2)^2
SB1=(1+a*(H1+H2))^(-1);dSB1=-((H1+H2)/(1+a*(H1+H2))^2);d2SB1=2*((H1+H2)^2/(1+a*(H1+H2))^3)
SB2=(1+a*(H1+H2))^(-2);dSB2=-2*((H1+H2)/(1+a*(H1+H2))^3);d2SB2=6*((H1+H2)^2/(1+a*(H1+H2))^4)
S1=(1+a*H1)^(-1);dS1=-(H1/(1+a*H1)^2);d2S1=2*(H1^2/(1+a*H1)^3)
S2=(1+a*H2)^(-1);dS2=-(H2/(1+a*H2)^2);d2S2=2*(H2^2/(1+a*H2)^3)
P11=(R+a)*R*SB2+(1-R)*R*SB1*(S1+S2)+((1-R)^2)*S1*S2
daP11=(R*SB2+(R^2+a*R)*dSB2+(R-R^2)*(dSB1*(S1+S2)+SB1*(dS1+dS2))+((1-R)^2)*(dS1*S2+S1*dS2))
d2aP11=(2*R*dSB2+(R^2+a*R)*d2SB2+
(R-R^2)*(d2SB1*(S1+S2)+2*dSB1*(dS1+dS2)+SB1*(d2S1+d2S2))+
((1-R)^2)*(d2S1*S2+2*dS1*dS2+S1*d2S2))
dRP11=((2*R+a)*SB2+(1-2*R)*SB1*(S1+S2)-2*(1-R)*S1*S2)
d2RP11=(2*SB2-2*SB1*(S1+S2)+2*S1*S2)
daRP11=(SB2+(2*R+a)*dSB2+(1-2*R)*(dSB1*(S1+S2)+SB1*(dS1+dS2))-2*(1-R)*(dS1*S2+S1*dS2))
P10=R*SB1+(1-R)*S1;daP10=R*dSB1+(1-R)*dS1;d2aP10=R*d2SB1+(1-R)*d2S1
P01=R*SB1+(1-R)*S2;daP01=R*dSB1+(1-R)*dS2;d2aP01=R*d2SB1+(1-R)*d2S2
dRP10=SB1-S1;dRP01=SB1-S2
daRP10=dSB1-dS1;daRP01=dSB1-dS2
d2a=sum((-2*((1-R)/a^3)*(LS1+LS2)+2*((1-R)/a^2)*(dLS1+dLS2)-((1-R)/a)*(d2LS1+d2LS2)-2*(R/a^3)*LSB+
2*(R/a^2)*dLSB-(R/a)*d2LSB)+cen1*cen2*((d2aP11/P11)-(daP11^2/P11^2))+
cen1*(1-cen2)*((d2aP10/P10)-(daP10^2/P10^2))+(1-cen1)*cen2*((d2aP01/P01)-(daP01^2/P01^2)))
d2R=sum(cen1*cen2*((d2RP11/P11)-(dRP11^2/P11^2))+cen1*(1-cen2)*(-(dRP10^2/P10^2))+
(1-cen1)*cen2*(-(dRP01^2/P01^2)))
daR=sum(((-1/a^2)*(LS1+LS2)+(1/a)*(dLS1+dLS2)+(1/a^2)*LSB-(1/a)*dLSB)+
cen1*cen2*((daRP11/P11)-((daP11*dRP11)/P11^2))+cen1*(1-cen2)*((daRP10/P10)-((daP10*dRP10)/P10^2))+
(1-cen1)*cen2*((daRP01/P01)-((daP01*dRP01)/P01^2)))
hesan=matrix(c(d2a,daR,daR,d2R),2,2)
-hesan
}

.bcgamfitcv<-function(X,Y,initfrailp,control){
time=Y[, 2];censor=Y[, 3]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)}
data.n1 <- nrow(X);data.n <-data.n1/2#### data.n is the number of pairs
risskset <- function(uniq_tim,x) ifelse(uniq_tim[1]<=x,1,0)
RI <- apply(as.array(uniq_tim),1,risskset,x=time)
RI <-t(RI)
numbevent <- function(uniq_tim,time,censor){
indic0=match(time,uniq_tim[1])*array(1:data.n1)
sum(censor[indic0[!is.na(indic0)]])}
n_eve<-apply(as.array(uniq_tim),1,numbevent,time,censor)
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
cph0 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = NULL, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph0$coefficients
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(x_bet)))
cen1=censor[indic1];cen2=censor[indic2]
if(length(initfrailp)==0){newtht<-c(0.5,0.5)}
if(length(initfrailp)>0){newtht<-initfrailp}
estim0=c(newtht)
W<-NULL
new.diff=1
e_frail<-.Expefrail(newtht=newtht,cen1=cen1,cen2=cen2,H1=H_bet_x[indic1],H2=H_bet_x[indic2])
W[indic1]<-e_frail$frail1;W[indic2]<-e_frail$frail2
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset =W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph1$coefficients
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet+W)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(x_bet)))
iter=0
repeat{
iter=iter+1
e_frail<-.Expefrail(newtht=newtht,cen1=cen1,cen2=cen2,H1=H_bet_x[indic1],H2=H_bet_x[indic2])
W[indic1]<-e_frail$frail1;W[indic2]<-e_frail$frail2
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset =W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph1$coefficients
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet+W)))
H0<-cumsum(c(n_eve/svexp_bet_xo))###obtain initial parameter for H0(t)
H_bet_x=c(H0[ind.haz]*exp(c(x_bet)))
fittr=do.call(nlminb, args=c(list(start=c(0.1,0.1), objective=.Llikgammacv,
gradient = .fdLlikgammacv,hessian =D2THTT,cen1=cen1,cen2=cen2,H1=H_bet_x[indic1],
H2=H_bet_x[indic2],lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par ##### obtain new estimates of sigma^2 and row
new.diff=max(abs(abs(c(newtht))-abs(estim0)))
estim0=c(newtht)
if((new.diff < control$tol)  |  (iter >= control$max.iter)) break}
if (iter > control$max.iter){warning("Ran out of iterations and did not converge")}
h0=diff(c(0,H0));nonzero_h0=h0[h0>0]
H=H0[ind.haz];h=h0[ind.haz];exp_bet_x=exp(c(x_bet))
lik=(sum(censor*c(x_bet))+sum(log(nonzero_h0))-fittr$objective+sum(censor))
DFF=try(D2THTT(theta=newtht,cen1,cen2,H1=H_bet_x[indic1],H2=H_bet_x[indic2])
,silent=TRUE)
if(length(attr(DFF,"class"))==0){
vcovth<-try(solve(DFF),silent=TRUE)
if(length(attr(vcovth,"class"))==0){
if(any(is.na(sqrt(diag(vcovth))))){
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}
colnames(vcovth) <- rownames(vcovth) <- c("theta","Row")}
if(length(attr(vcovth,"class"))!=0){
vcovth<-matrix(NA,2,2)
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")
colnames(vcovth) <- rownames(vcovth) <- c("theta","Row")}}
if(length(attr(DFF,"class"))!=0){
vcovth<-matrix(NA,2,2)
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")
colnames(vcovth) <- rownames(vcovth) <- c("theta","Row")}
vcov=cph1$var
colnames(vcov) <- rownames(vcov) <- colnames(X)
if(control$obt.se){
adj_se=try(.SEgammacv(sig=newtht[1],n_eve=n_eve,etime=uniq_tim,nonzero_h0=nonzero_h0,
exp_bet_xo=exp_bet_x[indx],Ho=H[indx],ho=h[indx],
censoro=censor[indx],timeo=time[indx],Xo=X[indx,]),silent=TRUE)
if(length(attr(adj_se,"class"))==0){adjse=c(adj_se,sqrt(diag(vcovth))[2])
if(any(is.na(adjse))){
adjse=c(sqrt(diag(vcov)),sqrt(diag(vcovth)))
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}}
if(length(attr(adj_se,"class"))!=0){
adjse=c(sqrt(diag(vcov)),sqrt(diag(vcovth)))
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}}
if(!control$obt.se){
adjse=c(sqrt(diag(vcov)),sqrt(diag(vcovth)))}
res <-list(coefficients=bet,frailparest= c(theta=newtht[1],Row=newtht[2]),
vcov2 = vcov,vcovth2= vcovth,loglilk0=cph0$loglik[1],loglilk=cph1$loglik[2],
Iloglilk=lik,cbasehaz=cbind(haz=H0,time=uniq_tim),X=X,time=time,censor=censor,
resid=cph1$residuals,lin.prid=cph1$linear.predictors,
frail=exp(W),stderr=adjse,iteration=iter,indx=indx,e.time=uniq_tim,n.event=n_eve,
converg = ifelse(new.diff< control$tol,TRUE,FALSE))
res$call <- match.call()
class(res) <- c(".bcgamfitcv")
res
}




.hesfunc = function(par,bet,censoro,Xo,RI,INC){
n=nrow(Xo)
i1=INC[,1];i2=INC[,2];i3=INC[,3]
W=par
W1<-W[i1]
W2<-W[i2]
x_bet<-c(Xo%*%bet+W)
expx_bet<-exp(c(Xo%*%bet+W))
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
as0=cumsum(c(censoro/svexp_bet_xo))
as1=((cumsum(c(censoro/(svexp_bet_xo^2))))*c(expx_bet^2)-(as0*c(expx_bet)))
as2=cumsum(c(censoro/svexp_bet_xo^2))
as3=as2[i3]*(expx_bet[i1])*(expx_bet[i2])
xx<-matrix(c(W1,W2,as1[i1],as1[i2],as3),(n/2),5)
xx
}

.llpenn= function (par,newtht,censoro,sumx,Xo,RI,INC){
n4=ncol(Xo)
n=nrow(Xo)
i1=INC[,1];i2=INC[,2];i3=INC[,3]
a=newtht[1];roh=newtht[2]
bet=par[array(1:n4)]
W=par[(1+n4):(n+n4)]
x_bet<-c(Xo%*%bet+W)
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
.llpenn=(-sum(censoro*(x_bet-log(svexp_bet_xo)))+
(1/(2*a*(1-roh^2)))*sum(W[i1]^2+W[i2]^2-2*roh*W[i1]*W[i2]))
Dw<-(-censoro+exp(c(x_bet))*cumsum(censoro/svexp_bet_xo))
Dw[i1]<-Dw[i1]+(W[i1]-roh*W[i2])*(1/(a*(1-roh^2)))
Dw[i2]<-Dw[i2]+(W[i2]-roh*W[i1])*(1/(a*(1-roh^2)))
attr(.llpenn,"gradient")<-c(as.vector(-sumx+colSums( c(censoro/svexp_bet_xo)*(RI%*%(Xo*c(exp(x_bet)))))),c(Dw))
.llpenn
}

.llppl= function (par,censoro,Xo,RI){
n4=ncol(Xo)
n=nrow(Xo)
bet=par[array(1:n4)]
W=par[(1+n4):(n+n4)]
x_bet<-c(Xo%*%bet+W)
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
logl=sum(censoro*(x_bet-log(svexp_bet_xo)))
logl
}


.SEpenal= function(par,newtht,censoro,Xo,RI,INC){
n4=ncol(Xo)
n=nrow(Xo)
i1=INC[,1];i2=INC[,2];i3=INC[,3]
a=newtht[1];roh=newtht[2]
bet=par[1:n4]
W=par[(1+n4):(n+n4)]
x_bet<-c(Xo%*%bet+W)
vexpx_bet<-c(exp(x_bet))
svexp_bet_xo=as.vector(RI%*%(exp(x_bet)))
dbexp_bet_xo=(RI%*%(vexpx_bet*Xo))
interacmat <- function(Xo,u){Xo*u}
resinteracmat<-apply(as.matrix(Xo),2,interacmat,u=Xo)
Xoxo<-matrix(resinteracmat,n,n4^2)
as0=cumsum(c(censoro/svexp_bet_xo))*vexpx_bet
ma0<-matrix(0,n,n)
ma0[lower.tri(ma0, diag = TRUE)]<-1
ma1<-t((ma0*vexpx_bet))
ma2<-ma1*(-cumsum(c(censoro/(svexp_bet_xo^2))))*vexpx_bet
digg<-diag(ma2)+as0
ma3<-t(ma2)
D2W<-ma3+ma2
digg[i1]<-digg[i1]+(1/(a*(1-roh^2)))
digg[i2]<-digg[i2]+(1/(a*(1-roh^2)))
diag(D2W)<-digg
D2W[cbind(i1,i2)]<-D2W[cbind(i1,i2)]-((roh)/(a*(1-roh^2)))
D2W[cbind(i2,i1)]<-D2W[cbind(i2,i1)]-((roh)/(a*(1-roh^2)))
dwbpart20<-(dbexp_bet_xo)*(censoro/(svexp_bet_xo^2))
cumsummat <- function(dwbpart20){cumsum(dwbpart20)}
dwbpart21<-apply(as.matrix(dwbpart20),2,cumsummat)
DWB<-(as0*Xo-dwbpart21*(vexpx_bet))
d2bpart1=matrix(c(colSums((censoro/svexp_bet_xo)*(RI%*%(vexpx_bet*Xoxo)))),n4,n4)
d2bpart2=t((censoro/svexp_bet_xo)*dbexp_bet_xo)%*%(dbexp_bet_xo*(censoro/svexp_bet_xo))
D2B<-d2bpart1-d2bpart2
MA=rbind(D2B,DWB);MB=rbind(t(DWB),D2W)
HES=cbind(MA,MB)
INVE<-solve(HES)
se=sqrt(diag(INVE))
se[1:n4]
}



.Lliklogncv= function (thetaa, xx){
xx=as.matrix(xx)
W1=xx[,1];W2=xx[,2];S1=xx[,3];S2=xx[,4];S3=xx[,5]
a = abs(thetaa[1]) #var
e = thetaa[2] #corr
n = nrow(xx)
P11<-(S1*S2-(1/(a*(1-e^2)))*(S1+S2)+
(1/(a^2*(1-e^2)^2))-S3^2-2*(e/(a*(1-e^2)))*S3-(e^2/(a^2*(1-e^2)^2)))
logli=(-(1/(2*a*(1-e^2)))*sum((W1^2+W2^2-2*e*W1*W2))-n*log(a)-(n/2)*log(1-e^2)+
(-1/2)*sum(log(P11)))
-logli}

.fdLliklogncv = function (thetaa, xx){
xx=as.matrix(xx)
W1=xx[,1];W2=xx[,2];S1=xx[,3];S2=xx[,4];S3=xx[,5]
a = abs(thetaa[1]) #var
e = thetaa[2] #corr
n=nrow(xx)
P11<-(S1*S2-(1/(a*(1-e^2)))*(S1+S2)+
(1/(a^2*(1-e^2)^2))-S3^2-2*(e/(a*(1-e^2)))*S3-(e^2/(a^2*(1-e^2)^2)))
daP11<-( (1/(a^2*(1-e^2)))*(S1+S2)-(2/(a^3*(1-e^2)^2))+2*(e/(a^2*(1-e^2)))*S3+
(2*e^2/(a^3*(1-e^2)^2)))
deP11<- ( -(2*a*e/(a*(1-e^2))^2)*(S1+S2)+((4*a^2*(1-e^2)*e)/(a^2*(1-e^2)^2)^2)+
2*(-1/(a*(1-e^2))-2*(e^2*a)/(a*(1-e^2))^2)*S3-
(e*2/(a^2*(1-e^2)^2)+(4*a^2*(1-e^2)*e^3)/(a^2*(1-e^2)^2)^2))
grr=c( ((1/(2*a^2*(1-e^2)))*sum((W1^2+W2^2-2*e*W1*W2))-n/a-(1/2)*sum(daP11/P11)),
(-((4*a*e)/(2*a*(1-e^2))^2)*sum((W1^2+W2^2-2*e*W1*W2))+(1/(a*(1-e^2)))*sum(W1*W2)+
n*(e/(1-e^2))-(1/2)*sum(deP11/P11)))
-grr
}



.sdLliklogncv = function (thetaa, xx){
xx=as.matrix(xx)
W1=xx[,1];W2=xx[,2];S1=xx[,3];S2=xx[,4];S3=xx[,5]
a = abs(thetaa[1]) #var
e = thetaa[2] #corr
n=nrow(xx)
P11<-(S1*S2-(1/(a*(1-e^2)))*(S1+S2)+
(1/(a^2*(1-e^2)^2))-S3^2-2*(e/(a*(1-e^2)))*S3-(e^2/(a^2*(1-e^2)^2)))
daP11<-( (1/(a^2*(1-e^2)))*(S1+S2)-(2/(a^3*(1-e^2)^2))+2*(e/(a^2*(1-e^2)))*S3+
(2*e^2/(a^3*(1-e^2)^2)))
deP11<- ( -(2*a*e/(a*(1-e^2))^2)*(S1+S2)+((4*a^2*(1-e^2)*e)/(a^2*(1-e^2)^2)^2)+
2*(-1/(a*(1-e^2))-2*(e^2*a)/(a*(1-e^2))^2)*S3-
(e*2/(a^2*(1-e^2)^2)+(4*a^2*(1-e^2)*e^3)/(a^2*(1-e^2)^2)^2))
d2aP11<-( -(2/(a^3*(1-e^2)))*(S1+S2)+(6/(a^4*(1-e^2)^2))-4*(e/(a^3*(1-e^2)))*S3-
(6*e^2/(a^4*(1-e^2)^2)))
daRP11<-( ((2*a^2*e)/(a^2*(1-e^2))^2)*(S1+S2)-(8*(a^3*(1-e^2)*e)/(a^3*(1-e^2)^2)^2)+
2*(1/(a^2*(1-e^2))+2*(a^2*e^2)/(a^2*(1-e^2))^2)*S3+
(4*e/(a^3*(1-e^2)^2)+8*(a^3*(1-e^2)*e^3)/(a^3*(1-e^2)^2)^2))
deP11<- ( -(2*a*e/(a*(1-e^2))^2)*(S1+S2)+((4*a^2*(1-e^2)*e)/(a^2*(1-e^2)^2)^2)+
2*(-1/(a*(1-e^2))-2*(e^2*a)/(a*(1-e^2))^2)*S3-
( e*2/(a^2*(1-e^2)^2)+(4*a^2*(1-e^2)*e^3)/(a^2*(1-e^2)^2)^2))
d2eP11<- ( -(2*a/(a*(1-e^2))^2 + 8*a^2*e^2*(a*(1-e^2))/(a*(1-e^2))^4)*(S1+S2)+
((4 * a^2 * (1 - e^2) - 4 * a^2 * (2 * e) * e)/(a^2 * (1 - e^2)^2)^2 +
(4 * a^2 * (1 - e^2) * e) * (2 * (a^2 * (2 * (2 * e * (1 - e^2))) * (a^2 * (1 - e^2)^2)))/((a^2 * (1 - e^2)^2)^2)^2) +
(-(2 * (a * (2 * e)/(a * (1 - e^2))^2 + (2 * (2 * e * a)/(a * (1 - e^2))^2 + 2 * (e^2 * a) * (2 * (a * (2 * e) * (a * (1 -
    e^2))))/((a * (1 - e^2))^2)^2))))*S3-
(2/(a^2 * (1 - e^2)^2) + e * 2 * (a^2 * (2 * (2 * e * (1 - e^2))))/(a^2 *
    (1 - e^2)^2)^2 + ((4 * a^2 * (1 - e^2) * (3 * e^2) - 4 *
    a^2 * (2 * e) * e^3)/(a^2 * (1 - e^2)^2)^2 + (4 * a^2 * (1 -
    e^2) * e^3) * (2 * (a^2 * (2 * (2 * e * (1 - e^2))) * (a^2 *
    (1 - e^2)^2)))/((a^2 * (1 - e^2)^2)^2)^2)))
D2a=((-(2 * (2 * a) * (1 - e^2)/(2 * a^2 * (1 - e^2))^2))*sum((W1^2+W2^2-2*e*W1*W2))+
n/a^2-(1/2)*sum(d2aP11/P11-daP11^2/P11^2))
DaR=((2 * a^2 * (2 * e)/(2 * a^2 * (1 - e^2))^2)*sum((W1^2+W2^2-2*e*W1*W2))-
(1/(a^2*(1-e^2)))*sum(W1*W2)-(1/2)*sum(daRP11/P11-(daP11*deP11)/P11^2))
D2R=( (-(4 * a/(2 * a * (1 - e^2))^2 + (4 * a * e) * (2 * (2 * a * (2 *
    e) * (2 * a * (1 - e^2))))/((2 * a * (1 - e^2))^2)^2))*sum((W1^2+W2^2-2*e*W1*W2))+
((2*a*e)/(a*(1-e^2))^2)*sum(W1*W2)+(a * (2 * e)/(a * (1 - e^2))^2)*sum(W1*W2)+
n*( 1/(1-e^2)+2*e^2/(1-e^2)^2)-(1/2)*sum(d2eP11/P11-deP11^2/P11^2))
hesi=matrix(c(D2a,DaR,DaR,D2R),2,2)
-hesi
}





.machinelern = function (newtht, newtht0,direct,lo,G){
curdirect=estimdiff=c(newtht-newtht0)
curdirect[curdirect<0]<-c(-1);curdirect[curdirect==0]<-0;curdirect[curdirect>0]<-1
prevdirect=direct
newtht0<-newtht
if(G==0){
for(j in 1:length(newtht)){
if(prevdirect[j]==1){
if(curdirect[j]==1){
lo[j]=lo[j]+2;newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}
if(curdirect[j]!=1){lo[j]=0.4*lo[j];newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}}
if(prevdirect[j]==-1){
if(curdirect[j]==-1){
lo[j]=lo[j]+2;newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}
if(curdirect[j]!=-1){lo[j]=0.4*lo[j];newtht0[j]<-newtht[j]+lo[j]*estimdiff[j]}}
if(prevdirect[j]==0){newtht0[j]<-newtht[j];lo[j]=0.5*lo[j]}}
if(any(newtht0<0.00000001)){
newtht0[newtht0<0.00000001]<-newtht[newtht0<0.00000001]
curdirect[newtht0<0.00000001]<-0
lo[newtht0<0.00000001]<-0}
if(newtht0[length(newtht0)]>0.9999999){
newtht0[length(newtht0)]<-newtht[length(newtht0)]
curdirect[length(newtht0)]<-0;lo[length(newtht0)]<-0}}
if(G==2){newtht0<-newtht}
list(newt=newtht0,direct=curdirect,lo=lo)}





.bclognfitcv<-function(X,Y,initfrailp,control){
time=Y[, 2];censor=Y[, 3]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)
if(any(is.na(ind.haz))){
tord_diff<-as.array(diff(c(0,timeo)))
id.zero_tord <- which(apply(((tord_diff<0.0000001)&(tord_diff>0)),1, all))
if(length(id.zero_tord)>0){time[indx[id.zero_tord]]<- time[indx[id.zero_tord-1]]
order=sort(time, decreasing = FALSE, index.return = TRUE);indx=order$ix;timeo=order$x}
uniq_tim<-unique(sort(time))
ind.haz=match(time,uniq_tim)}
data.n1 <- nrow(X);data.n <-data.n1/2#### data.n is the number of pairs
indic1<-2*array(1:data.n)-1;indic2<-2*array(1:data.n)
PPID<-array(1:data.n1)
pid=match(PPID,indx);i1<-pid[indic1];i2<-pid[indic2]
i3 <-pmin(i1,i2)
INC<-matrix(c(i1,i2,i3),data.n,3)
risskset <- function(timeo,x) ifelse(timeo[1]>=x,1,0)
RI <- apply(as.array(timeo),1,risskset,x=timeo)
numbevent <- function(uniq_tim,time,censor){
indic0=match(time,uniq_tim[1])*array(1:data.n1)
sum(censor[indic0[!is.na(indic0)]])}
n_eve<-apply(as.array(uniq_tim),1,numbevent,time,censor)
Xo=X[indx,]
Xo=as.matrix(Xo)
censoro=censor[indx]
cph0 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = NULL, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-cph0$coefficients
ncovar_coef=length(bet)
W<-par<-rep(0,(data.n1+ncovar_coef))
par[1:ncovar_coef]<-bet
sumx=as.vector(colSums(censoro*Xo))
newtht0=newtht<-NULL
if(length(initfrailp)==0){
newtht0<-newtht<-c(0.5,0.5)}
if(length(initfrailp)>0){
newtht0<-newtht<-initfrailp
newthtinit=newtht0}
fittrpen=do.call(nlm, args=c(list(f=.llpenn,p=par,newtht=newtht0,
censoro=censoro,sumx=sumx,Xo=Xo,RI=RI,INC=INC,fscale = control$fscale,
print.level = control$print.level, ndigit = control$ndigit,steptol = control$steptol,
iterlim= control$iterlim,gradtol = control$gradtol,check.analyticals = control$check.analyticals)))
par=fittrpen$estimate
xx<-.hesfunc(par=par[(1+ncovar_coef):(data.n1+ncovar_coef)],bet=par[1:ncovar_coef],
censoro=censoro,Xo=Xo,RI=RI,INC=INC)
fittr=do.call(nlminb, args=c(list(start=newtht0, objective=.Lliklogncv,
gradient = .fdLliklogncv,hessian =.sdLliklogncv,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par##### obtain new estimates of sigma^2 and row
old.diff=max(abs(newtht0-newtht))
direct<-newtht-newtht0
direct[direct<0]<--1
direct[direct>0]<-1
direct[direct==0]<-0
oldlik<-fittr$objective
estim0=newtht0<-c(newtht)
new.diff=1
lo=rep(5,length(newtht))
G=0
if(control$fastfit){
iter=0
repeat{
iter=iter+1
fittrpen=do.call(nlm, args=c(list(f=.llpenn,p=par,newtht=newtht0,
censoro=censoro,sumx=sumx,Xo=Xo,RI=RI,INC=INC,fscale = control$fscale,
print.level = control$print.level, ndigit = control$ndigit,steptol = control$steptol,
iterlim= control$iterlim,gradtol = control$gradtol,check.analyticals = control$check.analyticals)))
par=fittrpen$estimate
xx<-.hesfunc(par=par[(1+ncovar_coef):(data.n1+ncovar_coef)],bet=par[1:ncovar_coef],
censoro=censoro,Xo=Xo,RI=RI,INC=INC)
fittr=do.call(nlminb, args=c(list(start=newtht0, objective=.Lliklogncv,
gradient = .fdLliklogncv,hessian =.sdLliklogncv,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par##### obtain new estimates of sigma^2 and row
new.diff=abs(oldlik-(fittr$objective))
oldlik<-fittr$objective
newdif=max(abs(newtht-estim0))
if(iter>200){G=2}
newtht00<-.machinelern(newtht,newtht0,direct,lo,G)
newtht0<-newtht00$newt;direct<-newtht00$direct;lo<-newtht00$lo
estim0=c(newtht0)
if((newdif< control$toll)  |  (iter >=control$max.iter2)) break}}
if(!control$fastfit){
new.diff=1
iter=0
repeat{
iter=iter+1
fittrpen=do.call(nlm, args=c(list(f=.llpenn,p=par,newtht=newtht,
censoro=censoro,sumx=sumx,Xo=Xo,RI=RI,INC=INC,fscale = control$fscale,
print.level = control$print.level, ndigit = control$ndigit,steptol = control$steptol,
iterlim= control$iterlim,gradtol = control$gradtol,check.analyticals = control$check.analyticals)))
par=fittrpen$estimate
xx<-.hesfunc(par=par[(1+ncovar_coef):(data.n1+ncovar_coef)],bet=par[1:ncovar_coef],
censoro=censoro,Xo=Xo,RI=RI,INC=INC)
fittr=do.call(nlminb, args=c(list(start=newtht0, objective=.Lliklogncv,
gradient = .fdLliklogncv,hessian =.sdLliklogncv,xx=xx,lower = control$lower, upper = control$upper,
control=control$nlminb_control)))
newtht=fittr$par##### obtain new estimates of sigma^2 and row
new.diff=abs(oldlik-(fittr$objective))
oldlik<-fittr$objective
if((new.diff < control$toll)  |  (iter >=control$max.iter2)) break}}
if (iter >= control$max.iter2){warning("Ran out of iterations and did not converge")}
W<-par[(1+ncovar_coef):(data.n1+ncovar_coef)]
W[indx]<-W
cph1 <- survival::agreg.fit(x = X, y = Y, strata = NULL,
offset = W, init = NULL, control = survival::coxph.control(),
weights = NULL, method = "breslow", rownames = NULL)
bet<-par[1:ncovar_coef]
llppl<-.llppl(par=par,censoro=censoro,Xo=Xo,RI=RI)
lik=(-fittr$objective+llppl)
fishinfmat<-.sdLliklogncv(thetaa=newtht, xx=xx)
vcovth<-solve(fishinfmat)
if(any(is.nan(sqrt(diag(vcovth))))){
warning("Check the approperiatness of the model.")
warning("frailty parameters might be at the boundary of parameter space.")
warning("In obtaining SE Na is produced")}
colnames(vcovth) <- rownames(vcovth) <- c("theta","Row")
vcov=cph1$var
colnames(vcov) <- rownames(vcov) <- colnames(X)
adj_se=.SEpenal(par=par,newtht=newtht,censoro=censoro,Xo=Xo,RI=RI,INC=INC)
adjse=c(adj_se,sqrt(diag(vcovth)))
rissksett <- function(uniq_tim,x) ifelse(uniq_tim[1]<=x,1,0)
RI <- apply(as.array(uniq_tim),1,rissksett,x=time)
RI <-t(RI)
x_bet<-X%*%bet
svexp_bet_xo=as.vector(RI%*%(exp(x_bet+W)))
H0<-cumsum(c(n_eve/svexp_bet_xo))
res <-list(coefficients=bet,frailparest= c(theta=newtht[1],Row=newtht[2]),
vcov2 = vcov,vcovth= vcovth,loglilk0=cph0$loglik[1],loglilk=cph1$loglik[2],
penloglilk=-fittrpen$minimum,
Iloglilk=lik,cbasehaz=cbind(haz=H0,time=uniq_tim),X=X,time=time,censor=censor,
resid=cph1$residuals,lin.prid=cph1$linear.predictors,
frail=exp(W),stderr=adjse,
iteration=iter,indx=indx,e.time=uniq_tim,n.event=n_eve,
converg = ifelse(new.diff<control$tol,TRUE,FALSE))
res$call <- match.call()
class(res) <- c(".bclognfitcv")
res
}




.gendata1<-function(p.size,c.rate,fraildistrn,frail.par,bhaz.arg,covar.arg){
n<-p.size; n1=n*2;IID=array(1:n1);indic1=2*array(1:n)-1;indic2=2*array(1:n)
PID=1;e1=array(1:n);PID[indic1]=e1;PID[indic2]=e1
if(fraildistrn==c("gamma")){
lam0=lam1=lam2=1/(frail.par[1])
k0=frail.par[2]/(frail.par[1]);k1=k2=(1-frail.par[2])/frail.par[1]
y0=rgamma(n,shape=k0,scale=1/lam0);y1=rgamma(n,shape=k1,scale=1/lam1)
y2=rgamma(n,shape=k2,scale=1/lam2)
z1=(y0+y1);z2=(y0+y2)}
if(fraildistrn==c("lognormal")){
k0=frail.par[2]*(frail.par[1]);k1=k2=frail.par[1]-k0
y0=rnorm(n,mean=0,sd=sqrt(k0));y1=rnorm(n,mean=0,sd=sqrt(k1))
y2=rnorm(n,mean=0,sd=sqrt(k2))
z1=exp((y0+y1)); z2=exp((y0+y2))}
if(pmatch(c(covar.arg$types),c("B","U"))==1){
x11<-rbinom(n,size=covar.arg$size,prob=covar.arg$prob)
x12<-rbinom(n,size=covar.arg$size,prob=covar.arg$prob)}
if(pmatch(c(covar.arg$types),c("B","U"))==2){
x11<-runif(n,min=covar.arg$min,max=covar.arg$max)
x12<-runif(n,min=covar.arg$min,max=covar.arg$max)}
u1<-runif(n,  min=0, max=1);u2<-runif(n,  min=0, max=1)
if (bhaz.arg$distrn==c("weibull")){
T1 <- (-log(u1) / ((bhaz.arg$scale)*z1*exp(c(x11*covar.arg$coefs))))^(1/(bhaz.arg$shape))
T2 <- (-log(u2) / ((bhaz.arg$scale)*z2*exp(c(x12*covar.arg$coefs))))^(1/(bhaz.arg$shape))}
if (bhaz.arg$distrn==c("gompertz")){
T1 <- 1/(bhaz.arg$shape)*log(1-(bhaz.arg$shape)*log(u1)/((bhaz.arg$scale)*exp(c(x11*covar.arg$coefs))*z1))
T2 <- 1/(bhaz.arg$shape)*log(1-(bhaz.arg$shape)*log(u2)/((bhaz.arg$scale)*exp(c(x12*covar.arg$coefs))*z2))}
if (bhaz.arg$distrn==c("exponential")){
T1 <- ((-log(u1)) / ((bhaz.arg$rate)*z1*exp(c(x11*covar.arg$coefs))))
T2 <- ((-log(u2)) / ((bhaz.arg$rate)*z2*exp(c(x12*covar.arg$coefs))))}
cen1=cen2<-NULL
if(c.rate==0){t1<-T1;t2<-T2;cen1=cen2<-rep(1,n)}
if(c.rate>0){
order=sort(T1, decreasing = FALSE, index.return = TRUE);indx1=order$ix
order=sort(T2, decreasing = FALSE, index.return = TRUE);indx2=order$ix
cr1=round(c(0.05*n),digits=0);I5p=array(1:cr1)
q1005<-quantile(T1, probs = c(0.05));q1099<-quantile(T1, probs = c(0.9999))
q2005<-quantile(T2, probs = c(0.05));q2099<-quantile(T2, probs = c(0.9999))
cen1=rbinom(n, size=1, prob=c(1-c.rate/0.95));cen1[indx1[I5p]]<-1
cen2=rbinom(n, size=1, prob=c(1-c.rate/0.95));cen2[indx2[I5p]]<-1
CC01<-runif(n,  min=q1005, max=q1099);CC1<-c(((q1005+min(CC01))/2),CC01)
CC02<-runif(n,  min=q2005, max=q2099);CC2<-c(((q2005+min(CC02))/2),CC02)
C<-c(CC1,CC2);t1<-T1;t2<-T2;CL=1
for(j in 1:n){if(cen1[j]==0){CL=C[C<T1[j]]; if(length(CL)<=1){t1[j]<-0.95*T1[j]}
if(length(CL)>1){t1[j]<-sample(CL, size=1)}}
if(cen2[j]==0){CL=C[C<T2[j]];if(length(CL)<=1){t2[j]<-0.95*T2[j]}
if(length(CL)>1){t2[j]<-sample(CL, size=1)}}}}
X1<-time<-censor<-NULL
IID0=array(1:n)
time[indic1]=t1;time[indic2]=t2;censor[indic1]=cen1;censor[indic2]=cen2
X1[indic1]=x11;X1[indic2]=x12
data<- data.frame(IID=IID,PID=PID,time=time,censor=censor,X1=X1)
data}







