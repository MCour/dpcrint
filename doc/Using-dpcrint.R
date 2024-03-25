## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=200,
  fig.width=4,
  fig.height=3
)

## ----setup--------------------------------------------------------------------
library(dpcrint)
par(mar=c(3,3,1,1))

## ----datasim0-----------------------------------------------------------------
N=26000 # number of partitions
lambda=0.05 # average number of sequences in one partition
J=5 # number of targets
pobs=c(0.8,0.6,0.75,0.99,0.98) # amplification probability for each target
coverseq=defseqtypes(J) # sequence types definition in terms of target coverage
L=nrow(coverseq)
propl=numeric(L) # proportion of each sequence type
propl[1]=1 
# print(cbind(coverseq,propl)) 
ressim=SimulObs(N,lambda,propl,pobs)
ObsTypes=ressim$ObsTypes
NY=ressim$NY # number of each observation type in terms of target coverage
print(NY)
# identical(coverseq,ressim$coverseq)

## ----optpobs------------------------------------------------------------------
respobs=DirectOptPobs(NY,ObsTypes,lambda0=lambda)
print(respobs)
print(pobs)

## ----optpobslam---------------------------------------------------------------
respobslam=DirectOptPobs(NY,ObsTypes,evallam=TRUE)
print(respobslam$pobs)
print(pobs)
print(respobslam$lambda)
print(lambda)

## ----optpropl0, results='hide',cache=T----------------------------------------
respropl=OptProplEM(NY,ObsTypes,pobs=pobs,lambda0=lambda)
abline(h=propl,col=1:L,lty=2)
plot(as.vector(respropl$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)

## ----optpropllam0, results='hide',cache=T-------------------------------------
respropllam=OptProplEM(NY,ObsTypes,pobs=pobs,evallam=T)
abline(h=propl,col=1:L,lty=2)
plot(as.vector(respropllam$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)
plot(respropllam$LAMBDA,ylab="lambda",xlab="iterations")
abline(h=lambda)

## ----datasim------------------------------------------------------------------
propl=exp(0.99)
for(j in 2:J){
  propl=c(propl,rep(exp(0.99*j),j))
}
propl=propl/sum(propl) # proportion of each sequence type
# cbind(coverseq,propl) 
ressim=SimulObs(N,lambda,propl,pobs)
ObsTypes=ressim$ObsTypes
NY=ressim$NY # number of each observation type in terms of target coverage
print(NY)
# identical(coverseq,ressim$coverseq)

## ----optpropl, results='hide',cache=T-----------------------------------------
respropl=OptProplEM(NY,ObsTypes,pobs=pobs,lambda0=lambda)
abline(h=propl,col=1:L,lty=2)
plot(as.vector(respropl$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)

## ----optpropllam, results='hide',cache=T--------------------------------------
respropllam=OptProplEM(NY,ObsTypes,pobs=pobs,evallam=T)
abline(h=propl,col=1:L,lty=2)
plot(as.vector(respropllam$propl)~propl,ylab="estimated proportions",xlab="true proportions",pch=19)
abline(0,1,lty=2)
plot(respropllam$LAMBDA,ylab="lambda",xlab="iterations")
abline(h=lambda)

