#' EM algorithm estimating proportions of each sequence type
#'
#' @import stats
#' @import graphics
#'
#' @param NY numerical vector, contains the number of each observation type.
#' @param ObsTypes binary matrix, contains, for each observation type (in row), if each target is present (1) of absent (0).
#' @param LatTypes (optional) matrix, contains all the possibilities for the number of each sequence type inside a partition.
#' @param pobs (optional) vector, amplification probability for each target
#' @param lambda0 (optional) numeric, the average number of sequences in one partition of ddPCR. If not provided, it will be estimated together with the amplification probabilities.
#' @param prec (optional) integer, the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#' @param propl0 (optional) vector, initial value for the proportions of each possible sequence type
#' @param Nit (optional) integer, maximal number of iteration in the EM algorithm
#' @param plotopt (optional) bolean, should a parameter convergence plot be drawn?
#' @param evallam (optional) bolean, should lambda0 be used as is of should it be estimated along with amplification probabilities? Is automatically switched to TRUE if no lambda0 is provided.
#' @param convdiag (optional) numeric, under which limit should the algorithm consider it has converged with respect to proportions learning?
#' @param convdiaglam (optional) numeric, under which limit should the algorithm consider it has converged with respect to lambda learning?
#'
#' @return a list of class ProplEM containing
#' the successive values of the learned sequence types frequencies in the EM algorithm ($PROPL),
#' the final learned sequence types frequencies ($propl),
#' the target covered by each of those sequence types ($coverseq),
#' the average number of sequences in a partition ($lambda),
#' the successive values of the average number of sequences in a partition in the EM algorithm ($LAMBDA)
#' @export
OptProplEM=function(NY,ObsTypes,LatTypes=NULL,pobs=NULL,lambda0=NULL,prec=100,
                    propl0=NULL,Nit=300,plotopt=T,evallam=F,convdiag=1e-7,convdiaglam=1e-3){


  #### INITIALIZATION ####
  J=ncol(ObsTypes)
  coverseq=defseqtypes(J) # defining the sequence types (we wish to recover their frequencies)
  L=nrow(coverseq)
  if(is.null(pobs)){
    pobs=rep(0.99,J)
  }
  if(is.null(propl0)){
    lpropl=exp(rnorm(L,0,3))
    propl=lpropl/sum(lpropl)
  }else{propl=propl0}
  if(is.null(lambda0)){
    evallam=T
    lambda0=0.5
  }else{
    if(lambda0>1){
      warning("The algorithm may give biased results with such high concentration (latent space too big)")}
    if(evallam){lambda=rgamma(1,prec*lambda0,prec)
    print(lambda)
    }else{lambda=lambda0}
  }
  N=sum(NY)
  if(is.null(LatTypes)){
    if(evallam){
      LatTypes=LatLibraryLam(Nit=floor(N/2)*2,NK=floor(N/2)*6,lambda0=lambda0*2,prec=prec,L=L)
    }else{
      LatTypes=LatLibrary(Nit=floor(N/2)*2,NK=floor(N/2)*6,lambda=lambda0*2,L=L)
    }
  }
  LatTypesj=LatTypes%*%coverseq
  wobs=1-pobs
  Na=nrow(ObsTypes)
  Kmaxl=apply(LatTypes,2,max)
  # IF POBS IS KNOWN:
  PKA=PKAcompute(LatTypesj,ObsTypes,wobs)


  #### EM Algorithm ####
  PROPL=matrix(data=NA,ncol=L,nrow=Nit)
  if(evallam){LAMBDA=numeric(Nit)*NA}else{LAMBDA=NULL}
  for(it in 1:Nit){
    PROPL[it,] = propl
    if(evallam){LAMBDA[it]=lambda}
    if((it%%10)==0){print(paste("iter=",it,"/",Nit))}
    lambdal=lambda*propl

    # A priori probability of each possible latent scheme (within LatTypes)
    # (given the Poisson parameter and the sequences proportions):
    PK=PKcompute(LatTypes,lambdal)

    # Probability of each latent scheme given each observation type
    PPKA=PK*PKA
    # !! up to a constant !! => renormalization
    PPKA=t(t(PPKA)/apply(PPKA,2,sum))

    # E Step:
    EK=EKcompute(LatTypes,Kmaxl,PPKA)

    # M step:
    KYEK=NY%*%EK
    normpl=as.numeric(t(NY)%*%rowSums(EK))
    proplnew=KYEK/normpl
    # print(propl)

    if(evallam){
      lambdanew=normpl/sum(NY)
      if( (sqrt(mean((propl-proplnew)^2))<convdiag) & (abs(lambda-lambdanew) < convdiaglam) ) break
    }else{
      if(sqrt(mean((propl-proplnew)^2))<convdiag) break
    }


    propl=proplnew
    if(evallam){lambda=lambdanew}

  }

  if(plotopt){
    matplot(PROPL,type="l")
    # library(scales)
    # abline(h=proplT,lty=1:L,col=alpha(1:L,0.3),lwd=2)
  }

  res=list(propl=propl,PROPL=PROPL,lambda=lambda,LAMBDA=LAMBDA,coverseq=coverseq)
  class(res)="ProplEM"

  return(res)
}



#' Latente types library building
#'
#' @param Nit (integer), number of sequence type proportion generation.
#' @param NK (integer), number of latent vector generation per generated proportion. Nk has to be a multiple of Nit.
#' @param lambda (numeric), average number of sequences in one partition of ddPCR.
#' @param L (interger) number of sequence types
#'
#' @return a matrix containing all the possibilities for the number of each sequence type inside a partition.
#' @export
LatLibrary=function(Nit=1000,NK=100000,lambda,L){
  Kilsamp=matrix(data=NA,ncol=L,nrow=NK)
  if(NK%%Nit!=0){stop("NK HAS to be a multiple of Nit")}
  for(it in 1:Nit){
    # lpropl=rnorm(L-1,0,3)
    # propl=c(1/(1+sum(exp(lpropl))),exp(lpropl)/(1+sum(exp(lpropl))))
    lpropl=exp(rnorm(L,0,3))
    propl=lpropl/sum(lpropl)
    # print(round(propl,1))
    lambdal=lambda*propl
    for(l in 1:L){
      Kilsamp[(NK/Nit)*(it-1)+(1:(NK/Nit)),l]=rpois(NK/Nit,lambdal[l])
    }
  }
  Kilobs=apply(Kilsamp,1,paste0,collapse="-")
  LatTypesChar=sort(unique(Kilobs))
  M=length(LatTypesChar) # number of possible values of the target vector
  LatTypes=matrix(data=NA, ncol=l,nrow=M)
  for(m in 1:M){
    LatTypes[m,]=Kilsamp[which(Kilobs==LatTypesChar[m])[1],]
  }
  return(LatTypes)
}

#' Latente types library building
#'
#' @param Nit (integer), number of sequence type proportion generation.
#' @param NK (integer), number of latent vector generation per generated proportion. Nk has to be a multiple of Nit.
#' @param lambda0 (numeric), a starting point for average number of sequences in one partition of ddPCR.
#' @param prec (integer), the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#' @param L (interger) number of sequence types.
#'
#' @return a matrix containing all the possibilities for the number of each sequence type inside a partition.
#' @export
LatLibraryLam=function(Nit=1000,NK=100000,lambda0=0.45,prec=100,L){
  Kilsamp=matrix(data=NA,ncol=L,nrow=NK)
  lambdasamp=rgamma(Nit,prec*lambda0,prec)
  if(NK%%Nit!=0){stop("NK HAS to be a multiple of Nit")}
  for(it in 1:Nit){
    # lpropl=rnorm(L-1,0,3)
    # propl=c(1/(1+sum(exp(lpropl))),exp(lpropl)/(1+sum(exp(lpropl))))
    lpropl=exp(rnorm(L,0,3))
    propl=lpropl/sum(lpropl)
    # print(round(propl,1))
    #lambda=exp(rnorm(1,-2.5,1.5))
    lambdal=lambdasamp[it]*propl
    # print(lambdal)
    for(l in 1:L){
      Kilsamp[(NK/Nit)*(it-1)+(1:(NK/Nit)),l]=rpois(NK/Nit,lambdal[l])
    }
  }
  Kilobs=apply(Kilsamp,1,paste0,collapse="-")
  LatTypesChar=sort(unique(Kilobs))
  M=length(LatTypesChar) # number of possible values of the target vector
  LatTypes=matrix(data=NA,ncol=L,nrow=M)
  for(m in 1:M){
    LatTypes[m,]=Kilsamp[which(Kilobs==LatTypesChar[m])[1],]
  }
  # print(LatTypesChar[M])
  return(LatTypes)
}


#' Expected number of sequences per latent type per partition knowing the observation type
#'
#' @param LatTypes (matrix) contains all the possibilities for the number of each sequence type inside a partition.
#' @param Kmaxl (vector of integers) maximal number of sequences of each type in one partition of ddPCR.
#' @param PPKA (numeric matrix) probability of each latent scheme given each observation type.
#'
#' @return a matrix containing expected number of sequences per latent type per partition knowing the observation type.
#' @noRd
EKcompute=function(LatTypes,Kmaxl,PPKA){
  Na=ncol(PPKA)
  L=ncol(LatTypes)
  EK=matrix(data=NA,ncol=L,nrow=Na)
  for(a in 1:Na){
    for(l in 1:L){
      EK[a,l]=0
      for(b in 0:Kmaxl[l]){
        mb=which(LatTypes[,l]==b)
        EK[a,l]=EK[a,l]+b*sum(PPKA[mb,a])
      }
    }
  }
  return(EK)
}
