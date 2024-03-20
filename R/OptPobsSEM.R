#' SEM algorithm for amplification probability learning
#'
#' Learning of target-specific amplification probability when sequence types proportions are known
#'
#' @param NY (numerical vector) contains the number of each observation type.
#' @param ObsTypes (binary matrix) contains, for each observation type (in row), if each target is present (1) of absent (0).
#' @param propl (vector) sequence types frequencies.
#' @param coverseq (numeric matrix) targets covered by each of those sequence types.
#' @param LatTypes (optional) matrix, contains all the possibilities for the number of each sequence type inside a partition.
#' @param lambda0 (optional) numeric, the average number of sequences in one partition of ddPCR. If not provided, it will be estimated together with the amplification probabilities.
#' @param prec (optional) integer, the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#' @param evallam (optional) bolean, should lambda0 be used as is of should it be estimated along with amplification probabilities? Is automatically switched to TRUE if no lambda0 is provided.
#' @param Nit (optional) integer, maximal number of iteration in the SEM algorithm
#' @param plotopt (optional) bolean, should a parameter convergence plot be drawn?
#' @param pobsinit (optional, numerical vector) initial values for target-specific amplification probabilities
#'
#' @return target-specific amplification probabilities ($pobs),
#' their values all along the SEM ($POBS),
#' if lambda is not provided, its values all along the SEM ($LAMBDA),
#' the complete log-likelohood all along the SEM ($LLC)
#' @export
OptPobsSEM=function(NY,ObsTypes,propl,coverseq,
                    LatTypes=NULL,lambda0=NULL,prec=100,evallam=F,
                    Nit=30,plotopt=T,pobsinit=NULL){
  #### INITIALIZATION ####
  J=ncol(ObsTypes)
  L=nrow(coverseq)
  if(is.null(lambda0)){
    lambda=rgamma(1,prec*lambda0,prec)
    evallam=T
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
      LatTypes=LatLibraryPobsLam(Nit=floor(N/2)*2,NK=floor(N/2)*6,lambda0=lambda0*2,prec=prec,propl=propl)
    }else{
      LatTypes=LatLibraryPobs(Nit=floor(N/2)*2,NK=floor(N/2)*6,lambda=lambda0*2,propl=propl)
    }
  }
  if(evallam){
    res=OptLamPobs(NY,LatTypes=LatTypes,coverseq=coverseq,ObsTypes=ObsTypes,propl=propl,Nit=Nit,
                            plotopt=plotopt,pobsinit=pobsinit,lambda0=lambda,prec=prec)
  }else{
    res=OptPobs(NY,LatTypes=LatTypes,coverseq=coverseq,ObsTypes=ObsTypes,propl=propl,Nit=Nit,
                plotopt=plotopt,pobsinit=pobsinit,lambda=lambda)
  }
  return(res)
}


#' SEM algorithm for amplification probability (only) learning
#'
#' @param NY (numerical vector) contains the number of each observation type.
#' @param LatTypes (optional) matrix, contains all the possibilities for the number of each sequence type inside a partition.
#' @param coverseq (numeric matrix) targets covered by each of those sequence types.
#' @param ObsTypes (binary matrix) contains, for each observation type (in row), if each target is present (1) of absent (0).
#' @param propl (vector) sequence types frequencies.
#' @param lambda numeric, the average number of sequences in one partition of ddPCR.
#' @param Nit (optional) integer, maximal number of iteration in the SEM algorithm
#' @param plotopt (optional) bolean, should a parameter convergence plot be drawn?
#' @param pobsinit (optional, numerical vector) initial values for target-specific amplification probabilities
#'
#' @return target-specific amplification probabilities ($pobs),
#' their values all along the SEM ($POBS),
#' the complete log-likelohood all along the SEM ($LLC)
#' @noRd
OptPobs=function(NY,LatTypes,coverseq,ObsTypes,propl,lambda,Nit=30,plotopt=T,pobsinit=NULL){

  LatTypesj=LatTypes%*%coverseq
  L=ncol(LatTypes)
  J=ncol(ObsTypes)

  lambdal=lambda*propl
  PK=PKcompute(LatTypes,lambdal)

  #### INITIALIZATION ####
  if(is.null(pobsinit)){
    pobs=0.5+runif(J)/2
  }else{
    pobs=pobsinit
  }
  wobs=1-pobs

  #### SEM Algorithm ####
  POBS=matrix(data=NA,ncol=J,nrow=Nit)
  LLC=numeric(Nit)*NA
  for(it in 1:Nit){

    POBS[it,] = pobs


    #### SAMPLING ####

    # Probability of each observation
    # given each possible latent scheme:
    PKA=PKAcompute(LatTypesj,ObsTypes,wobs)

    # Probability of each latent scheme
    # given each observation type and a priori probabilities
    PPKA=PK*PKA
    # !! up to a constant !! => renormalization
    sumPPKA=apply(PPKA,2,sum)
    PPKA=t(t(PPKA)/sumPPKA)

    # Drawing a Kil sample from PPKA:
    # only keeping the coordinates of
    # the sample-observation combination:
    iKil=sampleiKil(PPKA,NY)

    lpobs=log((pobs-0.5)/wobs)


    #### MAXIMIZATION ####

    maxemlobs=function(lpobs){
      pobs=0.5*(1+1/(1+exp(-lpobs)))
      wobs=1-pobs
      # Probability of each observation
      # given each possible latent scheme:
      PKA=PKAcompute(LatTypesj,ObsTypes,wobs)
      LLc=LLcomputeNew2(PK,PKA,iKil)
      return(-LLc)
    }
    resopt=optim(par=lpobs,fn=maxemlobs)
    lpobsopt=resopt$par
    pobs=0.5*(1+1/(1+exp(-lpobsopt)))
    pobs[pobs==1] <- 0.999
    pobs[pobs==0.5] <- 0.501
    wobs=1-pobs
    if(it%%20==0){
      print(it)
      print(pobs)
    }
    LLc=-resopt$value
    LLC[it]=LLc

  }

  if(plotopt){
    plot(LLC)
    matplot(POBS,type="l")#ylim=range(rbind(POBS,pobsT)))
    # abline(h=pobsT,lty=1:(J),col=alpha(1:(J),0.3),lwd=2)
  }

  return(list(pobs=pobs,POBS=POBS,LLC=LLC))

}

#' SEM algorithm for amplification probability and lambda learning
#'
#' @param NY (numerical vector) contains the number of each observation type.
#' @param LatTypes (optional) matrix, contains all the possibilities for the number of each sequence type inside a partition.
#' @param coverseq (numeric matrix) targets covered by each of those sequence types.
#' @param ObsTypes (binary matrix) contains, for each observation type (in row), if each target is present (1) of absent (0).
#' @param propl (vector) sequence types frequencies.
#' @param Nit (optional) integer, maximal number of iteration in the SEM algorithm
#' @param plotopt (optional) bolean, should a parameter convergence plot be drawn?
#' @param pobsinit (optional, numerical vector) initial values for target-specific amplification probabilities
#' @param lambda0 numeric, the average number of sequences in one partition of ddPCR. If not provided, it will be estimated together with the amplification probabilities.
#' @param prec integer, the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#'
#' @return target-specific amplification probabilities ($pobs),
#' their values all along the SEM ($POBS),
#' lambda values all along the SEM ($LAMBDA),
#' the complete log-likelohood all along the SEM ($LLC)
#' @noRd
OptLamPobs=function(NY,LatTypes,coverseq,ObsTypes,propl,Nit=30,
                    plotopt=T,pobsinit=NULL,lambda0=0.45,prec=100){

  LatTypesj=LatTypes%*%coverseq
  L=ncol(LatTypes)
  J=ncol(ObsTypes)

  #### INITIALIZATION ####
  lambda=rgamma(1,prec*lambda0,prec)
  if(is.null(pobsinit)){
    pobs=0.5+runif(J)/2
  }else{
    pobs=pobsinit
  }
  wobs=1-pobs

  #### SEM Algorithm ####
  POBS=matrix(data=NA,ncol=J,nrow=Nit)
  LAMBDA=numeric(Nit)*NA
  LLC=numeric(Nit)*NA
  for(it in 1:Nit){

    POBS[it,] = pobs
    LAMBDA[it] = lambda

    lambdal=lambda*propl

    #### SAMPLING ####

    # # A priori probability of each possible
    # # latent scheme (within LatTypesChar)
    # # (given the Poisson parameter and the sequences proportions):
    PK=PKcompute(LatTypes,lambdal)

    # Probability of each observation
    # given each possible latent scheme:
    PKA=PKAcompute(LatTypesj,ObsTypes,wobs)

    # Probability of each latent scheme
    # given each observation type
    PPKA=PK*PKA
    # !! up to a constant !! => renormalization
    sumPPKA=apply(PPKA,2,sum)
    PPKA=t(t(PPKA)/sumPPKA)

    # Drawing a Kil sample from PPKA:
    # only keeping the coordinates of
    # the sample-observation combination:
    iKil=sampleiKil(PPKA,NY)

    #### MAXIMIZATION ####
    maxemlobs=function(llambdalpobs){
      llambda=llambdalpobs[1]
      lpobs=llambdalpobs[-1]
      lambda=1/(1+exp(-llambda))
      pobs=0.5*(1+1/(1+exp(-lpobs)))
      wobs=1-pobs
      lambdal=lambda*propl
      # Probability of each observation
      # given each possible latent scheme:
      PKA=PKAcompute(LatTypesj,ObsTypes,wobs)
      # A priori probability of each possible
      # latent scheme (within LatTypesChar)
      # (given the Poisson parameter and the sequences proportions):
      PK=PKcompute(LatTypes,lambdal)
      LLc=LLcomputeNew2(PK,PKA,iKil)
      return(-LLc)
    }
    lpobs=log((pobs-0.5)/wobs)
    llambda=log(lambda/(1-lambda))
    llambdapobs=c(llambda,lpobs)
    resopt=optim(par=llambdapobs,fn=maxemlobs)
    llambdapobsopt=resopt$par
    llambda=llambdapobsopt[1]
    lpobsopt=llambdapobsopt[-1]
    pobs=0.5*(1+1/(1+exp(-lpobsopt)))
    lambda=1/(1+exp(-llambda))
    pobs[pobs==1] <- 0.999
    pobs[pobs==0.5] <- 0.501
    wobs=1-pobs
    LLc=-resopt$value
    if(it%%25==0){
      print(it)
      print(pobs)
      print(lambda)
      print(LLc)
    }

    LLC[it]=LLc

  }

  if(plotopt){
    par(mfcol=c(1,3))
    # matplot(cbind(POBS,LAMBDA),type="l",ylim=range(cbind(POBS,LAMBDA,pobsT,lambdaT)))
    # abline(h=c(pobsT,lambdaT),lty=1:(J+1),col=alpha(1:(J+1),0.3),lwd=2)
    plot(LLC)
    plot(LAMBDA,type="l")
    matplot(POBS,type="l")#ylim=range(rbind(POBS,pobsT)))
    # abline(h=pobsT,lty=1:(J),col=alpha(1:(J),0.3),lwd=2)
  }

  return(list(pobs=pobs,POBS=POBS,LAMBDA=LAMBDA,LLC=LLC))

}


#' Drawing a Kil sample from PPKA
#'
#' only keeping the coordinates of the sample-observation combination
#'
#' @param PPKA (numeric matrix) probability of each latent scheme given each observation type.
#' @param NY (numerical vector) contains the number of each observation type.
#'
#' @return Kil sample from PPKA
#' @noRd
sampleiKil=function(PPKA,NY){
  A=ncol(PPKA)
  M=nrow(PPKA)
  iKil=matrix(data=NA,ncol=A,nrow=M)
  # a0=0
  for(a in 1:A){
    #   indexa=a0+(1:NY[a])
    #   iKil[indexa,1]=a
    #   iKil[indexa,2]=rmultinom(1,NY[a],PPKA[,a])
    #   a0=a0+NY[a]
    iKil[,a]=rmultinom(1,NY[a],PPKA[,a])
  }
  return(iKil)
}


#' Complete Log likelihood computation
#'
#' @param PK (numeric matrix) containing the a priori probability of each possible latent scheme.
#' @param PKA (numeric matrix) containing the probability of each observation given each possible latent scheme.
#' @param iKil coordinates of samples drawn at the S step
#'
#' @return complete-data log-likelihood
#' @noRd
LLcomputeNew2=function(PK,PKA,iKil){
  # # Probability of each latent scheme
  # # given each observation type
  # PPKA=PK*PKA
  # # !! up to a constant !! => renormalization
  # #sumPPKA=apply(PPKA,2,sum)
  # # PPKA[,sumPPKA==0]=1
  # PPKA=t(t(PPKA)/apply(PPKA,2,sum))
  # # PPKA=t(t(PPKA)/apply(PPKA,2,sum))
  #
  # PPKA0=PPKA
  # PPKA0[PPKA==0] <- 1
  # PKA0=PKA
  # PKA0[PKA==0] <- 1
  # lPKA=log(PKA0)
  lPKA=log(PKA^iKil)
  lPKA[is.infinite(lPKA)] = min(lPKA[!is.infinite(lPKA)])
  LLc1=sum(lPKA)

  # lPPKAIK=log(PPKA^iKil)
  # lPPKAIK[is.infinite(lPPKAIK)] = 0#(-100)
  # LLc1=sum(lPPKAIK)


  NK=apply(iKil,1,sum)
  no0=which(NK!=0)
  # names(NK)=LatTypesChar
  LLc2=sum(NK[no0]*log(PK[no0]))

  LLc=LLc1+LLc2

  return(LLc)
}

#' Latente types library building
#'
#' When proportions of sequence types are unknown
#'
#' @param Nit (integer), number of sequence type proportion generation.
#' @param NK (integer), number of latent vector generation per generated proportion. Nk has to be a multiple of Nit.
#' @param lambda (numeric) average number of sequences in one partition of ddPCR.
#' @param propl (vector) sequence types frequencies.
#'
#' @return a matrix containing all the possibilities for the number of each sequence type inside a partition.
#' @export
LatLibraryPobs=function(Nit=1000,NK=100000,lambda,propl){
  lambdal=lambda*propl
  L=length(propl)
  Kilsamp=matrix(data=NA,ncol=L,nrow=NK)
  if(NK%%Nit!=0){stop("NK HAS to be a multiple of Nit")}
  for(it in 1:Nit){
    #lpropl=rnorm(L-1,0,3)
    #propl=c(1/(1+sum(exp(lpropl))),exp(lpropl)/(1+sum(exp(lpropl))))
    # print(round(propl,1))
    #lambdal=lambda*propl
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
#' When proportions of sequence types and average number of sequences per partition are unknown
#'
#' @param Nit (integer), number of sequence type proportion generation.
#' @param NK (integer), number of latent vector generation per generated proportion. Nk has to be a multiple of Nit.
#' @param lambda0 (numeric), a starting point for average number of sequences in one partition of ddPCR.
#' @param prec (integer), the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#' @param propl (vector) sequence types frequencies.
#'
#' @return a matrix containing all the possibilities for the number of each sequence type inside a partition.
#' @export
LatLibraryPobsLam=function(Nit=1000,NK=100000,lambda0=0.45,prec=100,propl){
  L=length(propl)
  lambdasamp=rgamma(Nit,prec*lambda0,prec)
  Kilsamp=matrix(data=NA,ncol=L,nrow=NK)
  if(NK%%Nit!=0){stop("NK HAS to be a multiple of Nit")}
  for(it in 1:Nit){
    lambdal=lambdasamp[it]*propl
    #lpropl=rnorm(L-1,0,3)
    #propl=c(1/(1+sum(exp(lpropl))),exp(lpropl)/(1+sum(exp(lpropl))))
    # print(round(propl,1))
    #lambdal=lambda*propl
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
