#' Target amplification probabilities from ddPCR on full sequences
#'
#' Direct maximization of the Likelihood and of the
#' expected number of sequence per partition
#' with regard to target amplification probabilities
#' in a special case when we know that all sequences cover all the targets
#'
#' @import stats
#'
#' @param NY numerical vector containing the number of each observation type.
#' @param ObsTypes binary matrix containing, for each observation type (in row), if each target is present (1) of absent (0).
#' @param lambda0 (optional) numeric, the average number of sequences in one partition of ddPCR. If not provided, it will be estimated together with the amplification probabilities.
#' @param prec (optional) integer, the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#' @param evallam (optional) bolean, should lambda0 be used as is of should it be estimated along with amplification probabilities? Is automatically switched to TRUE if no lambda0 is provided.
#'
#' @return a vector containing the amplification probability for each target ; also returns the average number of sequences in one partition of ddPCR is evallam=TRUE or if lambda0 is not provided.
#' @export
DirectOptPobs=function(NY,ObsTypes,lambda0=NULL,prec=100,evallam=F){
  if(!is.null(lambda0)){
    if(lambda0>1){warning("The algorithm may give biased results with such high concentration (latent space too big)")}
    if(!evallam){lambda=lambda0}
  }else{
    evallam=T
    lambda0=0.5
    # if(!evallam){stop("if evallam=F (no learning of the expected numebr of sequences per partition),
    #                   a value should be provided for lambda0")}
  }
  if(evallam){
    res=DirectOptPobsL2(NY,ObsTypes,lambda0,prec)
  }else{
    res=DirectOptPobs2(NY,ObsTypes,lambda)
  }
  return(res)
}


#' Direct maximization of the Likelihood
#' with regard to target amplification probabilities
#'
#' @param NY numerical vector containing the number of each observation type.
#' @param ObsTypes binary matrix containing, for each observation type (in row), if each target is present (1) of absent (0).
#' @param lambda numeric, the average number of sequences in one partition of ddPCR.
#'
#' @return a vector containing the amplification probability for each target.
#'
#' @noRd
DirectOptPobs2=function(NY,ObsTypes,lambda){

  J=ncol(ObsTypes)
  Kmax=max(rpois(26000*4,lambda))
  Na=length(NY)

  Llik=function(lpobs){
    pobs=0.5*(1+1/(1+exp(-lpobs)))
    wobs=1-pobs
    lwobs=log(wobs)
    lltot=0
    for(a in 1:Na){
      obsa=ObsTypes[a,]
      m2obsa=2-obsa
      part1=matrix(data=NA,ncol=J,nrow=Kmax+1)
      for(k in 0:Kmax){
        matvalobsk=cbind(log(1-wobs^k),k*lwobs)
        part1[k+1,]=matvalobsk[cbind(1:J,m2obsa)]
      }
      lltot=lltot+NY[a]*logsumexp(apply(part1,1,sum)+dpois(0:Kmax,lambda,log =T))
      # print(lltot)
    }
    return(-lltot)
  }

  Ninit=10
  distvect=numeric(Ninit)*NA
  for(di in 1:Ninit){
    pobs=0.5+runif(J)/2
    wobs=1-pobs
    lpobs=log((pobs-0.5)/wobs)
    resopt=optim(lpobs,Llik)
    distvect[di]=resopt$value
    if(resopt$value<=min(distvect,na.rm=T)){
      lpobsopt=resopt$par
    }
  }

  pobs=0.5*(1+1/(1+exp(-lpobsopt)))
  # pobs=1/(1+exp(-lpobsopt))
  # pobs
  # pobsT

  return(pobs)

}

#' Direct maximization of the Likelihood
#' with regard to target amplification probabilities
#' and expected number of sequences in one partition
#'
#' @param NY numerical vector containing the number of each observation type.
#' @param ObsTypes binary matrix containing, for each observation type (in row), if each target is present (1) of absent (0).
#' @param lambda0 (optional) numeric, the average number of sequences in one partition of ddPCR. If not provided, it will be estimated together with the amplification probabilities.
#' @param prec (optional) integer, the higher the closer the estimated average number of sequences in one partition of ddPCR should be from lambda0 (if provided).
#'
#' @return a vector containing the amplification probability for each target ; also returns the average number of sequences in one partition of ddPCR.
#'
#' @noRd
DirectOptPobsL2=function(NY,ObsTypes,lambda0,prec=150){

  J=ncol(ObsTypes)
  Kmax=max(rpois(26000*4,lambda*1.5))
  Na=length(NY)

  Llik=function(llambdapobs){
    llambda=llambdapobs[1]
    lpobs=llambdapobs[-1]
    # lambda=1/(1+exp(-llambda))
    lambda=exp(llambda)
    pobs=0.5*(1+1/(1+exp(-lpobs)))
    wobs=1-pobs
    lwobs=log(wobs)
    lltot=0
    for(a in 1:Na){
      obsa=ObsTypes[a,]
      m2obsa=2-obsa
      part1=matrix(data=NA,ncol=J,nrow=Kmax+1)
      for(k in 0:Kmax){
        matvalobsk=cbind(log(1-wobs^k),k*lwobs)
        part1[k+1,]=matvalobsk[cbind(1:J,m2obsa)]
      }
      lltot=lltot+NY[a]*logsumexp(apply(part1,1,sum)+dpois(0:Kmax,lambda,log =T))
      # print(lltot)
    }
    return(-lltot)
  }

  Ninit=10
  distvect=numeric(Ninit)*NA
  for(di in 1:Ninit){
    pobs=0.5+runif(J)/2
    wobs=1-pobs
    lpobs=log((pobs-0.5)/wobs)
    lambda=rgamma(1,prec*lambda0,prec)
    # llambda=log(lambda/(1-lambda))
    llambda=log(lambda)
    llambdapobs=c(llambda,lpobs)
    resopt=optim(llambdapobs,Llik)
    distvect[di]=resopt$value
    if(resopt$value<=min(distvect,na.rm=T)){
      llambdapobsopt=resopt$par
    }
  }

  llambdapobs=llambdapobsopt
  llambda=llambdapobs[1]
  lpobs=llambdapobs[-1]
  # lambda=1/(1+exp(-llambda))
  lambda=exp(llambda)
  pobs=0.5*(1+1/(1+exp(-lpobs)))
  # pobs=1/(1+exp(-lpobsopt))
  # pobs
  # pobsT

  return(list(pobs=pobs,lambda=lambda))
}

#' logsumexp
#'
#' @param l vector for which we wish to compute log(sum(exp(l)))
#'
#' @return log(sum(exp(l))) with efficient/reliable computation
#'
#' @noRd
logsumexp=function(l){
  i=which.max(l); res=l[i]+log1p(sum(exp(l[-i]-l[i]))); if (is.nan(res)) res=-Inf; return(res);
}
