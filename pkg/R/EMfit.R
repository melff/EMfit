#' A Boilerplate for EM algorithms
#'
#' @description \code{EMfit} provides the boilerplate for EM algorithms and EM-NR accelerated EM algorithms (also known as the Louis-method of acceleration)
#'
#' @param psi.start A numeric vector with starting values
#' @param ll_cpl A function that computes the complete-data log-likelihood
#' @param mk_cpl_data A function that creates the 'complete data' from the 'observed' data
#' @param \dots Further arguments passed to \code{ll_cpl} and \code{mk_cpl_data}
#' @param eps Numeric; a criterion for convergence
#' @param maxiter Maximal number of iterations
#' @param maxiter.EM Maximal number of iterations after which the algorithm should switch to Newton-Raphson steps
#' @param maxiter.inner Maximal number of iterations for the M-step
#' @param Information Logical; should the observed-data information matrix be returned?
#' @param verbose Logical; should an interation trace be displayed?
#' @param verbose.inner Logical; should the iteration trace of the M-step be displayed (if applicable)?
#' @param show.psi Logical; should the current parameter value be displayed along with the iteration history?
#'
#' @return A list with the following components:
#'  \item{psi}{The MLE of the parameter vector}
#'  \item{logLik}{The maximized observed-data log-likelihood}
#'  \item{gradient}{The gradient of the log-likelihood function}
#'  \item{cplInfo}{The complete-data information matrix}
#'  \item{missInfo}{The missing-data information matrix}
#'  \item{obsInfo}{The observed-data information matrix. Use this to compute standard errors.}
#'  \item{converged}{A logical value, indicating whether the algorithm converged.}
#'  \item{psi.trace}{A matrix which contains the parameter values for each iteration}
#'  \item{logLik.trace}{A vector with the log-likelihood values of each iteration}
#'   
EMfit <- function(
  psi.start,
  ll_cpl,
  mk_cpl_data,
  ...,     # further arguments given to the functions
  eps=1e-7,
  maxiter=200,
  maxiter.EM=maxiter,
  maxiter.inner=25,
  Mstep=Mstep.default,
  Information=TRUE,
  verbose=TRUE,
  verbose.inner=FALSE,
  show.psi=FALSE
){

  psi <- psi.start

  cpl_data <- mk_cpl_data(psi,prev.data=NULL,...)
  weights.i <- cpl_data$weights.i # These should be frequency weights or such for each set of indiv obs
  i <- cpl_data$i  # indicates independent observations

  ll.ih <- ll_cpl(psi,cpl_data,...)

  LL.ih <- exp(ll.ih)
  LL.i <- rowsum(LL.ih,i)
  PPr.ih <- LL.ih/LL.i[i]
  logLik <- sum(weights.i*log(LL.i))
  if(verbose){
    cat("\nInitial log-likelihood:",logLik)
  }
  wPPr <- weights.i[i]*PPr.ih

  psi.trace <- matrix(nrow=length(psi),ncol=maxiter+1)
  psi.trace[,1] <- psi
  logLik.trace <- numeric(length=maxiter)
  #    if(verbose){
  #      cat("\nInitial psi:",psi)
  #    }
  EM.converged <- FALSE
  for(iter in 1:maxiter.EM){

    if(verbose){
      cat("\nEM Iteration ",iter,sep="")
    }

    last.logLik  <- logLik
    last.psi <- psi
    converged.inner <- FALSE

    ### M-step:
    psi <- Mstep(psi,cpl_data,wPPr,ll_cpl=ll_cpl,...,
                           maxiter=maxiter.innter,eps=eps,
                           verbose=verbose.inner)
    converged.inner <- attr(psi,"converged")

    psi.trace[,iter] <- psi

    ### E-step:
    cpl_data <- mk_cpl_data(psi,prev.data=NULL,...)
    weights.i <- cpl_data$weights.i # These should be frequency weights or such for each set of indiv obs
    i <- cpl_data$i  # indicates independent observations

    ll.ih <- ll_cpl(psi,cpl_data,...)
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]
    stopifnot(length(weights.i)==length(LL.i))
    logLik <- sum(weights.i*log(LL.i))
    logLik.trace[iter] <- logLik
    wPPr <- weights.i[i]*PPr.ih

    ### Check convergence:
    crit <- abs((logLik - last.logLik)/last.logLik)
    if(!is.finite(crit)){
      warning("Infinite log-likelihood - stepping back and bailing out ...")
      psi <- last.psi
      logLik <- last.logLik
      break
    }
    if(verbose){

      cat("\n\tLog-likelihood: ",logLik," criterion: ",crit,sep="")
    }

    if(crit < eps){

      EM.converged <- converged.inner

      if(verbose && converged.inner){

        cat(" - EM converged")
      }
      break
    }
  }

  if(!Information && maxiter.EM==maxiter)
    return(
      list(psi      = psi,
           logLik   = logLik
      ))

  last.iter.EM <- iter

  for(iter in (last.iter.EM+1):maxiter){ # NR acceleration and obs. info

    if(verbose &&!(maxiter.EM==maxiter)){
      cat("\nNR Iteration ",iter,sep="")
    }

    last.logLik <- logLik

    gradient <- attr(ll.ih,"gradient")
    if(!length(gradient)){
      Jcb.ih <- attr(ll.ih,"Jacobian")
      gradient <- colSums(Jcb.ih)
    }

    cplInfo <- attr(ll.ih,"cplInfo")
    missInfo <- attr(ll.ih,"missInfo")
    if(!length(missInfo)){
      Jcb.ih <- attr(ll.ih,"Jacobian")
      if(!length(Jcb.ih)) stop("need either missing information or Jacobian")
      grad.i <- rowsum(PPr.ih,i)
      missInfo <- crossprod(Jcb.ih/wPPr,Jcb.ih) - crossprod(grad.i/weights.i,grad.i)
      # Gradients and Jacobians are already weighted, we
      # need to compensate
    }

    obsInfo <- cplInfo - missInfo

    obsInfo.eigen <- eigen(obsInfo)
    last.psi <- psi
    if(any(obsInfo.eigen$values <= 0)) {
      #print(obsInfo.eigen$values)
      warning("observed Information matrix not positive definite -- using complete-data Information")
      psi <- c(psi + solve(cplInfo,gradient))
    }
    else {
      psi <- c(psi + solve(obsInfo,gradient))
    }

    cpl_data <- mk_cpl_data(psi,prev.data=NULL,...)
    weights.i <- cpl_data$weights.i # These should be frequency weights or such for each set of indiv obs
    i <- cpl_data$i  # indicates independent observations

    ll.ih <- ll_cpl(psi,cpl_data,...)
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]
    logLik <- sum(weights.i*log(LL.i))
    logLik.trace[iter] <- logLik
    wPPr <- weights.i[i]*PPr.ih

    psi.trace[,iter] <- psi
    crit <- (logLik - last.logLik)/abs(last.logLik)
    if(verbose){
      if(show.psi)cat(" - psi:",psi)
      cat(" -\tLog-likelihood: ",logLik," criterion: ",crit,sep="")
    }
    if(!is.finite(crit)){
      warning("Non-finite log-likelihood - stepping back and bailing out ...")
      psi <- last.psi
      logLik <- last.logLik
      break
    }
    if(abs(crit) < eps){

      NR.converged <- TRUE
      if(verbose){

        cat("\n... converged")
      }
      break
    }
    else if(crit < 0){
      warning("Cannot increase likelihood - stepping back and bailing out ...")
      psi <- last.psi
      logLik <- last.logLik
      break
    }
  }
  psi.trace <- psi.trace[,1:iter,drop=FALSE]
  logLik.trace <- logLik.trace[1:iter]

  list(psi      = psi,
       logLik   = logLik,
       gradient = gradient,
       cplInfo  = cplInfo,
       missInfo  = missInfo,
       obsInfo  = obsInfo,
       converged = if(Information) EM.converged else NR.converged,
       psi.trace=psi.trace,
       logLik.trace=logLik.trace
  )
}


#' Default Maximizer for the M-Step 
#' 
#' This function provides a standard Newton-Raphson maximizer of the complete data log-likelihood.
#'
#' @param psi A numeric vector with starting values
#' @param cpl_data A data structure with complete data
#' @param wPPr A vector with (weighted) posterior probabilities
#' @param ll_cpl  A function that computes the complete-data log-likelihood
#' @param ... Other argumets, passed on to \code{ll_cpl}
#' @param maxiter Maximal number of iterations
#' @param eps Numeric; a criterion for convergence
#' @param verbose Logical; should an interation trace be displayed?
#'
#' @return A paremeter value that maximizes the Q-function
Mstep.default <- function(psi,cpl_data,wPPr,ll_cpl,...,maxiter=25,eps=1e-7,
                               verbose=FALSE){


  last.psi <- psi
  converged <- FALSE
  ll.ih <- ll_cpl(psi,cpl_data,...)
  Q <- wPPr*ll.ih
  Q[wPPr==0]<-0
  Q <- sum(Q)

  for(iter in 1:maxiter){

    if(verbose)
      cat("\n  Inner iteration ",iter,sep="")

    last.Q <- Q

    gradient <- attr(ll.ih,"gradient")
    if(!length(gradient)){
      Jcb.ih <- attr(ll.ih,"Jacobian")
      gradient <- colSums(Jcb.ih)
    }

    cplInfo <- attr(ll.ih,"cplInfo")

    psi <- c(psi + solve(cplInfo,gradient))

    ll.ih <- ll_cpl(psi,cpl_data,...)

    Q <- wPPr*ll.ih
    Q[PPr.ih==0]<-0
    Q <- sum(Q)
    crit <- abs((Q-last.Q)/last.Q)

    if(verbose)
    {
      #cat(" - update:",psi.upd)
      #cat(" - Q:",Q)
      cat(" - criterion: ",crit,sep="")
    }
    if(!is.finite(crit)){
      warning("Non-finite Q-function - stepping back and bailing out ...")
      psi <- last.psi
      break
    }
    if(crit < eps) {

      if(verbose)
        cat(" - converged")
      converged <- TRUE
      break
    }
  }
  structure(psi,converged=converged)
}
