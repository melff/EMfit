#' A Boilerplate for EM algorithms
#'
#' @description \code{EMfit} provides the boilerplate for EM algorithms and EM-NR accelerated EM algorithms (also known as the Louis-method of acceleration)
#'
#' @param psi.start A numeric vector with starting values
#' @param mk_cpl_data A function that creates the 'complete data' from the 'observed' data. Its return value
#'        can be anything that can be used by the other functions given as arguments.
#'        The return value can also have a component "weights".
#' @param ll_cpl A function that computes the complete-data contributions to the log-likelihood.
#'        It should at least accept the arguments \code{psi}, the parameter vector, and
#'        \code{cpl_data}, the complete-data structure.
#'        It should return the complete-data contributions to the log-likelihood, with an attributed
#'        named "i" that indicates unconditionally independent groups of observations.
#' @param Mstep A function that conducts the M-step. It should accept the parameter vector as
#'        first argument, the complete-data structure as second argument, a vector 
#'        \code{wPPr} of posterior probabilities (weighted if applicable),
#'        and should accept anything that is passed via \dots
#' @param completeInfo A function that computes the 'complete-data' information matrix
#'        It should at least accept the arguments \code{psi}, the parameter vector,
#'        \code{cpl_data}, the complete-data structure, and \code{weights}, which are
#'        weights determined e.g. by posterior probabilities and a-priori weights.
#' @param Jacobian A function that computes the Jacobian of the log-likelihood function.
#'        It should return a matrix with the same number of rows as the length of the result
#'        of \code{ll_cpl} and the same number of colums as elements of \code{psi}
#'        and should have an attribute named "i" that indicates unconditionally independent groups of observations.
#' @param \dots Further arguments passed to \code{ll_cpl}, \code{mk_cpl_data}, etc.
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
  mk_cpl_data,
  ll_cpl,
  completeInfo,
  Jacobian,
  ...,     # further arguments given to the functions
  Mstep=Mstep.default,
  eps=1e-7,
  maxiter=200,
  maxiter.EM=maxiter,
  maxiter.inner=25,
  Information=TRUE,
  verbose=TRUE,
  verbose.inner=FALSE,
  show.psi=FALSE
){

  psi <- psi.start

  cpl_data <- mk_cpl_data(psi,prev.data=NULL,...)
  weights <- cpl_data$weights

  ll.ih <- ll_cpl(psi,cpl_data,...)
  i <- attr(ll.ih,"i")

  LL.ih <- exp(ll.ih)
  LL.i <- rowsum(LL.ih,i)
  PPr.ih <- LL.ih/LL.i[i]

  if(!length(weights))
    weights <- rep(1,length(LL.i))

  logLik <- sum(weights*log(LL.i))
  if(verbose){
    cat("\nInitial log-likelihood:",logLik)
  }
  wPPr <- weights[i]*PPr.ih

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
    if(missing(Mstep))
      psi <- Mstep.default(psi,cpl_data,wPPr,
                           ll_cpl=ll_cpl,
                           Jacobian=Jacobian,
                           completeInfo=completeInfo,...,
                           maxiter=maxiter.innter,eps=eps,
                           verbose=verbose.inner)
    else
      psi <- Mstep(psi,cpl_data,wPPr=wPPr,...,
                   maxiter=maxiter.inner,eps=eps,
                   verbose=verbose.inner)
    converged.inner <- attr(psi,"converged")

    psi.trace[,iter] <- psi

    ### E-step:
    cpl_data <- mk_cpl_data(psi,prev.data=NULL,...)
    weights <- cpl_data$weights

    ll.ih <- ll_cpl(psi,cpl_data,...)
    i <- attr(ll.ih,"i")
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]

    if(!length(weights))
      weights <- rep(1,length(LL.i))

    logLik <- sum(weights*log(LL.i))
    logLik.trace[iter] <- logLik
    wPPr <- weights[i]*PPr.ih

    ### Check convergence:
    crit <- abs((logLik - last.logLik)/last.logLik)
    if(!is.finite(crit)){
      warning("Infinite log-likelihood - stepping back and bailing out ...")
      psi <- last.psi
      logLik <- last.logLik
      break
    }
    if(verbose){

      cat("\tLog-likelihood: ",logLik," criterion: ",crit,sep="")
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

    cplInfo <- completeInfo(psi,cpl_data,weights=wPPr,...)
    Jcb <- Jacobian(psi,cpl_data,...)
    rowsum.Jcb.i <- rowsum(wPPr*Jcb,attr(Jcb,"i"))
    missInfo <- crossprod(Jcb,wPPr*Jcb) - crossprod(rowsum.Jcb.i)
    gradient <- colSums(rowsum.Jcb.i)

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
    weights <- cpl_data$weights

    ll.ih <- ll_cpl(psi,cpl_data,...)
    i <- attr(ll.ih,"i")
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]

    if(!length(weights))
      weights <- rep(1,length(LL.i))

    logLik <- sum(weights*log(LL.i))
    logLik.trace[iter] <- logLik
    wPPr <- weights[i]*PPr.ih

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

        cat("\n... converged\n")
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
#' @param completeInfo A function that computes the 'complete-data' information matrix
#' @param Jacobian A function that computes the Jacobian of the log-likelihood function
#' @param ... Other argumets, passed on to \code{ll_cpl}, etc.
#' @param maxiter Maximal number of iterations
#' @param eps Numeric; a criterion for convergence
#' @param verbose Logical; should an interation trace be displayed?
#'
#' @return A paremeter value that maximizes the Q-function
Mstep.default <- function(psi,cpl_data,
                          wPPr,
                          ll_cpl,
                          completeInfo,
                          Jacobian,
                          ...,
                          maxiter=25,eps=1e-7,
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

    cplInfo <- completeInfo(psi,cpl_data,weights=wPPr,...)
    Jcb <- Jacobian(psi,cpl_data,...)
    gradient <- crossprod(Jcb,wPPr)

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
        cat(" - converged\n")
      converged <- TRUE
      break
    }
  }
  structure(psi,converged=converged)
}
