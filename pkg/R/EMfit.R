EMfit <-
function(psi.start,
                  log_f,   
                  log_phi,
                  expd_data,
                  constraints=NULL,
                  enforce.constraints=FALSE,
                  ...,     # further arguments given to the functions
                  eps=1e-7,
                  maxiter=200,
                  maxiter.EM=maxiter,
                  maxiter.inner=maxiter,
                  Information=TRUE,
                  verbose=TRUE,
                  verbose.inner=FALSE,
                  show.psi=FALSE
                  ){
   
  QMat <- NULL
  if(length(constraints)){
     
     C <- constraints$lhs
     if(!length(C)) stop("no left-hand side for constraints given")
     d <- constraints$rhs
     if(!length(d)) d <- numeric(nrow(C))
     
     rstr <- restrictor(C,d)
     constraints.check <- sum(abs(C%*%psi.start-d))
     if(constraints.check>0) {
       if(enforce.constraints){
         
         warning("starting values do not meet constraints -- enforcing them")
         psi.start <- rstr$k + rstr$M%*%psi.start
       }
       else stop("starting values do not meet constraints")
     }
     QMat <- rstr$Q
   }
  psi <- psi.start
  
  f.phi.data <- expd_data(psi,prev.data=NULL,...)
  weights <- f.phi.data$weights # These should be frequency weights or such
  f.data <- f.phi.data$f.data
  phi.data <- f.phi.data$phi.data
  i <- f.phi.data$i
  weights.i <- weights[!duplicated(i)]
  
  log_f.ih <- log_f(f.data,psi,...) 
  log_phi.ih <- log_phi(phi.data,psi,...) 
  
  ll.ih <- log_f.ih + log_phi.ih
  
  LL.ih <- exp(ll.ih)
  LL.i <- rowsum(LL.ih,i)
  PPr.ih <- LL.ih/LL.i[i]
  
  if(!length(weights)) weights <- rep(1,length(PPr.ih))
  Q <- weights*PPr.ih*ll.ih
  Q[PPr.ih==0]<-0
  Q <- sum(Q)
  logLik <- sum(weights.i*log(LL.i))
  if(verbose){
    cat("\nInitial log-likelihood:",logLik)
  }
  
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
    
    last.ll.ih  <- ll.ih
    last.PPr.ih <- PPr.ih
    last.logLik  <- logLik
    last.psi <- psi
    converged.inner <- FALSE
    for(iiter in 1:maxiter.inner){
      
      if(verbose.inner)
        cat("\n  Inner iteration ",iiter,sep="")
      
      last.Q <- Q
      
      Jcb.phi.ih <- attr(log_f.ih,"Jacobian")
      Jcb.f.ih <- attr(log_phi.ih,"Jacobian")
      cplInfo.phi.ih <- attr(log_f.ih,"cplInfo") 
      cplInfo.f.ih <- attr(log_phi.ih,"cplInfo")
      
      Jcb.ih <- Jcb.f.ih + Jcb.phi.ih
      cplInfo.ih <- cplInfo.f.ih + cplInfo.phi.ih
      
      gradient <- crossprod(Jcb.ih,weights*PPr.ih)
      cplInfo <- colSums(weights*PPr.ih*cplInfo.ih)
      
      if(length(constraints)){
        psi <- c(psi + QMat%*%solve(crossprod(QMat,cplInfo)%*%QMat,
                                    crossprod(QMat,gradient)))
      }
      else
        psi <- c(psi + solve(cplInfo,gradient))
      
      log_f.ih <- log_f(f.data,psi,...)
      log_phi.ih <- log_phi(phi.data,psi,...)
      ll.ih <- log_f.ih + log_phi.ih
      
      Q <- weights*PPr.ih*ll.ih
      Q[PPr.ih==0]<-0
      Q <- sum(Q)
      crit.inner <- abs((Q-last.Q)/last.Q)
      
      if(verbose.inner)
      {
        if(show.psi)cat(" - psi:",psi)
        #cat(" - update:",psi.upd)
        #cat(" - Q:",Q)
        cat(" - criterion: ",crit.inner,sep="")
      }
      if(!is.finite(crit.inner)){
        #warning("Non-finite Q-function - stepping back and bailing out ...")
        psi <- last.psi
        logLik <- last.logLik
        break
      }
      if(crit.inner < eps) {
        
        if(verbose.inner)
          cat(" - converged")
        converged.inner <- TRUE
        break
      }
    }
    psi.trace[,iter] <- psi
    
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]
    Q <- weights*PPr.ih*ll.ih
    Q[PPr.ih==0]<-0
    Q <- sum(Q)
    stopifnot(length(weights.i)==length(LL.i))
    logLik <- sum(weights.i*log(LL.i))
    logLik.trace[iter] <- logLik
    
    crit.outer <- abs((logLik - last.logLik)/last.logLik)
    if(!is.finite(crit.outer)){
      warning("Non-finite Q-function - stepping back and bailing out ...")
      psi <- last.psi
      logLik <- last.logLik
      break
    }
    if(verbose){
      
      cat("\n\tLog-likelihood: ",logLik," criterion: ",crit.outer,sep="")
    }
    
    if(crit.outer < eps){
      
      EM.converged <- converged.inner
      
      if(verbose && converged.inner){
        
        cat(" - EM converged")
      }
      break
    }
    f.phi.data <- expd_data(psi,prev.data=f.phi.data,...)
    f.data <- f.phi.data$f.data
    phi.data <- f.phi.data$phi.data
    i <- f.phi.data$i
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
    
    Jcb.ih <- attr(log_f.ih,"Jacobian") + attr(log_phi.ih,"Jacobian")
    cplInfo.ih <- attr(log_f.ih,"cplInfo") + attr(log_phi.ih,"cplInfo")
    
    gradient <- crossprod(Jcb.ih,weights*PPr.ih)
    cplInfo <- colSums(weights*PPr.ih*cplInfo.ih)
    
    grad.i <- rowsum(PPr.ih*Jcb.ih,i)
    missInfo <- crossprod(Jcb.ih,weights*PPr.ih*Jcb.ih) - crossprod(grad.i,weights.i*grad.i)
    
    obsInfo <- cplInfo - missInfo
    
    if(length(constraints)){
      obsInfo.eigen <- eigen(crossprod(QMat,obsInfo)%*%QMat)
    }
    else
      obsInfo.eigen <- eigen(obsInfo)
    last.psi <- psi
    if(any(obsInfo.eigen$values <= 0)) {
      #print(obsInfo.eigen$values)
      warning("observed Information matrix not positive definite -- using complete-data Information")
      if(length(constraints)){
        psi <- c(psi + QMat%*%solve(crossprod(QMat,cplInfo)%*%QMat,
                                    crossprod(QMat,gradient)))
      }
      else
        psi <- c(psi + solve(cplInfo,gradient))
    }
    else 
      {
      if(length(constraints)){
        psi <- c(psi + QMat%*%solve(crossprod(QMat,obsInfo)%*%QMat,
                                    crossprod(QMat,gradient)))
      }
      else
        psi <- c(psi + solve(obsInfo,gradient))
    }
    
    log_f.ih <- log_f(f.data,psi,...)
    log_phi.ih <- log_phi(phi.data,psi,...)
    ll.ih <- log_f.ih + log_phi.ih
    
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]
    Q <- weights*PPr.ih*ll.ih
    Q[PPr.ih==0]<-0
    Q <- sum(Q)
    logLik <- sum(weights.i*log(LL.i))
    logLik.trace[iter] <- logLik
    
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
    f.phi.data <- expd_data(psi,prev.data=f.phi.data,...)
    f.data <- f.phi.data$f.data
    phi.data <- f.phi.data$phi.data
    i <- f.phi.data$i
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
       psi.trace=psi.trace,logLik.trace=logLik.trace,
       QMat = QMat
  )
}
  