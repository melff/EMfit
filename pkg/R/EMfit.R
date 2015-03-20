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
                  verbose.inner=FALSE
                  ){
   
  
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
  f.data <- f.phi.data$f.data
  phi.data <- f.phi.data$phi.data
  i <- f.phi.data$i
  
  log_f.ih <- log_f(f.data,psi,...) 
  log_phi.ih <- log_phi(phi.data,psi,...) 
  
  ll.ih <- log_f.ih + log_phi.ih
  
  LL.ih <- exp(ll.ih)
  LL.i <- rowsum(LL.ih,i)
  PPr.ih <- LL.ih/LL.i[i]
  
  Q <- PPr.ih*ll.ih
  Q[PPr.ih==0]<-0
  Q <- sum(Q)
  logLik <- sum(log(LL.i))
  
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
      
      gradient <- crossprod(Jcb.ih,PPr.ih)
      cplInfo <- colSums(PPr.ih*cplInfo.ih)
      
      if(length(constraints)){
        psi <- c(psi + QMat%*%solve(crossprod(QMat,cplInfo)%*%QMat,
                                    crossprod(QMat,gradient)))
      }
      else
        psi <- c(psi + solve(cplInfo,gradient))
      
      log_f.ih <- log_f(f.data,psi,...)
      log_phi.ih <- log_phi(phi.data,psi,...)
      ll.ih <- log_f.ih + log_phi.ih
      
      Q <- PPr.ih*ll.ih
      Q[PPr.ih==0]<-0
      Q <- sum(Q)
      crit.inner <- abs((Q-last.Q)/last.Q)
      
      if(verbose.inner)
      {
        #cat(" - psi:",psi)
        #cat(" - update:",psi.upd)
        #cat(" - Q:",Q)
        cat(" - criterion: ",crit.inner,sep="")
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
    Q <- PPr.ih*ll.ih
    Q[PPr.ih==0]<-0
    Q <- sum(Q)
    logLik <- sum(log(LL.i))
    logLik.trace[iter] <- logLik
    
    crit.outer <- abs((logLik - last.logLik)/last.logLik)
    
    if(verbose){
      
      if(verbose.inner)
        cat("\nIteration ",iter,": ",sep="")
      else
        cat(" -\t")
      cat("Log-likelihood: ",logLik," criterion: ",crit.outer,sep="")
    }
    
    if(crit.outer < eps){
      
      EM.converged <- TRUE
      break
      if(verbose){
        
        cat(" - EM converged")
      }
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
    
    gradient <- crossprod(Jcb.ih,PPr.ih)
    cplInfo <- colSums(PPr.ih*cplInfo.ih)
    
    grad.i <- rowsum(PPr.ih*Jcb.ih,i)
    missInfo <- crossprod(Jcb.ih,PPr.ih*Jcb.ih) - crossprod(grad.i)
    
    obsInfo <- cplInfo - missInfo
    
    if(length(constraints)){
      psi <- c(psi + QMat%*%solve(crossprod(QMat,obsInfo)%*%QMat,
                                  crossprod(QMat,gradient)))
    }
    else
      psi <- c(psi + solve(obsInfo,gradient))
    
    log_f.ih <- log_f(f.data,psi,...)
    log_phi.ih <- log_phi(phi.data,psi,...)
    ll.ih <- log_f.ih + log_phi.ih
    
    LL.ih <- exp(ll.ih)
    LL.i <- rowsum(LL.ih,i)
    PPr.ih <- LL.ih/LL.i[i]
    Q <- PPr.ih*ll.ih
    Q[PPr.ih==0]<-0
    Q <- sum(Q)
    logLik <- sum(log(LL.i))
    
    psi.trace[,iter] <- psi
    crit <- abs((logLik - last.logLik)/last.logLik)
    if(verbose){
      #cat(" - psi:",psi)
      cat(" -\tLog-likelihood: ",logLik," criterion: ",crit,sep="")
    }
    
    if(crit < eps){
      
      NR.converged <- TRUE
      if(verbose){
        
        cat("\n... converged")
      }
      break
    }
    f.phi.data <- expd_data(psi,prev.data=f.phi.data,...)
    f.data <- f.phi.data$f.data
    phi.data <- f.phi.data$phi.data
    i <- f.phi.data$i
  }
  
  list(psi      = psi,
       logLik   = logLik,
       gradient = gradient,
       cplInfo  = cplInfo,
       missInfo  = missInfo,
       obsInfo  = obsInfo,
       converged = if(Information) EM.converged else NR.converged,
       psi.trace=psi.trace
  )
}
  