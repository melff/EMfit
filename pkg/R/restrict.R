restrictor <- function(C,d=numeric(m),sign=7){

  ## Create a matrix that maps a short vector of
  ## linearly unrestricted parameters
  ## to a onger vector of linearly restricted
  ## parameters

  if(!is.matrix(C)) C <- t(as.vector(C))

  m <- nrow(C)
  n <- ncol(C)
  ginv.C <- t(solve(tcrossprod(C),C))

  M <- diag(n) - ginv.C %*% C
  QRM <- QR(M)
  list(
    Q=round(QRM$Q,sign),
    M=M,
    k=ginv.C%*%d,
    C=C,
    ginv.C=ginv.C,
    d=d
    )
}

QR <- function(M){
  ## standardised QR decomposition:
  ## diagonal elements of R are always
  ## positive, only "significant"
  ## columns in Q are returned
  
  qrM <- qr(M)
  rnk <- qrM$rank
  Q <- qr.Q(qrM)[,1:rnk,drop=FALSE]
  R <- qr.R(qrM)[1:rnk,,drop=FALSE]
  
  sgndR <- sign(diag(R))
  sgndR <- diag(x=sgndR,nrow=length(sgndR))
  Q <- Q%*%sgndR
  R <- sgndR%*%R
  
  list(Q=Q, R=R)
}

