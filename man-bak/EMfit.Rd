\name{EMfit}
\alias{EMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EMfit(psi.start, Y, X, Z, i, log_f, log_phi, U, ..., eps, maxiter = 200, maxiter.EM = maxiter, maxiter.inner = maxiter, Qmaximizer = NULL, verbose = TRUE, verbose.inner = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi.start}{
%%     ~~Describe \code{psi.start} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{log_f}{
%%     ~~Describe \code{log_f} here~~
}
  \item{log_phi}{
%%     ~~Describe \code{log_phi} here~~
}
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{maxiter.EM}{
%%     ~~Describe \code{maxiter.EM} here~~
}
  \item{maxiter.inner}{
%%     ~~Describe \code{maxiter.inner} here~~
}
  \item{Qmaximizer}{
%%     ~~Describe \code{Qmaximizer} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{verbose.inner}{
%%     ~~Describe \code{verbose.inner} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (psi.start, Y, X, Z, i, log_f, log_phi, U, ..., eps, 
    maxiter = 200, maxiter.EM = maxiter, maxiter.inner = maxiter, 
    Qmaximizer = NULL, verbose = TRUE, verbose.inner = FALSE) 
{
    psi <- psi.start
    log_f.i.u <- log_f(Y, X, i, U, psi, ...)
    log_phi.i.u <- log_phi(U, Z, i, psi, ...)
    ll.i.u <- log_f.i.u + log_phi.u
    log.w <- attr(U, "log-weights")
    LL.i.u <- if (length(log.w)) 
        exp(ll.i.u + log.w)
    else exp(ll.i.u)
    LL.i <- rowSums(LL.i.u)
    PPr.i.u <- LL.i.u/LL.i
    Q <- sum(PPr.i.u * ll.i.u)
    logLik <- sum(log(LL.i))
    psi.trace <- matrix(nrow = length(psi), ncol = maxiter + 
        1)
    psi.trace[, 1] <- psi
    logLik.trace <- numeric(length = maxiter)
    EM.converged <- FALSE
    for (iter in 1:maxiter.EM) {
        if (verbose) {
            cat("\nIteration ", iter, " - sample size = ", nsim, 
                sep = "")
        }
        last.ll.i.u <- ll.i.u
        last.PPr.i.u <- PPr.i.u
        last.logLik <- logLik
        if (is.function(Qmaximizer)) {
            psi <- QMaximizer(Y, X, i, Z, U, psi, log_f, log_phi, 
                verbose = verbose.inner, ...)
            log_f.i.u <- log_f(Y, X, i, U, psi, ...)
            log_phi.i.u <- log_phi(U, Z, i, psi, ...)
            ll.i.u <- log_f.i.u + log_phi.u
            Q <- sum(PPr.i.u * ll.i.u)
        }
        else {
            converged.inner <- FALSE
            for (iiter in 1:maxiter.inner) {
                if (verbose.inner) 
                  cat("\n  Inner iteration ", iter, sep = "")
                last.Q <- Q
                g.i.u <- attr(log_f.i.u, "gradient") + attr(log_phi.i.u, 
                  "gradient")
                g.i.u <- c(PPr.i.u) * g.i.u
                gradient <- colSums(g.i.u, dims = 2)
                negH.i.u <- attr(log_f.i.u, "negHessian")
                +attr(log_phi.i.u, "negHessian")
                negH.i.u <- c(PPr.i.u) * negH.i.u
                negHess <- colSums(negH.i.u, dims = 2)
                psi <- psi + solve(negHess, gradient)
                log_f.i.u <- log_f(Y, X, i, U, psi, ...)
                log_phi.i.u <- log_phi(U, Z, i, psi, ...)
                ll.i.u <- log_f.i.u + log_phi.u
                Q <- sum(PPr.i.u * ll.i.u)
                crit.inner <- abs((Q - last.Q)/Q)
                if (verbose.inner) 
                  cat("criterion: ", crit.inner, sep = "")
                if (crit.inner < eps) {
                  if (verbose.inner) 
                    cat(" - converged")
                  converged.inner <- TRUE
                  break
                }
            }
        }
        log.w <- attr(U, "log-weights")
        LL.i.u <- if (length(log.w)) 
            exp(ll.i.u + log.w)
        else exp(ll.i.u)
        LL.i <- rowSums(LL.i.u)
        PPr.i.u <- LL.i.u/LL.i
        Q <- sum(PPr.i.u * ll.i.u)
        logLik <- sum(log(LL.i))
        logLik.trace[iter] <- logLik
        crit.outer <- abs((logLik - last.logLik)/last.logLik)
        if (verbose) {
            if (verbose.inner) 
                cat("\nIteration ", iter, ": ", sep = "")
            cat("Log-likelihood: ", logLik, " criterion: ", crit, 
                sep = "")
        }
        if (crit.outer < eps) {
            EM.converged <- TRUE
            break
            if (verbose) {
                cat(" - EM converged")
            }
        }
    }
    if (!Information && maxiter.EM == maxiter) 
        return(list(psi = psi, logLik = logLik))
    maxiter.NR <- maxiter - maxiter.EM
    for (iter in 1:maxiter.NR) {
        if (verbose && !(maxiter.EM == maxiter)) {
            cat("\nIteration ", iter, " - sample size = ", nsim, 
                sep = "")
        }
        last.logLik <- logLik
        g.i.u <- attr(log_f.i.u, "gradient") + attr(log_phi.i.u, 
            "gradient")
        P.g.i.u <- c(PPr.i.u) * g.i.u
        gradient <- colSums(P.g.i.u, dims = 2)
        negH.i.u <- attr(log_f.i.u, "negHessian")
        +attr(log_phi.i.u, "negHessian")
        negH.i.u <- c(PPr.i.u) * negH.i.u
        cplInfo <- colSums(negH.i.u, dims = 2)
        G <- tprod(P.g.i.u, 1:2, g.i.u, 1:2)
        g.i <- dimSums(g.i.u, 2)
        misInfo <- G - crossprod(g.i)
        obsInfo <- cplInfo - misInfo
        if (maxiter.NR < 1) 
            break
        psi <- psi + solve(obsInfo, gradient)
        log_f.i.u <- log_f(Y, X, i, U, psi, ...)
        log_phi.i.u <- log_phi(U, Z, i, psi, ...)
        ll.i.u <- log_f.i.u + log_phi.u
        log.w <- attr(U, "log-weights")
        LL.i.u <- if (length(log.w)) 
            exp(ll.i.u + log.w)
        else exp(ll.i.u)
        LL.i <- rowSums(LL.i.u)
        PPr.i.u <- LL.i.u/LL.i
        logLik <- sum(log(LL.i))
        crit <- abs((logLik - last.logLik)/last.logLik)
        if (verbose) {
            if (verbose.inner) 
                cat("\nIteration ", iter, ": ", sep = "")
            cat(" criterion: ", crit, sep = "")
        }
        if (crit < eps) {
            NR.converged <- TRUE
            if (verbose) {
                cat(" - NR converged")
            }
            break
        }
    }
    list(psi = psi, logLik = logLik, gradient = gradient, cplInfo = cplInfo, 
        misInfo = misInfo, obsInfo = obsInfo, converged = if (Information) EM.converged else NR.converged, 
        psi.trace = psi.trace)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
