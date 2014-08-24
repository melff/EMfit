\name{MCEMfit}
\alias{MCEMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCEMfit(psi.start, Y, X, Z, i, log_f, log_phi, sim.U, ..., eps, maxiter = 200, maxiter.EM = maxiter, maxiter.inner = maxiter, Information = TRUE, nsim.start, check.sim.alpha = 0.05, check.sim.beta = 0.05, k.update = 3, max.nsim = Inf, verbose = TRUE, verbose.inner = FALSE)
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
  \item{sim.U}{
%%     ~~Describe \code{sim.U} here~~
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
  \item{Information}{
%%     ~~Describe \code{Information} here~~
}
  \item{nsim.start}{
%%     ~~Describe \code{nsim.start} here~~
}
  \item{check.sim.alpha}{
%%     ~~Describe \code{check.sim.alpha} here~~
}
  \item{check.sim.beta}{
%%     ~~Describe \code{check.sim.beta} here~~
}
  \item{k.update}{
%%     ~~Describe \code{k.update} here~~
}
  \item{max.nsim}{
%%     ~~Describe \code{max.nsim} here~~
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
function (psi.start, Y, X, Z, i, log_f, log_phi, sim.U, ..., 
    eps, maxiter = 200, maxiter.EM = maxiter, maxiter.inner = maxiter, 
    Information = TRUE, nsim.start, check.sim.alpha = 0.05, check.sim.beta = 0.05, 
    k.update = 3, max.nsim = Inf, verbose = TRUE, verbose.inner = FALSE) 
{
    psi <- psi.start
    nsim <- nsim.start
    z.alpha <- qnorm(1 - check.sim.alpha)
    z.beta <- qnorm(1 - check.sim.beta)
    U <- sim.U(log_f, log_phi, Y, X, i, Z, psi, NULL, nsim = nsim, 
        ...)
    log_f.i.u <- log_f(Y, X, i, U, psi, ...)
    log_phi.i.u <- log_phi(U, Z, i, psi, ...)
    ll.i.u <- log_f.i.u + log_phi.u
    log.h.u <- attr(U, "log-density")
    log.w <- ll.i.u - log.h.u
    w <- exp(log.w)
    sum.w <- rowSums(w)
    PPr.i.u <- w/sum.w
    psi.trace <- matrix(nrow = length(psi), ncol = maxiter + 
        1)
    psi.trace[, 1] <- psi
    nsim.trace <- numeric(length = maxiter)
    logLik.trace <- numeric(length = maxiter)
    EM.converged <- FALSE
    for (iter in 1:maxiter.EM) {
        if (verbose) {
            cat("\nIteration ", iter, " - sample size = ", nsim, 
                sep = "")
        }
        last.ll.i.u <- ll.i.u
        while (nsim <= max.nsim) {
            if (is.function(Qmaximizer)) {
                psi <- QMaximizer(Y, X, i, Z, U, psi, log_f, 
                  log_phi, verbose = verbose.inner, ...)
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
            delta.ll.i.u <- ll.i.u - last.ll.i.u
            delta.Q.i.u <- delta.ll.i.u * PPr.i.u
            delta.Q.i <- sum.delta.Q.i.u <- rowSums(delta.Q.i.u)
            var.delta.Q.i <- (sum.delta.Q.i.u/sum.w)^2(rowSums(delta.Q.i.u^2)/sum.delta.Q.i.u^2 - 
                2 * rowSums(delta.Q.i.u * w)/(sum.delta.Q.i.u * 
                  sum.w) + rowSums(w^2)/sum.w^2)
            se.delta.Q.i <- sqrt(var.delta.Q.i)
            z.delta.Q.i <- delta.Q.i/se.delta.Q.i
            p.delta.Q.i <- pnorm(p.delta.Q.i)
            nsim.suffic <- (all(1 - p.delta.Q.i < check.sim.alpha))
            if (nsim.suffic) {
                nsim <- nsim * max(1, (z.alpha + z.beta)^2/z.delta.Q.i^2)
            }
            else nsim <- nsim + ceiling(nsim/k.update)
            nsim <- min(nsim, max.nsim)
            U <- sim.U(log_f, log_phi, Y, X, i, Z, psi, NULL, 
                nsim = nsim, ...)
            log_f.i.u <- log_f(Y, X, i, U, psi, ...)
            log_phi.i.u <- log_phi(U, Z, i, psi, ...)
            ll.i.u <- log_f.i.u + log_phi.u
            log.h.u <- attr(U, "log-density")
            log.w <- ll.i.u - log.h.u
            w <- exp(log.w)
            sum.w <- rowSums(w)
            PPr.i.u <- w/sum.w
            if (nsim.suffic) 
                break
            else if (verbose) {
                cat("\nIteration ", iter, ": increased sample size to", 
                  nsim, sep = "")
            }
        }
        delta.Q <- sum(delta.Q.i)
        crit <- abs(delta.Q)/(Q + 0.01)
        ll.i <- log(sum.w)
        logLik <- sum(ll.i)
        nsim.trace[iter] <- nsim
        psi.trace[, iter + 1] <- psi
        logLik.trace[iter] <- logLik
        if (verbose) {
            if (verbose.inner) 
                cat("\nIteration ", iter, ": ", sep = "")
            cat("Log-likelihood: ", logLik, " criterion: ", crit, 
                sep = "")
        }
        if (crit < eps) {
            EM.converged <- TRUE
            break
            if (verbose) {
                cat(" - EM converged")
            }
        }
    }
    if (!Information && maxiter.EM == maxiter) 
        return(list(psi = psi, logLik = logLik, converged = EM.converged, 
            psi.trace = psi.trace))
    NR.converged <- FALSE
    EM.iter <- iter
    for (iter in seq(from = EM.iter + 1, to = maxiter)) {
        if (verbose && !(maxiter.EM == maxiter)) {
            cat("\nIteration ", iter, " - sample size = ", nsim, 
                sep = "")
        }
        last.ll.i.u <- ll.i.u
        last.logLik <- logLik
        while (nsim <= max.nsim) {
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
            delta.ll.i.u <- ll.i.u - last.ll.i.u
            ratio.L.i.u <- exp(delta.ll.i.u)
            w.ratio.L.i.u <- w * ratio.L.i.u
            w.ratio.L.i <- rowSums(w * ratio.L.i.u)
            delta.ll.i <- log(w.ratio.L.i) - log(sum.w)
            var.delta.ll.i <- rowSums(w.ratio.L.i.u^2)/w.ratio.L.i^2 - 
                2 * rowSums(w.ratio.L.i.u * w)/(w.ratio.L.i * 
                  sum.w) + rowSums(w^2)/sum.w^2
            se.delta.ll.i <- sqrt(var.delta.ll.i)
            z.delta.ll.i <- delta.ll.i/se.delta.ll.i
            p.delta.ll.i <- pnorm(p.delta.ll.i)
            nsim.suffic <- (all(1 - p.delta.Q.i < check.sim.alpha))
            if (nsim.suffic) {
                nsim <- nsim * max(1, (z.alpha + z.beta)^2/z.delta.Q.i^2)
            }
            else nsim <- nsim + ceiling(nsim/k.update)
            nsim <- min(nsim, max.nsim)
            U <- sim.U(log_f, log_phi, Y, X, i, Z, psi, NULL, 
                nsim = nsim, ...)
            log_f.i.u <- log_f(Y, X, i, U, psi, ...)
            log_phi.i.u <- log_phi(U, Z, i, psi, ...)
            ll.i.u <- log_f.i.u + log_phi.u
            log.h.u <- attr(U, "log-density")
            log.w <- ll.i.u - log.h.u
            w <- exp(log.w)
            sum.w <- rowSums(w)
            PPr.i.u <- w/sum.w
            if (nsim.suffic) 
                break
            else if (verbose) {
                cat("\nIteration ", iter, ": increased sample size to", 
                  nsim, sep = "")
            }
        }
        nsim.trace[iter] <- nsim
        psi.trace[, iter + 1] <- psi
        logLik.trace[iter] <- logLik
        ll.i <- log(sum.w)
        logLik <- sum(ll.i)
        if (Information) 
            break
        crit <- abs(logLik - last.logLik)/(last.logLik + 0.01)
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
        psi.trace = psi.trace, nsim.trace = nsim.trace)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
