\name{diagnose}
\alias{diagnose}
\title{Attempt to Diagnose Convergence Problems in a Fitted Mixed Model}
\description{
  \strong{EXPERIMENTAL.} For a given model, this function attempts to
  isolate potential causes of convergence problems. It checks (1) whether
  there are any unusually large coefficients; (2) whether there are any
  unusually scaled predictor variables; (3) if the Hessian (curvature of
  the negative log-likelihood surface at the MLE) is positive definite
  (i.e., whether the MLE really represents an optimum). For each case it
  tries to isolate the particular parameters that are problematic.
}
\usage{
diagnose(fit,
         eval_eps = 1e-5,
         evec_eps = 0.01,
         big_coef = 10,
         big_sd_log10 = 3,
         big_zstat = 5,
         check_coefs = TRUE,
         check_zstats = TRUE,
         check_hessian = TRUE,
         check_scales = TRUE,
         explain = TRUE)
}
\arguments{
  \item{fit}{a fitted \code{\linkS4class{merMod}} object (result of
    \code{\link{lmer}} or \code{\link{glmer}}).}
  \item{eval_eps}{numeric tolerance for detecting 'bad' (near-zero or
    negative) eigenvalues of the Hessian.}
  \item{evec_eps}{numeric tolerance for identifying 'bad' eigenvector
    elements (components with absolute values above this threshold are
    flagged).}
  \item{big_coef}{numeric threshold: coefficients with absolute value
    greater than this are flagged as unusually large.}
  \item{big_sd_log10}{numeric threshold on the log10 scale for flagging
    badly scaled predictors. For the default value of 3, predictors with
    standard deviation less than 1e-3 or greater than 1e3 are flagged.}
  \item{big_zstat}{numeric threshold: Z-statistics (estimate / std. error)
    with absolute value greater than this are flagged.}
  \item{check_coefs}{(logical) identify large-magnitude coefficients?
    For GLMMs using a unitless link function (log, logit, cloglog, or
    probit), both fixed-effect and random-effects (theta) parameters are
    checked; for LMMs (or GLMMs with other links) only random-effects
    (theta) parameters are checked.}
  \item{check_zstats}{(logical) identify parameters with unusually large
    Z-statistics (ratio of estimate to standard error)?  Identifies
    likely failures of Wald confidence intervals/p-values.}
  \item{check_hessian}{(logical) identify non-positive-definite Hessian
    components?  Requires that the model was fitted with
    \code{calc.derivs=TRUE} (the default for smaller models).}
  \item{check_scales}{(logical) identify predictors with unusually small
    or large standard deviations?}
  \item{explain}{(logical) print a detailed explanation of each test?}
}
\value{
  A logical value (returned invisibly) indicating whether the model
  appears to be \code{TRUE} (OK) or \code{FALSE} (problems detected).
}
\details{
  Problems in one category (e.g. complete separation) will generally
  also appear in downstream categories (e.g. non-positive-definite
  Hessians). Therefore, it is generally advisable to deal with problems
  in order: address complete-separation issues first, then re-run the
  diagnostics to see whether Hessian problems persist.

  The Hessian check requires that the model was fitted with
  \code{calc.derivs = TRUE} (which is the default for models with a
  small-to-moderate number of observations and parameters; see
  \code{\link{lmerControl}}).
}
\seealso{
  \code{\link{isSingular}}, \code{\link{allFit}},
  \code{\link{checkConv}}, \link{convergence}
}
\examples{
## well-behaved model - should print "model looks OK!"
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
diagnose(fm1)

## GLMM example
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
diagnose(gm1)
}
