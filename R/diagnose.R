##' Attempt to diagnose convergence problems in a fitted mixed model
##'
##' \strong{EXPERIMENTAL}. For a given model, this function attempts to isolate
##' potential causes of convergence problems. It checks (1) whether there are
##' any unusually large coefficients; (2) whether there are any unusually
##' scaled predictor variables; (3) if the Hessian (curvature of the
##' negative log-likelihood surface at the MLE) is positive definite
##' (i.e., whether the MLE really represents an optimum). For each
##' case it tries to isolate the particular parameters that are problematic.
##'
##' Problems in one category (e.g. complete separation) will generally
##' also appear in "downstream" categories (e.g. non-positive-definite
##' Hessians). Therefore, it is generally advisable to try to deal with
##' problems in order, e.g. address problems with complete separation first,
##' then re-run the diagnostics to see whether Hessian problems persist.
##'
##' @param fit a fitted \code{merMod} object (from \code{lmer}/\code{glmer})
##' @param eval_eps numeric tolerance for 'bad' eigenvalues
##' @param evec_eps numeric tolerance for 'bad' eigenvector elements
##' @param big_coef numeric tolerance for large coefficients
##' @param big_sd_log10 numeric tolerance for badly scaled parameters
##'   (log10 scale): for the default value of 3, predictor variables with
##'   standard deviation less than 1e-3 or greater than 1e3 will be flagged
##' @param big_zstat numeric tolerance for Z-statistic
##' @param check_coefs (logical) identify large-magnitude coefficients?
##'   For GLMMs using a unitless link function (log, logit, cloglog, or
##'   probit), both fixed-effect and random-effects (theta) parameters are
##'   checked, since the coefficients are on a unitless scale and values
##'   larger than \code{big_coef} are inherently suspicious.  For LMMs
##'   or GLMMs with other links, only random-effects (theta) parameters are
##'   checked, since fixed-effect magnitudes depend on the scale of the
##'   response. May produce false positives if predictor variables have
##'   extremely large scales.
##' @param check_zstats (logical) identify parameters with unusually large
##'   Z-statistics (ratio of estimate to standard error)?  Identifies
##'   likely failures of Wald confidence intervals/p-values.
##' @param check_hessian (logical) identify non-positive-definite Hessian
##'   components?  Requires that the model was fitted with
##'   \code{calc.derivs=TRUE} (the default for smaller models).
##' @param check_scales (logical) identify predictors with unusually small
##'   or large standard deviations?
##' @param explain (logical) provide a detailed explanation of each test?
##' @return a logical value (invisibly) indicating whether anything
##'   questionable was found (\code{FALSE} means problems were detected)
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
##' diagnose(fm1)
##' @export
diagnose <- function(fit,
                     eval_eps = 1e-5,
                     evec_eps = 0.01,
                     big_coef = 10,
                     big_sd_log10 = 3,
                     big_zstat = 5,
                     check_coefs = TRUE,
                     check_zstats = TRUE,
                     check_hessian = TRUE,
                     check_scales = TRUE,
                     explain = TRUE) {

    if (!is(fit, "merMod")) {
        stop("'fit' must be a merMod object (from lmer or glmer)")
    }

    prt_explain <- function(...) {
        if (explain) {
            s <- paste0(...)
            cat(strwrap(s), "\n", sep = "\n")
        }
        invisible(NULL)
    }

    ## link functions for which fixed-effect coefficients are unitless
    ## (and therefore comparable to big_coef without adjustment)
    unitless_links <- c("log", "cloglog", "logit", "probit")

    model_OK <- TRUE

    ## extract fixed effects
    ff <- fixef(fit)
    nn <- names(ff)

    ## extract random-effects parameters (theta) with names
    theta <- getME(fit, "theta")
    n_theta <- length(theta)

    ## --- Check large coefficients ---
    if (check_coefs) {
        ## For GLMMs with unitless link functions (log, logit, etc.),
        ## large fixed-effect coefficients may indicate complete separation.
        ## For LMMs (or GLMMs with identity link), fixed-effect magnitudes
        ## are on the scale of the response so we do not flag them.
        check_beta <- FALSE
        if (is(fit, "glmerMod")) {
            lnk <- family(fit)$link
            check_beta <- lnk %in% unitless_links
        }

        bigcoef <- if (check_beta) ff[abs(ff) > big_coef] else numeric(0)
        bigtheta <- theta[abs(theta) > big_coef]
        all_bigcoef <- c(bigcoef, bigtheta)

        if (length(all_bigcoef) > 0) {
            model_OK <- FALSE
            cat(sprintf("Unusually large coefficients (|x|>%g):\n\n", big_coef))
            print(all_bigcoef)
            cat("\n")
            prt_explain(
                "Large (in magnitude) random-effects parameters (Cholesky factor ",
                "elements) suggest the RE variance is very large or very small ",
                "relative to the residual variance; in the former case the model may ",
                "be overfitted, in the latter case unnecessary components may be ",
                "converging to zero on the constrained scale. Large fixed-effect ",
                "coefficients in binomial or Poisson models suggest ",
                "(quasi-)complete separation."
            )
        }
    }

    ## --- Check predictor scales ---
    if (check_scales) {
        ## use model.matrix() rather than getME("X") to get unscaled predictors
        X <- model.matrix(fit)
        sdvec <- apply(X, 2, sd)
        ## exclude zero-variance columns (typically intercepts)
        sdvec <- sdvec[sdvec > 0 & abs(log10(sdvec)) > big_sd_log10]
        if (length(sdvec) > 0) {
            model_OK <- FALSE
            cat(sprintf(
                "\nPredictors with unusually large or small standard deviations (|log10(sd)|>%g):\n\n",
                big_sd_log10
            ))
            print(sdvec)
            cat("\n")
            prt_explain(
                "Predictor variables with very narrow or wide ranges generally give ",
                "rise to parameters with very large or small magnitudes, which can ",
                "exacerbate numerical instability and may incorrectly suggest a ",
                "non-positive-definite Hessian."
            )
        }
    }

    ## --- Check Z-statistics ---
    if (check_zstats) {
        vv <- tryCatch(as.matrix(vcov(fit)), error = function(e) NULL)
        if (!is.null(vv) && nrow(vv) > 0) {
            se <- sqrt(diag(vv))
            z <- ff / se
            bigz <- z[!is.na(z) & abs(z) > big_zstat]
            if (length(bigz) > 0) {
                model_OK <- FALSE
                cat(sprintf("Unusually large Z-statistics (|x|>%g):\n\n", big_zstat))
                print(bigz)
                cat("\n")
                prt_explain(
                    "Large Z-statistics (estimate/std. error) suggest a possible ",
                    "failure of the Wald approximation - often associated with ",
                    "parameters near the boundary of their feasible range ",
                    "(e.g. random-effects standard deviations approaching 0). ",
                    "While Wald p-values and standard errors may be unreliable, ",
                    "likelihood ratio tests (e.g. from drop1()) are probably still OK."
                )
            }
        }
    }

    ## --- Check Hessian ---
    if (check_hessian) {
        H <- fit@optinfo$derivs$Hessian
        if (is.null(H)) {
            cat(paste0(
                "Hessian not available; skipping Hessian check.\n",
                "(Refit with calc.derivs=TRUE in lmerControl()/glmerControl()",
                " to enable this check.)\n\n"
            ))
        } else {
            eigs <- tryCatch(
                eigen(H, symmetric = TRUE),
                error = function(e) NULL
            )
            if (is.null(eigs)) {
                cat("Could not compute eigenvalues of Hessian; skipping Hessian check.\n\n")
            } else {
                max_eval <- max(abs(eigs$values))
                bad <- which(eigs$values / max_eval <= eval_eps)
                if (length(bad) > 0) {
                    model_OK <- FALSE
                    cat("Non-positive definite (NPD) Hessian\n\n")
                    prt_explain(
                        "The Hessian matrix represents the curvature of the ",
                        "log-likelihood surface at the MLE of the parameters. ",
                        "A non-positive-definite Hessian means that the likelihood ",
                        "surface is approximately flat (or upward-curving) in some ",
                        "directions, indicating the model is overfitted or poorly posed."
                    )
                    ## assign parameter names to rows/cols of the Hessian
                    ## the Hessian covers theta (RE) params and, if available,
                    ## also the fixed-effect (beta) params
                    h_names <- if (nrow(H) > n_theta) {
                        c(names(theta), nn)
                    } else {
                        names(theta)
                    }
                    if (length(h_names) == nrow(H)) {
                        cat(sprintf(
                            "maximum Hessian eigenvalue = %1.3g\n",
                            eigs$values[1]
                        ))
                        for (b in bad) {
                            cat(sprintf(
                                "Hessian eigenvalue %d = %1.3g (relative val = %1.3g)\n",
                                b, eigs$values[b], eigs$values[b] / eigs$values[1]
                            ))
                            bad_vec <- eigs$vectors[, b]
                            bad_elements <- which(abs(bad_vec) > evec_eps)
                            if (length(bad_elements) > 0) {
                                cat("   bad elements:",
                                    paste(h_names[bad_elements], collapse = ", "), "\n")
                            }
                        }
                        cat("\n")
                    }
                }
            }
        }
    }

    if (model_OK) cat("model looks OK!\n")
    return(invisible(model_OK))
}
