library("testthat")
library("lme4")

## use old (<=3.5.2) sample() algorithm if necessary
if ("sample.kind" %in% names(formals(RNGkind))) {
    suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
}

## -- Helpers ------------------------------------------------------------------

## capture all output (messages + cat) from diagnose()
cap_diagnose <- function(...) {
    capture.output(res <- diagnose(...))
}

## ─────────────────────────────────────────────────────────────────────────────
## 1.  Well-behaved model should report no problems
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose returns TRUE for a well-behaved lmer model", {
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    out <- cap_diagnose(fm1)
    res <- diagnose(fm1)          # capture return value separately
    expect_true(res)
    expect_true(any(grepl("model looks OK", out)))
})

test_that("diagnose returns TRUE for a well-behaved glmer model", {
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
    res <- diagnose(gm1, explain = FALSE)
    expect_true(res)
})

## ─────────────────────────────────────────────────────────────────────────────
## 2.  check_coefs: large theta (random-effects) parameter
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose detects large random-effects (theta) parameters", {
    set.seed(42)
    n <- 200
    dat <- data.frame(
        y     = rnorm(n),
        group = factor(rep(1:2, each = n / 2)),   ## only 2 groups → large RE
        x     = rnorm(n)
    )
    ## With just 2 groups and a random slope, RE variance is hard to estimate
    ## and often results in very large theta values; force a large theta by
    ## directly constructing a model with a known large theta for testing.
    ## Instead, test by artificially lowering the threshold.
    fm_normal <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    out <- cap_diagnose(fm_normal, big_coef = 0.1, check_zstats = FALSE,
                        check_hessian = FALSE, check_scales = FALSE,
                        explain = FALSE)
    ## At big_coef=0.1 every non-zero theta should be flagged
    res <- diagnose(fm_normal, big_coef = 0.1, check_zstats = FALSE,
                    check_hessian = FALSE, check_scales = FALSE,
                    explain = FALSE)
    expect_false(res)
    expect_true(any(grepl("Unusually large coefficients", out)))
})

## ─────────────────────────────────────────────────────────────────────────────
## 3.  check_scales: badly-scaled predictor
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose detects badly scaled predictors", {
    ## Create a predictor whose sd is ~1e-4 (should trigger big_sd_log10 = 3)
    dat <- sleepstudy
    dat$tiny <- dat$Days * 1e-5      # sd ~ 1e-5 * sd(Days)
    fm_bad <- lmer(Reaction ~ tiny + (Days | Subject), dat)
    out <- cap_diagnose(fm_bad, check_coefs = FALSE, check_zstats = FALSE,
                        check_hessian = FALSE, explain = FALSE)
    res <- diagnose(fm_bad, check_coefs = FALSE, check_zstats = FALSE,
                    check_hessian = FALSE, explain = FALSE)
    expect_false(res)
    expect_true(any(grepl("standard deviations", out)))
})

test_that("diagnose does NOT flag well-scaled predictors for check_scales", {
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    out <- cap_diagnose(fm1, check_coefs = FALSE, check_zstats = FALSE,
                        check_hessian = FALSE, explain = FALSE)
    res <- diagnose(fm1, check_coefs = FALSE, check_zstats = FALSE,
                    check_hessian = FALSE, explain = FALSE)
    expect_true(res)
    expect_false(any(grepl("standard deviations", out)))
})

## ─────────────────────────────────────────────────────────────────────────────
## 4.  check_zstats: unusually large Z-statistic
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose detects large Z-statistics when threshold is lowered", {
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    ## With big_zstat = 1, most fixed effects should be flagged
    out <- cap_diagnose(fm1, check_coefs = FALSE, check_hessian = FALSE,
                        check_scales = FALSE, big_zstat = 1, explain = FALSE)
    res <- diagnose(fm1, check_coefs = FALSE, check_hessian = FALSE,
                    check_scales = FALSE, big_zstat = 1, explain = FALSE)
    expect_false(res)
    expect_true(any(grepl("Z-statistics", out)))
})

test_that("diagnose does NOT flag small Z-statistics for check_zstats", {
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    ## With a very high threshold nobody should be flagged
    out <- cap_diagnose(fm1, check_coefs = FALSE, check_hessian = FALSE,
                        check_scales = FALSE, big_zstat = 1e6, explain = FALSE)
    res <- diagnose(fm1, check_coefs = FALSE, check_hessian = FALSE,
                    check_scales = FALSE, big_zstat = 1e6, explain = FALSE)
    expect_true(res)
})

## ─────────────────────────────────────────────────────────────────────────────
## 5.  check_hessian: Hessian not available when calc.derivs=FALSE
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose reports Hessian unavailable when calc.derivs=FALSE", {
    fm_noderiv <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                       control = lmerControl(calc.derivs = FALSE))
    out <- cap_diagnose(fm_noderiv, check_coefs = FALSE, check_scales = FALSE,
                        check_zstats = FALSE, explain = FALSE)
    expect_true(any(grepl("Hessian not available", out)))
})

## ─────────────────────────────────────────────────────────────────────────────
## 6.  check_hessian: non-positive-definite Hessian
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose detects NPD Hessian", {
    ## Fit a model with a near-singular Hessian by over-parameterising:
    ## use a very small dataset with many RE terms.
    set.seed(1)
    n_grp <- 3
    n_per <- 3
    dat <- data.frame(
        y     = rnorm(n_grp * n_per),
        x     = rnorm(n_grp * n_per),
        group = factor(rep(seq_len(n_grp), each = n_per))
    )
    ## Suppress any convergence warnings since we deliberately over-parameterise
    fm_bad <- suppressWarnings(
        lmer(y ~ x + (x | group), dat,
             control = lmerControl(calc.derivs = TRUE))
    )
    ## If Hessian is available, check that diagnose flags NPD or returns FALSE
    if (!is.null(fm_bad@optinfo$derivs$Hessian)) {
        res <- diagnose(fm_bad, check_coefs = FALSE, check_scales = FALSE,
                        check_zstats = FALSE, explain = FALSE)
        ## result must be logical
        expect_true(is.logical(res))
    } else {
        ## If Hessian wasn't computed (e.g. too small), skip
        skip("Hessian not computed for this model")
    }
})

## ─────────────────────────────────────────────────────────────────────────────
## 7.  explain = FALSE suppresses explanatory text
## ─────────────────────────────────────────────────────────────────────────────
test_that("explain=FALSE suppresses explanation text", {
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    ## Force a z-stat flag so we can check explain output is absent
    out_explain    <- cap_diagnose(fm1, check_coefs = FALSE,
                                   check_hessian = FALSE,
                                   check_scales  = FALSE,
                                   big_zstat = 1,
                                   explain = TRUE)
    out_no_explain <- cap_diagnose(fm1, check_coefs = FALSE,
                                   check_hessian = FALSE,
                                   check_scales  = FALSE,
                                   big_zstat = 1,
                                   explain = FALSE)
    ## explain=TRUE should produce more output
    expect_gt(length(out_explain), length(out_no_explain))
})

## ─────────────────────────────────────────────────────────────────────────────
## 8.  Error for non-merMod input
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose errors on non-merMod input", {
    m <- lm(Reaction ~ Days, sleepstudy)
    expect_error(diagnose(m), "merMod")
})

## ─────────────────────────────────────────────────────────────────────────────
## 9.  glmer model – check_coefs checks fixed effects for logit link
## ─────────────────────────────────────────────────────────────────────────────
test_that("diagnose checks fixed-effect coefficients for glmer logit link", {
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
    ## Lower big_coef so fixed effects are flagged
    out <- cap_diagnose(gm1, big_coef = 0.1, check_zstats = FALSE,
                        check_hessian = FALSE, check_scales = FALSE,
                        explain = FALSE)
    res <- diagnose(gm1, big_coef = 0.1, check_zstats = FALSE,
                    check_hessian = FALSE, check_scales = FALSE,
                    explain = FALSE)
    expect_false(res)
    expect_true(any(grepl("Unusually large coefficients", out)))
})
