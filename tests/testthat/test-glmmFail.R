data("sleepstudy", package = "lme4")

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gSim(), a general simulation function ...

set.seed(101)
dBc <- gSim(family=binomial(link="cloglog"), nbinom = 1) # {0,1} Binomial

## m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##             family = binomial, data = cbpp)
#context("Errors and warnings from glmer")
test_that("glmer", {
    expect_error(glmer(y ~ 1 + (1|block), data=dBc, family=binomial(link="cloglog")),
                 "Response is constant")
    expect_error(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                       family = binomial, data = cbpp, REML=TRUE),
                   "unused argument.*REML")
    ## glmer(family=gaussian) should now return glmerMod (using full GLMM machinery)
    m3 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m5 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy, family=gaussian)
    expect_is(m3, "glmerMod")
    expect_is(m5, "glmerMod")
    ## fixed effects should be approximately equal to lmer results
    expect_equal(fixef(m3), fixef(m4), tolerance = 1e-4)
    expect_equal(fixef(m3), fixef(m5))
})

