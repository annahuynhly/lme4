if (.Platform$OS.type != "windows") {
    library(lme4)
    library(testthat)

    ## glmer(*, gaussian) now uses full GLMM machinery and returns glmerMod
    m3 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m5 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                                 family=gaussian)
    expect_is(m3, "glmerMod")
    expect_is(m5, "glmerMod")
    expect_equal(fixef(m3),fixef(m5))
    ## fixed effects should be approximately equal to lmer results
    expect_equal(fixef(m3), fixef(m4), tolerance = 1e-4)
    expect_equal(VarCorr(m3), VarCorr(m5))
    expect_equal(getME(m3,"theta"), getME(m5,"theta"))

    ## glmer() - poly() + interaction
    if (requireNamespace("mlmRev")) withAutoprint({
        data(Contraception, package="mlmRev")
        ## ch := with child
        Contraception <- within(Contraception, ch <- livch != "0")
        ## gmC1 <- glmer(use ~ poly(age,2) + ch + age:ch + urban + (1|district),
        ##               Contraception, binomial)
### not a 'warning' per se {cannot suppressWarnings(.)}:
###    fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
### also printed with print(): labeled as  "fit warnings"

        ## ==> from ../R/modular.R  chkRank.drop.cols()
        ## --> Use   control = glmerControl(check.rankX = "ignore+drop.cols"))
        ## because further investigation shows "the problem" is really already
        ##     in model.matrix():
        set.seed(101)
        dd <- data.frame(ch = c("Y","N")[1+rbinom(12, 1, 0.7)], age = rlnorm(12, 16))
        colnames(mm1 <- model.matrix( ~ poly(age,2) + ch + age:ch, dd))
        ## "(Int.)" "poly(age, 2)1" "poly(age, 2)2" "chY" "chN:age" "chY:age"      
        ## If we make the poly() columns to regular variables, can interact:
        d2 <- within(dd, { p2 <- poly(age,2); ageL <- p2[,1]; ageQ <- p2[,2]; rm(p2)})
        ## then, we can easily get what want
        (mm2 <- model.matrix( ~ ageL+ageQ + ch + ageL:ch, d2))
        ## actually even more compactly now ("drawback": 'ageQ' at end):
        (mm2. <- model.matrix( ~ ageL*ch + ageQ, d2))
        cn2 <- colnames(mm2)
        stopifnot(identical(mm2[,cn2], mm2.[,cn2]))
    })
} ## skip on windows (for speed)
