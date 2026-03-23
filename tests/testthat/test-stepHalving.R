load(system.file("testdata","survdat_reduced.Rda",package="lme4"))

test_that('Step-halving works properly', {
  # this example is known to require step-halving (or at least has in the past
  # required step-halving)
  form <- survprop~(1|nobs)
  m <- glmer(form,weights=eggs,data=survdat_reduced,family=binomial,nAGQ=1L)
  expect_that(m, is_a("glmerMod"))
})

test_that('PIRLS step-halving NaN recovery for binomial(link="log")', {
  # Regression test for GitHub issue: "PIRLS Step-Halving Failure"
  # binomial(link="log") can cause mu > 1 during PIRLS, which produces NaN
  # working weights.  NaN working weights lead to NaN in the Cholesky
  # factorisation (Inf/Inf = NaN in off-diagonal elements), and hence NaN
  # delu/delb.  The step-halving formula (olddelu + NaN)/2 = NaN, so naive
  # halving never recovers.  The fix resets delu to the previous valid
  # values before attempting further halvings.
  df.odd <- data.frame(
    Subject     = rep(seq(1, 24, 2), each = 4),
    Condition   = c("A", "B", "A", "B"),
    Block       = c(1, 2, 3, 4),
    # High Yes counts push the log-link predictor close to 0,
    # making it easy for the optimiser to temporarily explore mu > 1.
    Yes = c(9, 10, 8, 10, 10, 9, 9, 10,  8,  9, 10, 9,
            10,  8, 9, 10,  9, 8, 10, 9, 10, 10,  9, 8)
  )
  df.odd$No <- 10L - df.odd$Yes

  df.even <- data.frame(
    Subject   = rep(seq(2, 24, 2), each = 4),
    Condition = c("A", "B", "A", "B"),
    Block     = c(1, 2, 3, 4),
    Yes = c(8, 9, 10, 9, 10, 8,  9,  8, 10,  9, 8, 10,
            9, 10, 8, 9,  8, 10, 9, 10,  9,  8, 10,  9)
  )
  df.even$No <- 10L - df.even$Yes

  df <- rbind(df.odd, df.even)
  df$Counterbalance <- df$Subject %% 2

  # This used to throw "Error: PIRLS loop resulted in NaN value"
  m <- suppressWarnings(
    glmer(cbind(Yes, No) ~ Condition + Block + Counterbalance + (1 | Subject),
          family = binomial(link = "log"),
          data = df)
  )
  expect_is(m, "glmerMod")
})
