cc_multiclogit <- function(df1,
                           ccVar,
                           adj,
                           stratum,
                           x_vars,
                           ci_adj = 1,
                           scaling = TRUE) {
  # Description:
  #   Fits clogit models and returns a table with results.
  # Args:
  #   df1: Data frame.
  #   ccVar: String containing name of case_control variable.
  #   adj:  Vector containing names of variables in df1 to adjust for.
  #   stratum:  String containing name of strata variable used for matching.
  #   x_vars: Vector containing names of x-vars.
  #   scaling: scale variables before running regression (TRUE/FALSE).
  #   ci_adj: Correct 95 # CI:s for multiple comparison.
  # Dependencies:
  check_packages(
    bioc_packages = c(""),
    cran_packages = c("survival")
  )
  library(survival)
  # Returns:
  #   table1: Table with clogit stats.

  # Scale numerical vars
  if (scaling == TRUE) {
    df1 <- rapply(
      df1,
      function(x) {
        scale(x)
      },
      classes = "numeric",
      how = "replace"
    )
  }

  # Convert cc var to numeric (case = 1, ctrl = 0)
  df1[, ccVar] <- as.numeric(df1[, ccVar]) - 1

  df1_x <- names(df1[x_vars])
  n <- length(df1_x)

  n_incl <- rep(NA, n)
  logit_or <- rep(NA, n)
  ci_low <- rep(NA, n)
  ci_high <- rep(NA, n)
  logit_p <- rep(NA, n)

  # Run separate models for x-vars
  for (j in seq(1, n)) {
    form_c <- as.formula(
      paste(ccVar, " ~ ", df1_x[j], adj, "+ strata(", stratum, ")", sep = "")
    )
    fit_c <- survival::clogit(form_c, df1)
    sum_c <- summary(fit_c)

    n_incl[j] <- sum_c$n
    logit_or[j] <- round(sum_c$coefficients[1, 2], 2)
    ci_low[j] <- round(exp(confint(fit_c, level = 1 - 0.05/ci_adj)[1, 1]), 2)
    ci_high[j] <- round(exp(confint(fit_c, level = 1 - 0.05/ci_adj)[1, 2]), 2)
    logit_p[j] <- sum_c$coefficients[1, 5]
  }

  # Make table
  table1 <- data.frame(
    name = as.character(df1_x),
    n = as.numeric(n_incl),
    OR = as.numeric(logit_or),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    p = as.numeric(logit_p)
  )
  table1$name <- as.character(table1$name)

  # Return data
  return(table1)
}
