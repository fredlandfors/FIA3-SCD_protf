cc_dstats_main <- function(df1,
                           ccVar,
                           denom_case,
                           denom_ctrl,
                           n_digits,
                           use_var_metadata = FALSE,
                           old_names = NULL,
                           new_names = NULL,
                           new_unit = NULL,
                           p_vals = FALSE) {
  # Description:
  #   Make a table with descriptive stats
  # Args:
  #   df1: Data frame with variables to include in table.
  #   ccVar: Case_control variable.
  #   denom_case: Case var name.
  #   denom_ctrl: Ctrl var name.
  #   n_digits: Number of digits.
  #   use_var_metadata: True/false.
  #   old_names: Vector of names in input df.
  #   new_names: Vector of names for table output.
  #   new_unit: Vector of new units to use in output df.
  #   p_vals: Return P-values (true/false.)
  # Dependencies:
  check_packages(
    bioc_packages = c(""),
    cran_packages = c("")
  )
  # Returns:
  #   table1: Table with descriptive stats by case-control status.

  #----------------------------------------------------------------------------
  # Split data
  #----------------------------------------------------------------------------

  df1_split <- split(df1, df1[ccVar])

  df1_case <- df1_split[[denom_case]]
  df1_ctrl <- df1_split[[denom_ctrl]]

  #----------------------------------------------------------------------------
  # Column 1: variable names column
  #----------------------------------------------------------------------------
  varName <- names(df1)

  if (use_var_metadata == TRUE) {
    df_cha <- data.frame(
      old_n = old_names,
      new_n = new_names,
      new_u = new_unit
    )
  }

  # Change var names
  for (i in 1:length(varName)) {
    # Numeric
    if (is.numeric(df1_case[, i]) == TRUE) {
      if (use_var_metadata == TRUE) {
        names_df <- subset(df_cha, old_n == varName[i])
        varName[i] <- paste0(names_df$new_n, " (", names_df$new_u, ")")
      }
    }

    # Factors
    if (is.factor(df1_case[, i]) == TRUE) {
      if (use_var_metadata == TRUE) {
        names_df <- subset(df_cha, old_n == varName[i])
        varName[i] <- paste0(names_df$new_n)
      }

      # Add factor levels to name
      fct_names <- levels(na.omit(df1[, i]))

      varName[i] <- paste0(
        varName[i],
        "&  ",
        paste0(fct_names, collapse = "&  ")
      )
    }
  }


  #----------------------------------------------------------------------------
  # Column 2: n CASE
  #----------------------------------------------------------------------------
  nCase <- seq(1, ncol(df1_case))
  for (i in 1:length(nCase)) {
    nCase[i] <- length(na.omit(df1_case[, i]))
  }
  #----------------------------------------------------------------------------
  # Column 3: mean+-SD or n+% var, CASE
  #----------------------------------------------------------------------------
  case_mean_or_percent <- rep("", ncol(df1_case))
  for (i in 1:length(case_mean_or_percent)) {
    # Numeric
    if (is.numeric(df1_case[, i]) == TRUE & nCase[i] > 0) {
      cc_M <- signif(mean(na.omit(df1_case[, i])), n_digits)
      cc_SD <- signif(sd(na.omit(df1_case[, i])), n_digits)
      case_mean_or_percent[i] <- paste(cc_M, " (SD: \u00B1 ", cc_SD, ")", sep = "")
    }

    # Factors
    if (is.factor(df1_case[, i]) == TRUE) {
      # Calc n fct levels
      fct_levels <- length(levels(na.omit(df1[, i])))
      fct_names <- levels(na.omit(df1[, i]))

      # .. with >= 2 levels
      if (fct_levels >= 2) {
        calcPerc <- as.numeric(rep(0, fct_levels))
        calcN <- as.numeric(rep(0, fct_levels))
        stringPerc <- as.character(rep("filler", fct_levels))

        for (m in seq_len(fct_levels)) {
          calcN[m] <- table(df1_case[, i])[m]
          calcPerc[m] <- 100 * (table(df1_case[, i])[m] / length(na.omit(df1_case[, i])))
          calcPerc[m] <- round(calcPerc[m], 1)
        }

        for (m in seq_len(fct_levels)) {
          stringPerc[m] <- paste0(" &", calcN[m], " (", calcPerc[m], " %)")
        }

        case_mean_or_percent[i] <- paste0(stringPerc, collapse = "")
      }
    }
  }

  #----------------------------------------------------------------------------
  # Column 4: median+IQR+range, CASE
  #----------------------------------------------------------------------------
  case_median <- rep("", ncol(df1_case))
  for (i in 1:length(case_median)) {
    # Numeric
    if (is.numeric(df1_case[, i]) == TRUE & nCase[i] > 0) {
      cc_median <- signif(summary(df1_case[, i])["Median"], n_digits)
      cc_iqr1 <- signif(summary(df1_case[, i])["1st Qu."], n_digits)
      cc_iqr3 <- signif(summary(df1_case[, i])["3rd Qu."], n_digits)
      cc_min <- signif(summary(df1_case[, i])["Min."], n_digits)
      cc_max <- signif(summary(df1_case[, i])["Max."], n_digits)
      case_median[i] <- paste0(
        cc_median, " (IQR: ", cc_iqr1, " - ", cc_iqr3,
        "; range: ", cc_min, " - ", cc_max, ")"
      )
    }
  }

  #----------------------------------------------------------------------------
  # Column 5: n CTRL
  #----------------------------------------------------------------------------
  nCtrl <- seq(1, ncol(df1_ctrl))
  for (i in 1:length(nCtrl)) {
    nCtrl[i] <- length(na.omit(df1_ctrl[, i]))
  }

  #----------------------------------------------------------------------------
  # Column 6: mean+-SD or n+% var, CTRL
  #----------------------------------------------------------------------------
  ctrl_mean_or_percent <- rep("", ncol(df1_ctrl))
  for (i in 1:length(ctrl_mean_or_percent)) {
    # Numeric
    if (is.numeric(df1_ctrl[, i]) == TRUE & nCtrl[i] > 0) {
      cc_M <- signif(mean(na.omit(df1_ctrl[, i])), n_digits)
      cc_SD <- signif(sd(na.omit(df1_ctrl[, i])), n_digits)
      ctrl_mean_or_percent[i] <- paste(cc_M, " (SD: \u00B1 ", cc_SD, ")", sep = "")
    }

    # Factors
    if (is.factor(df1_ctrl[, i]) == TRUE) {
      # Calc n fct levels
      fct_levels <- length(levels(na.omit(df1[, i])))
      fct_names <- levels(na.omit(df1[, i]))

      # With >= 2 levels
      if (fct_levels >= 2) {
        calcPerc <- as.numeric(rep(0, fct_levels))
        calcN <- as.numeric(rep(0, fct_levels))
        stringPerc <- as.character(rep("filler", fct_levels))

        for (m in seq_len(fct_levels)) {
          calcN[m] <- table(df1_ctrl[, i])[m]
          calcPerc[m] <- 100 * (table(df1_ctrl[, i])[m] / length(na.omit(df1_ctrl[, i])))
          calcPerc[m] <- round(calcPerc[m], 1)
        }

        for (m in seq_len(fct_levels)) {
          stringPerc[m] <- paste0(" &", calcN[m], " (", calcPerc[m], " %)")
        }
      }

      ctrl_mean_or_percent[i] <- paste0(stringPerc, collapse = "")
    }
  }


  #----------------------------------------------------------------------------
  # Column 7: median+IQR+range, CTRL
  #----------------------------------------------------------------------------
  ctrl_median <- rep("", ncol(df1_ctrl))
  for (i in 1:length(ctrl_median)) {
    # Numeric
    if (is.numeric(df1_ctrl[, i]) == TRUE & nCtrl[i] > 0) {
      cc_median <- signif(summary(df1_ctrl[, i])["Median"], n_digits)
      cc_iqr1 <- signif(summary(df1_ctrl[, i])["1st Qu."], n_digits)
      cc_iqr3 <- signif(summary(df1_ctrl[, i])["3rd Qu."], n_digits)
      cc_min <- signif(summary(df1_ctrl[, i])["Min."], n_digits)
      cc_max <- signif(summary(df1_ctrl[, i])["Max."], n_digits)
      ctrl_median[i] <- paste0(
        cc_median, " (IQR: ", cc_iqr1, " - ", cc_iqr3,
        "; range: ", cc_min, " - ", cc_max, ")"
      )
    }
  }

  #----------------------------------------------------------------------------
  ###### P-value
  #----------------------------------------------------------------------------
  pVal <- rep(NA, ncol(df1))
  pVal_mann <- rep(NA, ncol(df1))
  idVar <- names(df1)

  for (i in 1:length(pVal)) {

    # Numeric
    if (is.numeric(df1[, i]) == TRUE &&
      names(df1[i]) != names(df1[ccVar])) {
      if (length(unique(na.omit(df1[df1[, ccVar] == denom_case, i]))) > 0 &&
        length(unique(na.omit(df1[df1[, ccVar] == denom_ctrl, i]))) > 0) {
        formula1 <- formula(paste(idVar[i], " ~ ", ccVar, sep = ""))
        test1 <- t.test(formula1, data = df1)
        mann1 <- wilcox.test(formula1, data = df1, alternative = "two.sided")
        pVal[i] <- signif(test1$p.value, n_digits)
        pVal_mann[i] <- signif(mann1$p.value, n_digits)
      }
    }

    # Factor
    if (is.factor(df1[, i]) == TRUE && names(df1[i]) != names(df1[ccVar])) {
      # .. with 2 levels
      if (length(unique(na.omit(df1[df1[, ccVar] == denom_case, i]))) == 2 &&
        length(unique(na.omit(df1[df1[, ccVar] == denom_ctrl, i]))) == 2) {
        test1 <- chisq.test(x = df1[, i], y = df1[, ccVar])
        pVal[i] <- signif(test1$p.value, n_digits)
      }
    }
  }

  #----------------------------------------------------------------------------
  # Make table
  #----------------------------------------------------------------------------
  table1 <- data.frame(
    variable = varName,
    cases = nCase,
    mean1 = case_mean_or_percent,
    median1 = case_median,
    ctrls = nCtrl,
    mean2 = ctrl_mean_or_percent,
    median2 = ctrl_median,
    p = pVal,
    mann = pVal_mann,
    stringsAsFactors = FALSE
  )

  table1_append <- data.frame()

  for (i in seq_len(nrow(table1))) {
    if (is.factor(df1[, i]) == TRUE) {
      new_vars <- strsplit(table1[i, "variable"], "&")[[1]][-1]
      new_mean1 <- strsplit(table1[i, "mean1"], "&")[[1]][-1]
      new_mean2 <- strsplit(table1[i, "mean2"], "&")[[1]][-1]
      n_fact <- length(new_vars)
      n_appended <- nrow(table1_append)

      append1 <- data.frame(
        index = seq(i + n_appended + 1, i + n_fact + n_appended),
        variable = new_vars,
        cases = rep("", n_fact),
        mean1 = new_mean1,
        median1 = rep("", n_fact),
        ctrls = rep("", n_fact),
        mean2 = new_mean2,
        median2 = rep("", n_fact),
        p = rep("", n_fact),
        mann = rep("", n_fact),
        stringsAsFactors = FALSE
      )

      table1_append <- rbind(table1_append, append1)
    }
  }

  insertRow <- function(existingDF, newrow, r) {
    existingDF[seq(r + 1, nrow(existingDF) + 1), ] <- existingDF[seq(r, nrow(existingDF)), ]
    existingDF[r, ] <- newrow
    existingDF
  }

  table1_fin <- table1

  for (i in seq_len(nrow(table1_append))) {
    table1_fin <- insertRow(
      table1_fin,
      table1_append[i, -1],
      table1_append[i, "index"]
    )
  }

  table1_fin2 <- vapply(
    table1_fin,
    function(x) {
      for (i in seq_len(nrow(table1_fin))) {
        x[i] <- ifelse(
          grepl("&", x[i]),
          strsplit(x[i], split = "&")[[1]][1],
          x[i]
        )
      }
      return(x)
    },
    FUN.VALUE = character(nrow(table1_fin))
  )

  table1_fin2[is.na(table1_fin2)] <- ""

  table1_fin2 <- data.frame(table1_fin2, stringsAsFactors = FALSE)

  # Return P-values (yes/no)
  if (p_vals == FALSE) {
    table1_fin2 <- table1_fin2[1:7]
  }

  #----------------------------------------------------------------------------
  # Final return
  #----------------------------------------------------------------------------
  return(table1_fin2)
}
