cc_dstats_supp <- function(df1,
                           ccVar,
                           denom_case,
                           denom_ctrl,
                           signif_digits = 3,
                           include_var_meta = FALSE,
                           var_meta_df = NULL,
                           var_meta_merge_by = NULL) {
  # Description:
  #   Make a descriptive stats table.
  # Args:
  #   df1: Data frame with variables to include in table.
  #   ccVar: Case_control variable.
  #   denom_case: Case var name.
  #   denom_ctrl: Ctrl var name.
  #   signif_digits: Number of significant digits to return.
  #   include_var_meta: Add variable metadata (TRUE/FALSE).
  #   var_meta_df: Data frame with variable metadata.
  #   var_meta_merge_by: String with var name to merge meta and summary data by.
  # Dependencies:
  check_packages(
    bioc_packages = c(""),
    cran_packages = c("")
  )
  # Returns:
  #   tableS1: Table with descriptive stats by case-control status.

  # Split data
  if (include_var_meta == TRUE) {
    df1 <- df1[c(ccVar, var_meta_df[var_meta_merge_by][, 1])]
  }

  df1_split <- split(df1, df1[ccVar])

  df1_case <- df1_split[[denom_case]]
  df1_ctrl <- df1_split[[denom_ctrl]]

  # Filter numerical vars
  df1_list <- list(
    df1 = Filter(is.numeric, df1),
    df1_case = Filter(is.numeric, df1_case),
    df1_ctrl = Filter(is.numeric, df1_ctrl)
  )

  # Apply summary stats
  df1_summary <- lapply(
    df1_list,
    function(x) {
      y <- data.frame(
        n = apply(x, 2, function(x2) {
          length(na.omit(x2))
        }),
        mean = apply(x, 2, mean, na.rm = TRUE),
        sd = apply(x, 2, sd, na.rm = TRUE),
        median = apply(x, 2, function(x2) {
          summary(x2)["Median"]
        }),
        IQR1 = apply(x, 2, function(x2) {
          summary(x2)["1st Qu."]
        }),
        IQR3 = apply(x, 2, function(x2) {
          summary(x2)["3rd Qu."]
        }),
        min = apply(x, 2, function(x2) {
          summary(x2)["Min."]
        }),
        max = apply(x, 2, function(x2) {
          summary(x2)["Max."]
        })
      )
      y <- sapply(y, function(x3) {
        signif(x3, digits = signif_digits)
      })
      return(y)
    }
  )

  df1_re_merged <- data.frame(
    name1 = names(df1_list$df1),
    df1_summary[[1]],
    df1_summary[[2]],
    df1_summary[[3]]
  )

  # Merge with var metadata
  if (include_var_meta == TRUE) {
    t_all <- merge(
      x = var_meta_df,
      y = df1_re_merged,
      by.x = var_meta_merge_by,
      by.y = var_meta_merge_by,
      all.y = TRUE,
      sort = FALSE
    )

    # Add top layer
    head1 <- c(
      rep("", ncol(var_meta_df)),
      "Full cohort",
      rep("", ncol(df1_summary[[1]]) - 1),
      denom_case,
      rep("", ncol(df1_summary[[1]]) - 1),
      denom_ctrl,
      rep("", ncol(df1_summary[[1]]) - 1)
    )

    tableS1 <- rbind(head1, colnames(t_all), t_all)
  } else {
    tableS1 <- df1_re_merged
  }

  # Return table
  return(tableS1)
}
