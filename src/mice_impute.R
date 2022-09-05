impute_by_MICE <- function(df,
                           vars_include,
                           method1,
                           set_seed = 123,
                           scaling = FALSE) {
  # Description:
  #   Performs missing values imputation on clinical data using MICE.
  # Args:
  #   df: Data frame.
  #   vars_include: Variables to include in MICE model.
  #   method: Method to impute by.
  #   scaling: If true, performs scaling before imputation.
  # Dependencies:
  check_packages(
    bioc_packages = c(""),
    cran_packages = c( "naniar", "mice", "ggplot2")
  )
  # Returns:
  #   return_list: List of objects.
  # Constants:
  # Error handling:
  # Define lists:

  return_list <- list()
  return_list$data <- list()
  return_list$functions <- list()
  return_list$plots <- list()

  # Define functions:

  return_list$functions$show_missing <- function(x) {
    get_perc <- function(z) {
      sum(is.na(z)) / length(z) * 100
    }
    y <- sort(apply(x, 2, get_perc), decreasing = FALSE)
    return(y)
  }

  #----------------------------------------------------------------------------
  # Scale data
  #----------------------------------------------------------------------------
  if (scaling == TRUE) {
    df <- scale(Filter(is.numeric, df))
  }

  #----------------------------------------------------------------------------
  # Missing values analysis
  #----------------------------------------------------------------------------
  return_list$data$data_import <- df
  return_list$plots$vismissplot <- naniar::vis_miss(df)
  return_list$plots$upsetplot <- naniar::gg_miss_upset(df)

  return_list$data$percent_missing <- return_list$functions$show_missing(df)

  #----------------------------------------------------------------------------
  # Perform MICE imputation on clinical data
  #----------------------------------------------------------------------------
  mice_data_in <- df[vars_include]
  
  return_list$data$mice_data_in <- mice_data_in

  mice_fit <- mice::mice(
    mice_data_in,
    m = 5,
    maxit = 50,
    method = method1,
    seed = set_seed
  )
  return_list$data$mice_fit <- mice_fit

  return_list$plots$mice_meansd <- plot(mice_fit, layout = c(2, 9))
  return_list$plots$mice_density <- mice::densityplot(mice_fit)
  return_list$plots$mice_stripplot <- mice::stripplot(mice_fit)

  mice_out <- mice::complete(mice_fit)
  
  clinical_imputed <- data.frame(
    df[!names(df) %in% names(mice_out)],
    mice_out
  )
  
  return_list$data$mice_data_out <- clinical_imputed
  
  #----------------------------------------------------------------------------
  # FIN RETURN
  #----------------------------------------------------------------------------
  return(return_list)
  
}
