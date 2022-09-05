mtrx_glm <- function(df,
                     col_vars,
                     row_vars_all,
                     adj_vars) {
  # Description:
  #   Function for standardized beta/logistic regression matrix, depending on
  #   whether y is linear or binary. Scales variables.
  # Args:
  #   df: Data frame with variables to include in table.
  #   col_vars: columns variable.
  #   row_vars_all: row variables.
  #   adj_vars: adjustment variables.
  # Dependencies:
  check_packages(
    bioc_packages = c(""),
    cran_packages = c("")
  )
  # Returns:
  #   glm_est_matrix: Matrix of beta/OR coefficients.
  
  row_vars <- names(df[names(df) %in% row_vars_all])
  df <- df[c(col_vars, row_vars, adj_vars)]
  
  df <- data.frame(rapply(df, scale, classes = "numeric", how = "replace"))
  
  glm_est_matrix <- matrix(ncol = length(col_vars), nrow = length(row_vars))
  colnames(glm_est_matrix) <- col_vars
  rownames(glm_est_matrix) <- row_vars
  
  glm_p_matrix <- matrix(ncol = length(col_vars), nrow = length(row_vars))
  colnames(glm_p_matrix) <- col_vars
  rownames(glm_p_matrix) <- row_vars
  
  for (j in seq(1, length(col_vars))) {
    if (is.numeric(df[, j])) {
      for (i in seq(1, length(row_vars))) {
        formula_ij <- as.formula(paste0(
          colnames(glm_est_matrix)[j], " ~ ", rownames(glm_est_matrix)[i], " + ",
          paste(adj_vars, collapse = " + ")
        ))
        fit_ij <- lm(formula = formula_ij, data = df)
        sum_ij <- summary(fit_ij)
        glm_est_matrix[i, j] <- sum_ij$coefficients[2 , "Estimate"]
        glm_p_matrix[i, j] <- sum_ij$coefficients[2 , "Pr(>|t|)"]
      }
    }
    if (is.factor(df[, j]) & length(levels(df[, j])) == 2) {
      for (i in seq(1, length(row_vars))) {
        formula_ij <- as.formula(paste0(
          colnames(glm_est_matrix)[j], " ~ ", rownames(glm_est_matrix)[i], " + ",
          paste(adj_vars, collapse = " + ")
        ))
        fit_ij <- glm(formula = formula_ij, data = df,
                      family = binomial(link = "logit"))
        sum_ij <- summary(fit_ij)
        glm_est_matrix[i, j] <- exp(sum_ij$coefficients[2 , "Estimate"])
        glm_p_matrix[i, j] <- sum_ij$coefficients[2 , "Pr(>|z|)"]
      }
    }
  }
  glm_est_matrix <- round(glm_est_matrix, 2)
  coef <- data.frame(glm_est_matrix, stringsAsFactors = TRUE)
  p_vals <- data.frame(glm_p_matrix, stringsAsFactors = TRUE)
  glm_list <- list(coef = coef, p_vals = p_vals)
  return(glm_list)
}
