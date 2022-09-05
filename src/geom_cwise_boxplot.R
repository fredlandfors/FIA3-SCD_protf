geom_cwise_boxplot <- function(x,
                              c,
                              id = as.character(seq(1, nrow(x)))) {
  # Description:
  #   Plots a column-wise boxplot
  # Args:
  #   x: PCA x numerical data.
  #   c: Variable to color by (obs! as df with one column).
  #   id: Sample ID var (obs! as df with one column).
  # Dependencies:
  check_packages(
    cran_packages = c("ggplot2")
  )
  # Returns:
  #   p2: The boxplot
  
  # Make data frame
  gboxplot_data <- data.frame(
    id_var = id,
    colour_by = c,
    x
  )

  gboxplot_data <- rapply(
    gboxplot_data,
    scale,
    classes = "numeric",
    how = "replace"
  )

  long_data <- reshape(
    gboxplot_data,
    varying = names(x),
    timevar = "variable",
    idvar = "id",
    v.names = "conc",
    direction = "long"
  )

  long_data$variable2 <- vapply(
    long_data$variable,
    function(z) {
      names(x)[z]
    },
    character(1)
  )

  long_data$variable2 <- factor(
    long_data$variable2,
    levels = names(x),
    ordered = TRUE
  )

  p1 <- ggplot(long_data, aes(x = variable2, y = conc, group = variable2)) +
    geom_boxplot()

  out1 <- round(unlist(ggplot_build(p1)[["data"]][[1]][["outliers"]]), 8)
  long_data$conc <- round(long_data$conc, 8)

  long_data$is_outlier[long_data$conc %in% out1] <- "yes"

  long_data$id_outlier <- long_data$id_var
  long_data$id_outlier[!long_data$conc %in% out1] <- NA

  p2 <- ggplot(long_data, aes(x = variable2, y = conc, group = variable2)) +
    geom_boxplot() +
    geom_text(
      aes(
        label = id_outlier,
        colour = colour_by
      ),
      size = 2,
      na.rm = TRUE,
      hjust = -0.5
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
    ) +
    ylab("Standardized concentration") +
    xlab("Variable ID")
  p2
}
