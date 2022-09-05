geom_rwise_boxplot <- function(x,
                              id = as.character(seq(1, nrow(x))),
                              fill1 = "none") {
  # Description:
  #   Plots a row-wise boxplot
  # Args:
  #   x: PCA x numerical data.
  #   id: Sample ID var (obs! as df with one column).
  #   fill1: column to fill plot by.
  # Dependencies:
  check_packages(
    cran_packages = c("ggplot2")
  )
  # Returns:
  #   p2: The boxplot
  
  # Get df and scale numerical vars
  gboxplot_data <- rapply(
    x,
    scale,
    classes = "numeric",
    how = "replace"
  )

  gboxplot_data <- data.frame(
    t(gboxplot_data)
  )

  colnames(gboxplot_data) <- id

  gboxplot_data <- cbind(rownames(gboxplot_data), gboxplot_data)
  names(gboxplot_data)[1] <- "variable_name"

  long_data <- reshape(
    gboxplot_data,
    varying = names(gboxplot_data)[-1],
    timevar = "sample_no",
    idvar = "variable_no",
    v.names = "conc",
    direction = "long"
  )

  long_data$sample_no2 <- rep(NA, nrow(long_data))

  long_data$sample_no2 <- sapply(
    long_data$sample_no,
    function(x) {
      id[x]
    }
  )

  p1 <- ggplot(long_data, aes(x = sample_no2, y = conc, group = sample_no2)) +
    geom_boxplot()

  out1 <- round(unlist(ggplot_build(p1)[["data"]][[1]][["outliers"]]), 8)
  long_data$conc <- round(long_data$conc, 8)

  long_data$is_outlier[long_data$conc %in% out1] <- "yes"

  long_data$id_outlier <- long_data$variable_name
  long_data$id_outlier[!long_data$conc %in% out1] <- NA

  if ("none" %in% fill1) {
    p2 <- ggplot(long_data, aes(x = sample_no2, y = conc, group = sample_no2)) +
      geom_boxplot() +
      geom_text(
        aes(
          label = id_outlier,
          colour = id_outlier
        ),
        size = 2,
        na.rm = TRUE,
        hjust = -0.3
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      ) +
      ylab("Standardized concentration") +
      xlab("Sample ID")
    p2
  } else {
    long_data$col_fill <- vapply(
      long_data$sample_no,
      function(x) {
        fill1[x]
      },
      character(1)
    )

    p2 <- ggplot(long_data, aes(x = sample_no2, y = conc, group = sample_no2)) +
      geom_boxplot(aes(fill = col_fill)) +
      geom_text(
        aes(
          label = id_outlier,
          colour = id_outlier
        ),
        size = 2,
        na.rm = TRUE,
        hjust = -0.3
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      ) +
      ylab("Standardized concentration") +
      xlab("Sample ID")
    p2
  }
}
