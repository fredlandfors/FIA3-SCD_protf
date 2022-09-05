plot_pca <- function(x,
                     c,
                     id = data.frame(id = seq(1, nrow(x)), stringsAsFactors = FALSE),
                     opt_ellipse = "auto",
                     n_pc = 9,
                     return_fit = TRUE) {
  # Description:
  #   Performs PCA analysis of MS data colour for pre-analytical factor.
  # Args:
  #   x: PCA x numerical data.
  #   c: Variable to color by (obs! as df with one column).
  #   id: Sample ID var (obs! as df with one column).
  #   opt_ellipse: "auto", "groupwise" or "globalonly". Bivariate t distribution
  #                ellipse.
  #   n_pc: Number of principal components to plot
  #   return_fit: Return fit data frame (true/false).
  # Dependencies:
  check_packages(
    cran_packages = c("ggplot2")
  )
  # Returns:
  #   pca_fit: PCA model.
  #   ..._plot: Variance explained, scores and loadings plot.
  
  # Fit PCA
  pca_data <- x
  pca_fit <- prcomp(pca_data, center = TRUE, scale. = TRUE)

  # p1: Variance explained barplot
  p1_data <- t(data.frame(importance = summary(pca_fit)$importance[2, ]))
  p1_data <- data.frame(p1_data)

  p1_longdata <- reshape(
    p1_data,
    varying = grep("^PC", names(p1_data), value = TRUE),
    timevar = "PC",
    idvar = "id",
    v.names = "explained_variance",
    direction = "long"
  )

  p1 <- ggplot(p1_longdata, aes(x = PC, y = explained_variance)) +
    geom_col()

  # p2: PC1 vs PCn scores plot
  if (inherits(c[, 1], "Date")) {
    c[, 1] <- sapply(
      c,
      function(x) {
        round(((as.numeric(x) - 6159) / 365.25 + 1986.86242), 0)
      }
    )
  }

  if (is.character(c[, 1]) && length(unique(c[, 1])) > 10) {
    c[, 1] <- as.numeric(c[, 1])
  }

  if (is.factor(c[, 1]) && length(unique(c[, 1])) > 10) {
    c[, 1] <- as.numeric(c[, 1])
  }

  p2_data <- data.frame(
    id,
    c,
    pca_fit$x,
    stringsAsFactors = FALSE #OBS
  )

  p2_longdata <- reshape(
    p2_data,
    varying = grep("^PC", names(p2_data), value = TRUE)[-1],
    timevar = "PC",
    idvar = "id",
    v.names = "pc.score",
    direction = "long"
  )
  p2_longdata$PC <- p2_longdata$PC + 1
  p2_longdata <- subset(p2_longdata, PC <= n_pc)

  if (opt_ellipse == "auto" && length(unique(c[, 1])) <= 3) {
    opt_ellipse <- "groupwise"
  }
  else {
    opt_ellipse <- "globalonly"
  }

  if (opt_ellipse == "groupwise") {
    p2 <- ggplot(p2_longdata, aes(x = PC1, y = pc.score)) +
      geom_point(aes_string(colour = names(c))) +
      geom_text(aes_string(label = names(id), colour = names(c)),
        size = 2, hjust = -0.5
      ) +
      stat_ellipse(type = "t", size = 0.5) +
      stat_ellipse(aes_string(colour = names(c)),
        type = "t", linetype = 2, size = 0.5
      ) +
      facet_wrap(~PC, labeller = label_both)
  }

  if (opt_ellipse == "globalonly") {
    p2 <- ggplot(p2_longdata, aes(x = PC1, y = pc.score)) +
      geom_point(aes_string(colour = names(c))) +
      geom_text(aes_string(label = names(id), colour = names(c)),
        size = 2, hjust = -0.5
      ) +
      stat_ellipse(type = "t", size = 0.5) +
      facet_wrap(~PC, labeller = label_both)
  }

  if (!opt_ellipse %in% c("groupwise", "globalonly")) {
    stop("opt_ellipse must be specified.")
  }

  # p3: Variable loadings plot
  p3_data <- data.frame(
    pca_fit$rotation,
    var_name = rownames(pca_fit$rotation)
  )

  p3_longdata <- reshape(
    p3_data,
    varying = grep("^PC", names(p3_data), value = TRUE)[-1],
    timevar = "PC",
    idvar = "id",
    v.names = "var_loading",
    direction = "long"
  )

  p3_longdata$PC <- p3_longdata$PC + 1
  p3_longdata <- subset(p3_longdata, PC <= n_pc)

  p3 <- ggplot(p3_longdata, aes(x = PC1, y = var_loading)) +
    geom_point() +
    geom_text(aes(label = var_name),
      size = 2, hjust = -0.2
    ) +
    geom_hline(yintercept = 0, size = 0.3, linetype = 2) +
    geom_vline(xintercept = 0, size = 0.3, linetype = 2) +
    scale_x_continuous(limits = c(
      min(p3_longdata$PC1) - sd(p3_longdata$PC1),
      max(p3_longdata$PC1) + sd(p3_longdata$PC1)
      )
    ) +
    scale_y_continuous(limits = c(
      min(p3_longdata$var_loading) - sd(p3_longdata$var_loading),
      max(p3_longdata$var_loading) + sd(p3_longdata$var_loading)
      )
    ) + 
    facet_wrap(~ PC, labeller = label_both, strip.position = "left") 

  # p4: PC1 vs PC2
  p4_data <- data.frame(
    id,
    c,
    pca_fit$x
  )

  if (opt_ellipse == "groupwise") {
    p4 <- ggplot(p4_data, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(colour = names(c))) +
      geom_text(aes_string(label = names(id), colour = names(c)),
        size = 2, hjust = -0.5
      ) +
      stat_ellipse(type = "t", size = 0.5) +
      stat_ellipse(aes_string(colour = names(c)),
        type = "t", linetype = 2, size = 0.5
      )
  }

  if (opt_ellipse == "globalonly") {
    p4 <- ggplot(p4_data, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(colour = names(c))) +
      geom_text(aes_string(label = names(id), colour = names(c)),
        size = 2, hjust = -0.5
      ) +
      stat_ellipse(type = "t", size = 0.5)
  }

  if (!opt_ellipse %in% c("groupwise", "globalonly")) {
    stop("opt_ellipse must be specified.")
  }

  # Return
  if (return_fit == TRUE) {
    return_list <- list(
      fit = pca_fit,
      var_expl = p1,
      scores_plot = p2,
      loadings_plot = p3,
      pc1vs2_plot = p4
    )
  } else {
    return_list <- list(
      pc1vs2_plot = p4,
      scores_plot = p2,
      var_expl = p1,
      loadings_plot = p3
    )
  }
  

  return(return_list)
}
