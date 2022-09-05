geom_volc <- function(df1,
                      df2) {
  # Description:
  #   Plots a volcano plot using model statistics
  # Args:
  #   df1-2: Data frame with model statistics.
  # Dependencies:
  check_packages(
    cran_packages = c("ggplot2")
  )
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(wesanderson)
  library(ggsignif)
  # Returns:
  #   p1: The volcano plot
  
  # Get data ----
  volcano_table_1 <- df1
  volcano_table_1$neg_log10_p <- abs(log10(volcano_table_1$p))
  volcano_table_1$color1 <- "Crude"
  volcano_table_1$label1 <- volcano_table_1$name
  
  volcano_table_2 <- df2
  volcano_table_2$neg_log10_p <- abs(log10(volcano_table_2$p))
  volcano_table_2$color1 <- "Full"
  volcano_table_2$label1 <- ""
  
  # Merge ----
  volcano_table_3 <- rbind(volcano_table_1, volcano_table_2)
  volcano_table_3$name3 <- rep(volcano_table_1$name, 2)
  volcano_table_3$name3[volcano_table_3$name3 == ""] <- seq(1, length(volcano_table_3$name3[volcano_table_3$name3 == ""]))
  
  # Calc p_thresh (Bonferroni) ----
  p_thresh_prot <- abs(
    log10(
      0.05/nrow(volcano_table_3[volcano_table_3$color1 == "Crude", ]) 
    )
  )
  
  # Draw labels if Crude estimate meets FWER > 5 % criteria.  ----
  volcano_table_3$label2 <- base::ifelse(
    volcano_table_3$neg_log10_p < p_thresh_prot,
    "",
    volcano_table_3$label1
  )
  
  # Draw non-significant values in plot as black points ----
  volcano_table_3$Model <- rep(NA, nrow(volcano_table_3))
  half_df <-  length(volcano_table_3$color1[volcano_table_3$color1 == "Crude"])
  
  for (i in seq_along(volcano_table_3$Model)) {
    if (volcano_table_3[i, "neg_log10_p"] > p_thresh_prot) {
      if (volcano_table_3[i, "color1"] == "Crude") {
        volcano_table_3[i, "Model"] <- "Crude"
      }
      if (volcano_table_3[i, "color1"] == "Full") {
        volcano_table_3[i, "Model"] <- "Full"
      }
    }
    
    if (volcano_table_3[i, "neg_log10_p"] < p_thresh_prot) {
      if (volcano_table_3[i, "color1"] == "Crude") {
        volcano_table_3[i, "Model"] <- "black1"
      }
      if (volcano_table_3[i, "color1"] == "Full") {
        if (volcano_table_3[i - half_df, "Model"] == "Crude") {
          volcano_table_3[i, "Model"] <- "Full"
        } else {
          volcano_table_3[i, "Model"] <- "black1"
        }
      }
    }
  }
  
  ## Volcano plot graphical parameters ----
  cols_volcano1 <- c("Crude" = "#FF0000", "Full" = "#5BBCD6", "black1" = "black")
  
  # Draw plot ----
  volcano1 <- ggplot(volcano_table_3, aes(x = OR, y = neg_log10_p, label = label2, group = name3)) +
    geom_hline(
      yintercept = p_thresh_prot, 
      linetype = "dashed", color = "black", alpha = 0.5, size = 0.5
    ) +
    annotate("text", label = "FWER = 5 %",
             x = -0.1, y = p_thresh_prot + 0.2, 
             size = 1.8, colour = "black", alpha = 0.5, fontface = "bold",
             family = "Helvetica") +
  geom_line(linetype = "dashed", color = "black", alpha = 0.5, size = 0.4) +
    geom_point(aes(color = Model), alpha = 1, size = 1.2) +
    geom_text_repel(data = subset(volcano_table_3, OR > 1.58),
                    size          = 2.5,
                    box.padding   = 0.1,
                    point.padding = 0.1,
                    force         = 0.2,
                    nudge_y       = 0.1,
                    nudge_x       = 0.5,
                    fontface      = "bold",
                    segment.size  = 0.5,
                    segment.color = "black",
                    direction     = "both"
    ) +
    geom_text_repel(data = subset(volcano_table_3, OR < 1.58),
                    size          = 2.5,
                    box.padding   = 0.1,
                    point.padding = 0.1,
                    force         = 1,
                    nudge_y       = 0.3,
                    nudge_x       = -0.5,
                    fontface      = "bold",
                    segment.size  = 0.5,
                    segment.color = "black",
                    direction     = "both",
                    max.iter      = 1000000
    ) +
    ylab("-log10(P-value)") +
    scale_y_continuous(
      breaks = seq(0, 9, 1),
      limits = c(0, 9)
    ) +
    xlab("OR per 1-SD increase in plasma protein levels") + 
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5),
      limits = c(-0.5, 2.5)
    ) +
    theme_classic() +
    theme(
      legend.justification = c(-0.05, 1), 
      legend.position = c(0, 1),
      legend.title = element_text(
        color = "black",
        size = 8,
        family = "Helvetica",
        face = "bold"
      ),
      legend.text = element_text(
        color = "black",
        size = 8,
        family = "Helvetica",
        face = "bold"
      ),
      axis.text = element_text(
        color = "black",
        size = 8,
        family = "Helvetica",
        face = "plain"
      ),
      axis.title = element_text(
        color = "black",
        size = 10,
        family = "Helvetica",
        face = "bold"
      ),
      plot.title = element_text(
        color = "black",
        size = 11,
        family = "Helvetica",
        face = "bold",
        hjust = 0.5
      )
    ) +
    scale_color_manual(values = cols_volcano1,
                       breaks = c("Crude", "Full")) + 
    ggtitle("Plasma: \n Protein-SCD effect estimate vs. P-value")
}
  