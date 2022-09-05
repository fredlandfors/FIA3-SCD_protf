# Figure 1A: Volcano plot ----
# Check dependencies
check_packages(
  cran_packages = c("ggrepel", "ggsci", "wesanderson", "ggplot2",
                    "ggsignif", "circlize", "grid", "cowplot", "gridGraphics"),
  bioc_packages = c("ComplexHeatmap")
)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(wesanderson)
library(ggsignif)
library(ComplexHeatmap)

# Get data
volcano_table_1 <- tableS2_1
volcano_table_1$neg_log10_p <- abs(log10(volcano_table_1$p))
volcano_table_1$color1 <- "Crude"
volcano_table_1$label1 <- volcano_table_1$name

volcano_table_2 <- tableS2_2
volcano_table_2$neg_log10_p <- abs(log10(volcano_table_2$p))
volcano_table_2$color1 <- "Full"
volcano_table_2$label1 <- ""

# Merge
volcano_table_3 <- rbind(volcano_table_1, volcano_table_2)
volcano_table_3$name3 <- rep(volcano_table_1$name, 2)
volcano_table_3$name3[volcano_table_3$name3 == ""] <- seq(1, length(volcano_table_3$name3[volcano_table_3$name3 == ""]))

# Calc p_thresh (Bonferroni)
p_thresh_prot <- abs(
  log10(
    0.05/nrow(volcano_table_3[volcano_table_3$color1 == "Crude", ]) 
  )
)
p_thresh_prot_1 <- abs(
  log10(
    0.01/nrow(volcano_table_3[volcano_table_3$color1 == "Crude", ]) 
  )
)
p_thresh_prot_01 <- abs(
  log10(
    0.001/nrow(volcano_table_3[volcano_table_3$color1 == "Crude", ]) 
  )
)

# Draw labels if Crude estimate meets FWER > 5 % criteria. 
volcano_table_3$label2 <- base::ifelse(
  volcano_table_3$neg_log10_p < p_thresh_prot,
  "",
  volcano_table_3$label1
)

# Draw non-significant values in plot as black points
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

# Draw plot
volcano1 <- ggplot(volcano_table_3, aes(x = OR, y = neg_log10_p, label = label2, group = name3)) +
  geom_hline(
    yintercept = p_thresh_prot, 
    linetype = "dashed", color = "black", alpha = 0.5, size = 0.5
  ) +
  annotate("text", label = "FWER = 5 %",
           x = -0.1, y = p_thresh_prot + 0.2, 
           size = 1.8, colour = "black", alpha = 0.5, fontface = "bold",
           family = "Helvetica") +
  # geom_hline(
  #   yintercept = p_thresh_prot_1, 
  #   linetype = "dashed", color = "black", alpha = 0.5, size = 0.5
  # ) +
  # annotate("text", label = "FWER = 1 %",
  #          x = -0.1, y = p_thresh_prot_1 + 0.2, 
  #          size = 1.8, colour = "black", alpha = 0.5, fontface = "bold",
  #          family = "Helvetica") +
  # geom_hline(
  #   yintercept = p_thresh_prot_01, 
  #   linetype = "dashed", color = "black", alpha = 0.5, size = 0.5
  # ) +
  # annotate("text", label = "FWER = 0.1 %",
  #          x = -0.1, y = p_thresh_prot_01 + 0.2, 
  #          size = 1.8, colour = "black", alpha = 0.5, fontface = "bold",
  #          family = "Helvetica") +
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
                     breaks = c("Crude", "Full")) 
#+ ggtitle("Plasma: \n Protein-SCD effect estimate vs. P-value")


# Figure 1B: Boxviolin plot ----
# Get variables
boxviolin_vars <- c(
  c("Case_control"),
  names(preproc$raw_org$data_out)[names(preproc$raw_org$data_out) %in% preproc$raw_org$varmeta_data$name1],
  names(preproc$raw_inf$data_out)[names(preproc$raw_inf$data_out) %in% preproc$raw_inf$varmeta_data$name1]
)

# Get data
boxviolin_data_1 <- merge(
  x = preproc$raw_fia$data_out,
  y = preproc$raw_inf$data_out,
  by = "Patient_code"
)

boxviolin_data_2 <- merge(
  x = boxviolin_data_1,
  y = preproc$raw_org$data_out,
  by = "Patient_code"
)

boxviolin_data_2 <- boxviolin_data_2[boxviolin_vars]
boxviolin_data_2 <- rapply(
  boxviolin_data_2,
  scale,
  classes = "numeric", how = "replace"
)

# Set factor levels
levels(boxviolin_data_2$Case_control) <- c("Control", "SCD")


## Boxviolin graphical parameters ----
boxviolin_colors <- c("Control" = "#F2AD00", "SCD" = "#F2AD00")
boxviolin_axis.text.y = element_text(
  color = "black",
  size = 8,
  family = "Helvetica",
  face = "plain"
)
boxviolin_axis.text.x = element_text(
  color = "black",
  size = 10,
  family = "Helvetica",
  face = "bold",
  angle = 45,
  vjust = 1,
  hjust = 1
)
boxviolin_axis.title = element_text(
  color = "black",
  size = 10,
  family = "Helvetica",
  face = "bold"
)
boxviolin_strip.text = element_text(
  color = "black",
  size = 10,
  family = "Helvetica",
  face = "bold"
)
boxviolin_ggtitle = element_text(
  color = "black",
  size = 11,
  family = "Helvetica",
  face = "bold",
  hjust = 0.5
)

# Draw plots
boxviolin1 <- ggplot(
  data = boxviolin_data_2,
  aes(x = Case_control, y = LTA4H)
) +
  geom_violin(alpha = 0) +
  geom_boxplot(alpha = 0) +
  geom_jitter(aes(color = Case_control), size = 0.5, alpha = 0.6, width = 0.1) +
  ylab("LTA4H (SD)") +
  scale_y_continuous(
    limits = c(
      min(boxviolin_data_2$HGF) - sd(boxviolin_data_2$HGF),
      max(boxviolin_data_2$HGF) + sd(boxviolin_data_2$HGF)
    ),
    breaks = seq(-3, 6, 1)
  ) +
  xlab("") +
  theme_classic() +
  theme(
    axis.text.y = boxviolin_axis.text.y,
    axis.text.x = boxviolin_axis.text.x,
    axis.title = boxviolin_axis.title,
    strip.text = boxviolin_strip.text,
    legend.position = "none",
    plot.title = boxviolin_ggtitle
  ) +
  geom_signif(
    comparisons = list(c("Control", "SCD")),
    map_signif_level = TRUE
  ) +
  scale_color_manual(
    values = boxviolin_colors,
    breaks = c("Control", "SCD")
  ) + 
  ggtitle("LTA4H")

boxviolin2 <- ggplot(
  data = boxviolin_data_2,
  aes(x = Case_control, y = HGF)
) +
  geom_violin(alpha = 0) +
  geom_boxplot(alpha = 0) +
  geom_jitter(aes(color = Case_control), size = 0.5, alpha = 0.6, width = 0.1) +
  ylab("HGF (SD)") +
  scale_y_continuous(
    limits = c(
      min(boxviolin_data_2$HGF) - sd(boxviolin_data_2$HGF),
      max(boxviolin_data_2$HGF) + sd(boxviolin_data_2$HGF)
    ),
    breaks = seq(-3, 6, 1)
  ) +
  xlab("") +
  theme_classic() +
  theme(
    axis.text.y = boxviolin_axis.text.y,
    axis.text.x = boxviolin_axis.text.x,
    axis.title = boxviolin_axis.title,
    strip.text = boxviolin_strip.text,
    legend.position = "none",
    plot.title = boxviolin_ggtitle
  ) +
  geom_signif(
    comparisons = list(c("Control", "SCD")),
    map_signif_level = TRUE
  ) +
  scale_color_manual(
    values = boxviolin_colors,
    breaks = c("Control", "SCD")
  ) + 
  ggtitle("HGF")


# Figure 1C: Correlation heatmap with covariates ----
## Heatmap graphical parameters ----
heatmap_column_title_gp = gpar(fontsize = 10, fontfamily = "Helvetica", fontface = "bold")
heatmap_column_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica", fontface = "bold")
heatmap_row_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica", fontface = "bold")
heatmap_cell_fun_gp = gpar(fontsize = 5, fontfamily = "Helvetica", fontface = "bold")
heatmap_color = "#FF0000"

## Data wrangling ----

# Get vars
heatmap_clinvars <- c(
  "Case_control", "Sex", "Age", "BMI", "Chol_tot", "Sbt_VIP",
  "Dbt_VIP", "Glc_0h", "Smoker_2fct", "Education_2fct"
)
heatmap_mtrx_colvars <- c(
  "BMI", "Chol_tot", "ApoB100", "ApoA1", "Lpa", "CRP",
  "Sbt_VIP", "Dbt_VIP", "Glc_0h", "Glc_2h",
  "Smoker_2fct", "Education_2fct"
)
heatmap_mtrx_colnames <- c(
  "BMI", "Cholesterol", "ApoB-100", "ApoA-1", "Lp(a)",
  "CRP", "SBT", "DBT", "Glucose", "OGTT",
  "Smoking", "Education"
)
heatmap_mtrx_adjvars <- c("Age", "Sex")

heatmap_biom <- subset(volcano_table_1, neg_log10_p > p_thresh_prot)[, 1] 

heatmap_vars <- unique(
  c("Patient_code", heatmap_clinvars, heatmap_mtrx_colvars, heatmap_biom)
)

# Get data
heatmap_data_1 <- merge(
  x = preproc$raw_fia$data_out,
  y = preproc$raw_inf$data_out,
  by = "Patient_code"
)
heatmap_data_2 <- merge(
  x = heatmap_data_1,
  y = preproc$raw_org$data_out,
  by = "Patient_code"
)
heatmap_data_2 <- heatmap_data_2[heatmap_vars]

# Log 2 transform non-normal distributed clinical vars
heatmap_data_2$Glc_0h <- log(heatmap_data_2$Glc_0h, base = 2)
heatmap_data_2$Glc_2h <- log(heatmap_data_2$Glc_2h, base = 2)
heatmap_data_2$CRP <- log(heatmap_data_2$CRP + 0.03, base = 2)
heatmap_data_2$Lpa <- log(heatmap_data_2$Lpa, base = 2)

## I. Controls ----
heatmap_data_ctrl <- subset(heatmap_data_2, Case_control == "Ctrl")

# Source func
source("./src/mtrx_glm.R")

# Calculate matrix
heatmap_mtrx_ctrl <- mtrx_glm(
  heatmap_data_ctrl,
  col_vars = heatmap_mtrx_colvars,
  row_vars_all = heatmap_biom,
  adj_vars = heatmap_mtrx_adjvars
)

# Log 10 transform P-values
heatmap_mtrx_ctrl$p_vals <- sapply(
  heatmap_mtrx_ctrl$p_vals,
  function(x) {
    y <- -log10(x)
    return(y)
  }
)
row.names(heatmap_mtrx_ctrl$p_vals) <- row.names(heatmap_mtrx_ctrl$coef)

# Change col names for ctrls
colnames(heatmap_mtrx_ctrl$coef) <- heatmap_mtrx_colnames
colnames(heatmap_mtrx_ctrl$p_vals) <- heatmap_mtrx_colnames

# Source heatmap script
source("./src/mtrx_heatmap.R")

# Draw heatmap
heatmap_p_thresh <- -log10(
  0.05 / (nrow(heatmap_mtrx_ctrl$coef) * ncol(heatmap_mtrx_ctrl$coef))
)
ht1 <- mtrx_heatmap(
  text_df = heatmap_mtrx_ctrl$coef,
  p_df = heatmap_mtrx_ctrl$p_vals,
  name1 = "Neglog10p_ctrl",
  title1 = "I. Controls",
  column_title_gp1 = heatmap_column_title_gp,
  column_names_gp1 = heatmap_column_names_gp,
  row_names_gp1 = heatmap_row_names_gp,
  cell_fun_gp1 = heatmap_cell_fun_gp,
  signif_thresh = heatmap_p_thresh,
  max_p = 15,
  color = heatmap_color,
  column_split = TRUE,
  column_split_cat1 = "(Std. beta)",
  column_split_cat2 = "(OR)",
  column_split_n = 10
)

## II. Cases ----
heatmap_data_case <- subset(heatmap_data_2, Case_control == "Case")

# Calculate matrix
heatmap_mtrx_case <- mtrx_glm(
  heatmap_data_case,
  col_vars = heatmap_mtrx_colvars,
  row_vars_all = heatmap_biom,
  adj_vars = heatmap_mtrx_adjvars
)

# Log 10 transform P-values
heatmap_mtrx_case$p_vals <- sapply(
  heatmap_mtrx_case$p_vals,
  function(x) {
    y <- -log10(x)
    return(y)
  }
)
row.names(heatmap_mtrx_case$p_vals) <- row.names(heatmap_mtrx_case$coef)

# Change col names for cases
colnames(heatmap_mtrx_case$coef) <- heatmap_mtrx_colnames
colnames(heatmap_mtrx_case$p_vals) <- heatmap_mtrx_colnames

# Draw heatmap
heatmap_p_thresh <- -log10(
  0.05 / (nrow(heatmap_mtrx_case$coef) * ncol(heatmap_mtrx_case$coef))
)
ht2 <- mtrx_heatmap(
  text_df = heatmap_mtrx_case$coef,
  p_df = heatmap_mtrx_case$p_vals,
  name1 = "Neglog10p_case",
  title1 = "II. Cases",
  column_title_gp1 = heatmap_column_title_gp,
  column_names_gp1 = heatmap_column_names_gp,
  row_names_gp1 = heatmap_row_names_gp,
  cell_fun_gp1 = heatmap_cell_fun_gp,
  signif_thresh = heatmap_p_thresh,
  max_p = 15,
  color = heatmap_color,
  column_split = TRUE,
  column_split_cat1 = "(Std. Beta)",
  column_split_cat2 = "(OR)",
  column_split_n = 10
)

## Combine heatmaps ----
lgd1 <- Legend(
  col_fun = circlize::colorRamp2(
    c(0, heatmap_p_thresh, 15), c("white", "white", heatmap_color)
  ),
  title = "-log10(P-value)",
  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica"),
  labels_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica"),
  direction = "vertical"
)

grob <- grid.grabExpr(
  draw(
    ht1 + ht2,
    main_heatmap = "Neglog10p_case",
    ht_gap = unit(1, "cm"),
    show_heatmap_legend = FALSE,
    annotation_legend_list = lgd1,
    annotation_legend_side = "right",
    #column_title = "Plasma: Risk factor vs. protein",
    column_title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Helvetica")
  )
)

# Figure 1: Combined ----
figure1 <- cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(
      NULL,
      ggplot() +
        ggtitle("Plasma: \n Protein-SCD effect estimate vs. P-value") + 
        theme(
          plot.title = element_text(
            color = "black",
            size = 11,
            family = "Helvetica",
            face = "bold",
            hjust = 0.5
          )
        ),
      volcano1,
      NULL,
      nrow = 4, ncol = 1,
      rel_heights = c(4, 8, 85, 3)
    ),
    cowplot::plot_grid(
      NULL,
      ggplot() +
        ggtitle("Plasma: SCD vs. protein levels") +
        theme(plot.title = boxviolin_ggtitle), 
      NULL,
      cowplot::plot_grid(boxviolin2, boxviolin1, ncol = 2),
      NULL, 
      nrow = 5, ncol = 1,
      rel_heights = c(4, 4, 4, 82, 6)
    ),
    rel_widths = c(50, 50),
    nrow = 1,
    labels = c("A", "B")
  ),
  cowplot::plot_grid(
    NULL,
    ggplot() +
      ggtitle("Plasma: Risk factor vs. protein levels") + 
      theme(
        plot.title = element_text(
          color = "black",
          size = 11,
          family = "Helvetica",
          face = "bold",
          hjust = 0.5
        )
    ),
    grob,
    nrow = 3, ncol = 1,
    rel_heights = c(5, 2, 93)
  ),
  nrow = 2,
  rel_heights = c(1, 1),
  labels = c("", "C")
)

figure2 <- cowplot::plot_grid(
  NULL,
  ggplot() +
    ggtitle("Plasma: Risk factor vs. protein levels") + 
    theme(
      plot.title = element_text(
        color = "black",
        size = 11,
        family = "Helvetica",
        face = "bold",
        hjust = 0.5
      )
    ),
  grob,
  nrow = 3, ncol = 1,
  rel_heights = c(5, 2, 93)
)


