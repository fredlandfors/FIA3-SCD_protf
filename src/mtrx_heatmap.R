mtrx_heatmap <- function(text_df,
                         p_df,
                         name1,
                         title1,
                         signif_thresh = "null",
                         column_title_gp1 = gpar(fontsize = 10),
                         column_names_gp1 = gpar(fontsize = 8),
                         row_names_gp1 = gpar(fontsize = 8),
                         cell_fun_gp1 = gpar(fontsize = 5),
                         max_p = max(na.omit(p_df)),
                         color = "red",
                         na_color = "grey",
                         column_split = FALSE,
                         column_split_cat1 = "Std. Beta",
                         column_split_cat2 = "OR",
                         column_split_n = 2
                         ) {
  # Description:
  #   Plots a heatmap of  data from cov_matrix functions
  # Args:
  #   text_df: matrix with values to plot as text in each cell.
  #   p_df: matrix to colour by.
  #   name1: Name of plot
  #   title1: Title of plot
  #   signif_thresh: Significance threshold
  #   column_title_gp1: Graphical parameters.
  #   column_names_gp1: Graphical parameters.
  #   row_names_gp1: Graphical parameters.
  #   cell_fun_gp1: Graphical parameters.
  #   max_p: Largest P-value on scale.
  #   color: Heat color.
  #   na_color: NA color.
  #   heatmap_split_cat: Split categories.
  #   heatmap_split_n: After column n.
  # Dependencies:
  check_packages(
    bioc_packages = c("ComplexHeatmap", "circlize"),
    cran_packages = c("")
  )
  # Returns:
  #   table1: Table with descriptive stats by case-control status.
  library(ComplexHeatmap)
  library(circlize)

  # Calculate Bonferroni threshold
  if (signif_thresh == "null") {
    signif_thresh <- abs(log10(0.05 / (ncol(p_df) * nrow(p_df))))
  }

  # Set legend colours after signif threshold
  col_fun <- circlize::colorRamp2(
    c(0, signif_thresh, max_p), c("white", "white", color)
  )

  # Draw heatmap
  p_df <- as.matrix(p_df)
  
  if (column_split == TRUE) {
    ha_labels <- c(column_split_cat2, column_split_cat1)
    
    ha <- HeatmapAnnotation(
      # foo1 = anno_block(
      #   gp = gpar(
      #     fill = "black", col = "black"
      #   ),
      #   height = unit(0.01, "mm")
      # ),
      foo2 = anno_block(
        gp = gpar(
          fill = "white", col = "white"
        ),
        labels = ha_labels,
        labels_gp = gpar(
          fontsize = 7, fontfamily = "Helvetica", fontface = "bold"
        ),
        height = unit(3, "mm")
      )
    )
    ht <- Heatmap(
      p_df,
      col = col_fun,
      na_col = na_color,
      name = name1,
      
      column_title = title1,
      column_title_gp = column_title_gp1,
      
      cluster_columns = FALSE,
      column_names_side = "top",
      column_names_gp = column_names_gp1,
      column_split = factor(
        c(
          rep(column_split_cat1, column_split_n),
          rep(column_split_cat2, ncol(p_df) - column_split_n)
        )
      ),
      bottom_annotation = ha,
      
      cluster_rows = FALSE,
      row_names_side = "left",
      row_names_gp = row_names_gp1,
      
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(text_df[i, j], x, y, gp = cell_fun_gp1)
      },
    )
  }
  if (column_split == FALSE) {
    ht <- Heatmap(
      p_df,
      col = col_fun,
      na_col = na_color,
      name = name1,
      
      column_title = title1,
      column_title_gp = column_title_gp1,
      
      cluster_columns = FALSE,
      column_names_side = "top",
      column_names_gp = column_names_gp1,
      
      cluster_rows = FALSE,
      row_names_side = "left",
      row_names_gp = row_names_gp1,
      
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          ifelse(is.na(text_df[i, j]), " ", text_df[i, j]),
          x, y, gp = cell_fun_gp1
        )
      },
    )
  }

  return(ht)
}
