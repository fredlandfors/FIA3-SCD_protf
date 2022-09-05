import_olink <- function(path_sample,
                         path_var_meta,
                         panel1,
                         plate_id_var,
                         QC_warning_var,
                         rm_QC_warning,
                         LOD_thresh = 0.05,
                         clin_data = NULL,
                         pat_code_var = NULL,
                         stratum_var = NULL) {
  # Description:
  #   Performs data cleaning of clinical  data.
  # Args:
  #   path_sample: Path for sample data.
  #   path_varmeta: Path for variable metadata.
  #   panel1: Which panel, e.g. "Organ damage".
  #   plate_id_var: Plate id variable
  #   QC_warning_var: Quality control warning variable.
  #   rm_QC_warning: Remove samples with olink QC warning (true/false).
  #   LOD_thresh: Remove variables over x threshold.
  #   clin_data: Clinical df to use for removal of duplicate samples.
  #   pat_code_var: clin_data patient code var name as string.
  #   stratum_var: clin_data stratum var name as string.
  # Dependencies:
  check_packages(
    bioc_packages = c(""),
    cran_packages = c("")
  )
  # Returns:
  #   return_list: List of objects.
  # Constants:
  # Error handling:
  # Define lists and functions:

  return_list <- list()

  #----------------------------------------------------------------------------
  # Import data
  #----------------------------------------------------------------------------

  # Import variable metadata:
  olink_meta <- read.csv2(
    path_var_meta,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = rep("character", 7)
  )
  
  # Import sample data:
  olink_sample <- read.csv2(
    path_sample,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = olink_meta[, "data_type"],
  )

  #----------------------------------------------------------------------------
  # Clean 
  #----------------------------------------------------------------------------
  # Assign appropriate names to var metadata
  olink_meta$name1 <- make.names(olink_meta$name1)
  names(olink_sample) <- olink_meta$name1
  
  # Return variable meta data
  return_list$varmeta_data <- subset(olink_meta, panel == panel1)

  # Filter samples with QC warnings if rm_QC_warning == TRUE
  if (rm_QC_warning == TRUE) {
    if (!is.null(QC_Warning_var)) {
      row_index_select1 <- olink_sample[, QC_Warning_var] == "Pass" 
      olink_sample <- olink_sample[row_index_select1, ]
    } else {
      stop("Error in import_olink: must supply QC_warning variable.")
    }
  }
  
  # Subset only the complete case_control sets if clinical data set is present
  if (!is.null(clin_data)) {
    olink_sample_complete <- merge(
      x = clin_data[c(pat_code_var, stratum_var)],
      y = olink_sample,
      by = pat_code_var
    )
    check_table <- table(olink_sample_complete[, stratum_var])
    row_index_select2 <- olink_sample_complete[, stratum_var] %in% names(check_table[check_table == 2])
    olink_sample <- olink_sample_complete[row_index_select2, ]
  }

  # Filter variables with values below LOD % threshold
  b_LOD_v <- seq(0, LOD_thresh * 100)
  b_LOD_v_char <- sapply(b_LOD_v, function(x) {
    paste0(x, "%")
  })
  olink_meta2 <- subset(
    olink_meta,
    below_LOD_freq %in% b_LOD_v_char | is.na(below_LOD_freq)
  )
  olink_cleanvars <- olink_meta2[, "name1"]
  olink_sample2 <- olink_sample[olink_cleanvars]

  # Filter by olink panel
  olink_sample3 <- olink_sample2[c(
    "Patient_code",
    plate_id_var,
    QC_warning_var,
    olink_meta2$name1[olink_meta2$panel == panel1]
  )]
  return_list$data_out <- olink_sample3

  # Return data
  return(return_list)
}
