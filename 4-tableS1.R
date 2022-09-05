# Table S1: Descriptive stats suppl. ----
# Check dependencies
check_packages(
  cran_packages = c(""),
  bioc_packages = c("")
)

# Get vars
tableS1_vars <- c(
  names(preproc$raw_org$data_out)[names(preproc$raw_org$data_out) %in% preproc$raw_org$varmeta_data$name1],
  names(preproc$raw_inf$data_out)[names(preproc$raw_inf$data_out) %in% preproc$raw_inf$varmeta_data$name1]
)

# Get data
tableS1_data_1 <- merge(
  x = preproc$raw_fia$data_out,
  y = preproc$raw_inf$data_out,
  by = "Patient_code"
)

tableS1_data_2 <- merge(
  x = tableS1_data_1,
  y = preproc$raw_org$data_out,
  by = "Patient_code"
)

tableS1_meta_1 <- rbind(
  preproc$raw_inf$varmeta_data,
  preproc$raw_org$varmeta_data
)

tableS1_meta_2 <- subset(
  tableS1_meta_1,
  name1 %in% tableS1_vars
)

# Source func
source("./src/cc_dstats_supp.R")

# Make table
tableS1 <- cc_dstats_supp(
  df1 = tableS1_data_2,
  ccVar = "Case_control",
  denom_case = "Case",
  denom_ctrl = "Ctrl",
  signif_digits = 3,
  include_var_meta = TRUE,
  var_meta_df = tableS1_meta_2,
  var_meta_merge_by = "name1"
)