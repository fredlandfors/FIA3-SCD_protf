#-------------------------------------------------------------------------------
# Description:
#   Runs the preprocessing scripts for the FIA3-SCD-prot project
#-------------------------------------------------------------------------------

# Check dependencies
source("./src/check_packages.R")

check_packages(
  cran_packages = c("mice", "naniar", "ggplot2"),
  bioc_packages = c("")
)

# Session info
session_info <- list()
session_info$preproc <- sessionInfo()

# Get data paths
data_paths <- read.csv2(file = "./data_paths.csv", stringsAsFactors = FALSE)

# Set seed
seed = 123

#-------------------------------------------------------------------------------
# 1. Import FIA data
#-------------------------------------------------------------------------------

# Run import script
source("./src/import_fia.R")
preproc <- list()
preproc$raw_fia <- import_clinical(
  path_sample = subset(data_paths, data == "clin_sampledata", select = absolute_path)[, 1],
  path_var_meta = subset(data_paths, data == "clin_metadata", select = absolute_path)[, 1]
)

# Impute missing values using multiple chained equations
source("./src/mice_impute.R")

# Select vars to include in imputation
imp_fia_vars <- c(
  "Age", "Sex", "Glc_0h", "Glc_2h", "ApoB100", "ApoA1", "BMI", "Sbt_VIP",
  "Dbt_VIP", "Chol_tot", "Smoker_2fct", "Diabetes_Q", "Education_2fct",
  "Lpa", "CRP"
)

# Run imputation function
preproc$mice_fia <- impute_by_MICE(
  preproc$raw_fia$data_out,
  vars_include = imp_fia_vars,
  method1 = "pmm",
  set_seed = seed,
  scaling = FALSE
)

# Fin out from fia preproc
preproc$imp_fia <- preproc$mice_fia$data$mice_data_out

#-------------------------------------------------------------------------------
# 2. Import olink data
#-------------------------------------------------------------------------------

# Run import script
source("./src/import_olink.R")

preproc$raw_inf <- import_olink(
  path_sample = subset(data_paths, data == "pea_sampledata", select = absolute_path)[, 1],
  path_var_meta = subset(data_paths, data == "pea_metadata", select = absolute_path)[, 1],
  panel1 = "Olink INFLAMMATION(v.3021)",
  plate_id_var = "Plate_ID_1",
  QC_warning_var = "QC_Warning_1",
  rm_QC_warning = FALSE,
  LOD_thresh = 0.05,
  clin_data = preproc$raw_fia$data_out,
  pat_code_var = "Patient_code",
  stratum_var = "Set_FIA3"
)

preproc$raw_org <- import_olink(
  path_sample = subset(data_paths, data == "pea_sampledata", select = absolute_path)[, 1],
  path_var_meta = subset(data_paths, data == "pea_metadata", select = absolute_path)[, 1],
  panel1 = "Olink ORGAN DAMAGE(v.3311)",
  QC_warning_var = "QC_Warning_2",
  plate_id_var = "Plate_ID_2",
  rm_QC_warning = FALSE,
  LOD_thresh = 0.05,
  clin_data = preproc$raw_fia$data_out,
  pat_code_var = "Patient_code",
  stratum_var = "Set_FIA3"
)

#-------------------------------------------------------------------------------
# 3. Save to rData for use in rmarkdown
#-------------------------------------------------------------------------------
save.image(subset(data_paths, data == "markdown", select = absolute_path)[, 1])
