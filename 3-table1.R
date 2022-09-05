# Table 1: Descriptive stats table ----

# Check dependencies
check_packages(
  cran_packages = c(""),
  bioc_packages = c("")
)

# Get vars
table1_vars <- c(
  "Case_control",
  "years_TO_scd",     
  "years_AT_scd",
  "Sex",
  "Chol_tot",
  "Age",
  "ApoB100",
  "ApoA1",
  "BMI",
  "Sbt_VIP",
  "Dbt_VIP",
  "Smoker_2fct",
  "Glc_0h",            
  "Glc_2h",
  "CRP",
  "Lpa",
  "Diabetes_Q",
  "lab_diabetes_tpq",
  "Fast_sample",
  "Education_2fct",        
  "SCD_type",
  "SCD_timetodeath",
  "BP_drug",
  "Heart_drug",
  "BZoHist_drug",
  "PPI_drug",
  "Lipid_drug",
  "Diabet_diet",
  "Diabet_pill",
  "Diabet_insulin"
)

# Get data
table1_data <- merge(
  x = preproc$raw_org$data_out["Patient_code"],
  y = preproc$raw_fia$data_out[c("Patient_code", table1_vars)],
  by = "Patient_code"
)
table1_meta <- preproc$raw_fia$varmeta_data
  
# Source func
source("./src/cc_dstats_main.R")

# Make table
table1 <- cc_dstats_main(
  df1 = table1_data[table1_vars],
  ccVar = "Case_control",
  denom_case = "Case",
  denom_ctrl = "Ctrl",
  n_digits = 3,
  use_var_metadata = TRUE,
  old_names = table1_meta$name2,
  new_names = table1_meta$name3,
  new_unit = table1_meta$unit,
  p_vals = FALSE
)