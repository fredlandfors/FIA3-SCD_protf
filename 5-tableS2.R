# Table S2: Conditional logistic regression models 1 & 2 ----
# Check dependencies
check_packages(
  cran_packages = c("survival"),
  bioc_packages = c("")
)

# Get vars
tableS2_vars <- c(
  names(preproc$raw_org$data_out)[names(preproc$raw_org$data_out) %in% preproc$raw_org$varmeta_data$name1],
  names(preproc$raw_inf$data_out)[names(preproc$raw_inf$data_out) %in% preproc$raw_inf$varmeta_data$name1]
)

# Get data
tableS2_data_1 <- merge(
  x = preproc$imp_fia,
  y = preproc$raw_inf$data_out,
  by = "Patient_code"
)

tableS2_data_2 <- merge(
  x = tableS2_data_1,
  y = preproc$raw_org$data_out,
  by = "Patient_code"
)

tableS2_meta_1 <- rbind(
  preproc$raw_inf$varmeta_data,
  preproc$raw_org$varmeta_data
)

tableS2_meta_2 <- subset(
  tableS2_meta_1,
  name1 %in% tableS2_vars
)

# Source function
source("./src/cc_multiclogit.R")

# Model 1: fit conditional logistic regression
tableS2_1 <- cc_multiclogit(
      df1 = tableS2_data_2, 
      ccVar = "Case_control",
      adj = "",
      stratum = "Set_FIA3",
      x_vars = tableS2_vars,
      scaling = TRUE,
      ci_adj = 122
)

# Model 2: fit conditional logistic regression with adjustment
tableS2_2 <- cc_multiclogit(
      df1 = tableS2_data_2, 
      ccVar = "Case_control",
      adj = " + BMI + Chol_tot + Sbt_VIP + Glc_0h + Smoker_2fct",
      stratum = "Set_FIA3",
      x_vars = tableS2_vars,
      scaling = TRUE,
      ci_adj = 122
)

# Merge
tableS2_xlsx <- list()

tableS2_xlsx$Model1 <- merge(
  x = tableS2_meta_2,
  y = tableS2_1,
  by.x = "name1",
  by.y = "name",
  sort = FALSE
)

tableS2_xlsx$Model2 <- merge(
  x = tableS2_meta_2,
  y = tableS2_2,
  by.x = "name1",
  by.y = "name",
  sort = FALSE
)