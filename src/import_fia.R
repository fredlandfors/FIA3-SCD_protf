import_clinical <- function(path_sample, path_var_meta) {
  # Description:
  #   Performs data cleaning of data from the FIA-VHU-MONICA database.
  # Args:
  #   path_sample: Path for sample data.
  #   path_meta: Path for variable metadata.
  #   method: Which method was used e.g. "LC pos target".
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
  # IMPORT
  #----------------------------------------------------------------------------

  # Import variable metadata:
  clin_meta <- read.csv2(
    path_var_meta,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = rep("character", 7)
  )
  length_include <- seq(1, length(clin_meta$name1[clin_meta$name1 != ""])) # deals with added vars in metadata
  return_list$varmeta_data <- clin_meta

  # Import sample data:
  clin_sample <- read.csv2(
    path_sample,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = clin_meta[length_include, "data_type"],
    na.strings = c("", 5555, 6666, 7777, 8888, 9999)
  )
  return_list$clin_import <- clin_sample

  #----------------------------------------------------------------------------
  # CLEAN
  #----------------------------------------------------------------------------

  # Rename variables:
  check_rename <- data.frame(
    var1 = names(clin_sample),
    var2 = clin_meta$name1[length_include],
    var3 = clin_meta$name2[length_include]
  )
  return_list$check_rename <- check_rename
  names(clin_sample) <- clin_meta$name2[length_include]
  
  # Remove duplicates (some FIA2 controls later became cases in FIA3)
  check_table <- table(clin_sample$Patient_code)
  check_table_names <- names(check_table[check_table > 1])
  
  # Subset the cases of the duplicates.
  clin_sample <- subset(
    clin_sample, 
    (Patient_code %in% check_table_names & Case_control == "F") |
    !Patient_code %in% check_table_names
  )
  
  return_list$renamed_data <- clin_sample

  # Variable cleaning:
  # Case_control
  # Changing case/control status variable from "F"/"K" to "Case"/"Ctrl" 
  clin_sample$Case_control <- as.numeric(clin_sample$Case_control)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (clin_sample$Case_control[i] == 2) {
      clin_sample$Case_control[i] <- 0
    }
  }
  clin_sample$Case_control <- as.factor(clin_sample$Case_control)
  levels(clin_sample$Case_control) <- c("Ctrl", "Case")

  # Sex
  # Changing sex status variable from "2"/"1" to "female"/"male" 
  clin_sample$Sex <- as.numeric(clin_sample$Sex)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (clin_sample$Sex[i] == 2) {
      clin_sample$Sex[i] <- 0
    }
  }
  clin_sample$Sex <- as.factor(clin_sample$Sex)
  levels(clin_sample$Sex) <- c("female", "male")

  # SCD_type
  # Setting NA's as controls, if they are a control
  clin_sample$SCD_type <- as.numeric(clin_sample$SCD_type)
  table(clin_sample$SCD_type)
  table(as.numeric(clin_sample$SCD_type))
  table(as.factor(as.numeric(clin_sample$SCD_type)))
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (is.na(clin_sample$SCD_type[i]) == TRUE && clin_sample$Case_control[i] == "Ctrl") {
      clin_sample$SCD_type[i] <- 7
    }
  }
  table(clin_sample$SCD_type)
  clin_sample$SCD_type <- as.factor(clin_sample$SCD_type)
  table(clin_sample$SCD_type)
  levels(clin_sample$SCD_type) <- c("LBBB", "NSTEMI", "unknown", "PM", "SCD", "STEMI", "ctrl")
  table(clin_sample$SCD_type)

  # Fast_questionnaire
  # Changing  variable from "2"/"1" to "no"/"yes" for readability
  clin_sample$Fast_questionnaire <- as.numeric(clin_sample$Fast_questionnaire)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (is.na(clin_sample$Fast_questionnaire[i]) == FALSE) {
      if (clin_sample$Fast_questionnaire[i] == 2) {
        clin_sample$Fast_questionnaire[i] <- 0
      }
    }
  }
  clin_sample$Fast_questionnaire <- as.factor(clin_sample$Fast_questionnaire)
  levels(clin_sample$Fast_questionnaire) <- c("no", "yes")

  # Fast_sample
  # Changing  variable from 0/1/2/3 to "0-4h", "4-6h", "6-8h", ">8h" 
  clin_sample$Fast_sample <- as.factor(clin_sample$Fast_sample)
  levels(clin_sample$Fast_sample) <- c("0-4h", "4-6h", "6-8h", ">8h")

  # Diabetes_Q
  # Changing  variable from "2"/"1" to "no"/"yes" 
  clin_sample$Diabetes_Q <- as.numeric(clin_sample$Diabetes_Q)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (is.na(clin_sample$Diabetes_Q[i]) == FALSE) {
      if (clin_sample$Diabetes_Q[i] == 2) {
        clin_sample$Diabetes_Q[i] <- 0
      }
    }
  }
  clin_sample$Diabetes_Q <- as.factor(clin_sample$Diabetes_Q)
  levels(clin_sample$Diabetes_Q) <- c("no", "yes")
  
  # Smoker_2fct
  # Changing  variable from Rökstatus:
  # 1 = current smoker, 2 = Ex smoker, 3 = Never-smoker,  4 = Occasional smoker, 5 = Former occasional smoker
  # to
  # 0 = no (including ex smokers), 1 = yes (including occasional smokers)
  clin_sample$Smoker_2fct <- as.numeric(clin_sample$Smoker)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (is.na(clin_sample$Smoker_2fct[i]) == FALSE) {
      # run recode
      if (clin_sample$Smoker[i] == 1) {
        clin_sample$Smoker_2fct[i] <- 1
      }
      if (clin_sample$Smoker[i] == 2) {
        clin_sample$Smoker_2fct[i] <- 0
      }
      if (clin_sample$Smoker[i] == 3) {
        clin_sample$Smoker_2fct[i] <- 0
      }
      if (clin_sample$Smoker[i] == 4) {
        clin_sample$Smoker_2fct[i] <- 1
      }
      if (clin_sample$Smoker[i] == 5) {
        clin_sample$Smoker_2fct[i] <- 0
      }
    }
  }
  clin_sample$Smoker_2fct <- factor(clin_sample$Smoker_2fct)
  levels(clin_sample$Smoker_2fct) <- c("no", "yes")

  # Smoker_3fct
  # Changing  variable from Rökstatus:
  # 1 = current smoker, 2 = Ex smoker, 3 = Never-smoker,  4 = Occasional smoker, 5 = Former occasional smoker
  # to
  # 0 = never, 1 = ex (including ex occasional), 2 = yes (including occasional)
  clin_sample$Smoker_3fct <- as.numeric(clin_sample$Smoker)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (is.na(clin_sample$Smoker_3fct[i]) == FALSE) {
      # set placeholders
      if (clin_sample$Smoker_3fct[i] == 1) {
        clin_sample$Smoker_3fct[i] <- 11
      }
      if (clin_sample$Smoker_3fct[i] == 2) {
        clin_sample$Smoker_3fct[i] <- 12
      }
      # run recode
      if (clin_sample$Smoker_3fct[i] == 3) {
        clin_sample$Smoker_3fct[i] <- 0
      }
      if (clin_sample$Smoker_3fct[i] == 4) {
        clin_sample$Smoker_3fct[i] <- 2
      }
      if (clin_sample$Smoker_3fct[i] == 5) {
        clin_sample$Smoker_3fct[i] <- 1
      }
      # recode placeholders
      if (clin_sample$Smoker_3fct[i] == 11) {
        clin_sample$Smoker_3fct[i] <- 2
      }
      if (clin_sample$Smoker_3fct[i] == 12) {
        clin_sample$Smoker_3fct[i] <- 1
      }
    }
  }
  clin_sample$Smoker_3fct <- factor(clin_sample$Smoker_3fct, ordered = TRUE)
  levels(clin_sample$Smoker_3fct) <- c("never", "ex", "yes")

  # Smoker_5fct
  # Changing  variable from Rökstatus:
  # 1 = current smoker, 2 = Ex smoker, 3 = Never-smoker,  4 = Occasional smoker, 5 = Former occasional smoker
  # to
  # 0 = never, 1 = ex occasional, 2 = ex, 3 = occasional smoker, 4 = yes
  clin_sample$Smoker_5fct <- as.numeric(clin_sample$Smoker)
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (is.na(clin_sample$Smoker_5fct[i]) == FALSE) {
      # run recode
      if (clin_sample$Smoker[i] == 1) {
        clin_sample$Smoker_5fct[i] <- 4
      }
      if (clin_sample$Smoker[i] == 2) {
        clin_sample$Smoker_5fct[i] <- 2
      }
      if (clin_sample$Smoker[i] == 3) {
        clin_sample$Smoker_5fct[i] <- 0
      }
      if (clin_sample$Smoker[i] == 4) {
        clin_sample$Smoker_5fct[i] <- 3
      }
      if (clin_sample$Smoker[i] == 5) {
        clin_sample$Smoker_5fct[i] <- 1
      }
    }
  }
  clin_sample$Smoker_5fct <- factor(clin_sample$Smoker_5fct, ordered = TRUE)
  levels(clin_sample$Smoker_5fct) <- c("never", "ex occ", "ex", "occ", "yes")

  # Systolic blood pressure
  # Blood pressure is measured when sitting in MONICA and lying down in VIP before 2009-09-01
  # All samples are from earlier than 2009-09-01
  # Therefore sbt and dbt have to recalculated for MONICA/VIP subjects
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (clin_sample[i, "Cohort"] == "MO") {
      if (clin_sample[i, "Sex"] == "female") {
        if (clin_sample[i, "Age"] <= 45) {
          clin_sample[i, "Sbt_VIP"] <- 8.669 + (0.919 * clin_sample[i, "Sbt_MO"]) #
          clin_sample[i, "Dbt_VIP"] <- 5.784 + (0.890 * clin_sample[i, "Dbt_MO"]) #
        }
        if (clin_sample[i, "Age"] > 45 & clin_sample[i, "Age"] <= 55) {
          clin_sample[i, "Sbt_VIP"] <- 16.051 + (0.859 * clin_sample[i, "Sbt_MO"]) #
          clin_sample[i, "Dbt_VIP"] <- 13.566 + (0.798 * clin_sample[i, "Dbt_MO"]) #
        }
        if (clin_sample[i, "Age"] > 55) {
          clin_sample[i, "Sbt_VIP"] <- 9.999 + (0.914 * clin_sample[i, "Sbt_MO"]) #
          clin_sample[i, "Dbt_VIP"] <- 7.992 + (0.870 * clin_sample[i, "Dbt_MO"]) #
        }
      }
      if (clin_sample[i, "Sex"] == "male") {
        if (clin_sample[i, "Age"] <= 45) {
          clin_sample[i, "Sbt_VIP"] <- 24.595 + (0.792 * clin_sample[i, "Sbt_MO"]) #
          clin_sample[i, "Dbt_VIP"] <- 17.282 + (0.753 * clin_sample[i, "Dbt_MO"]) #
        }
        if (clin_sample[i, "Age"] > 45 & clin_sample[i, "Age"] <= 55) {
          clin_sample[i, "Sbt_VIP"] <- 9.850 + (0.910 * clin_sample[i, "Sbt_MO"]) #
          clin_sample[i, "Dbt_VIP"] <- 12.363 + (0.812 * clin_sample[i, "Dbt_MO"]) #
        }
        if (clin_sample[i, "Age"] > 55) {
          clin_sample[i, "Sbt_VIP"] <- 7.763 + (0.936 * clin_sample[i, "Sbt_MO"]) #
          clin_sample[i, "Dbt_VIP"] <- 9.029 + (0.864 * clin_sample[i, "Dbt_MO"]) #
        }
      }
    }
  }
  for (i in 1:m) {
    if (clin_sample[i, "Cohort"] == "VIP") {
      if (clin_sample[i, "Sex"] == "female") {
        if (clin_sample[i, "Age"] <= 45) {
          clin_sample[i, "Sbt_MO"] <- 19.922 + (0.830 * clin_sample[i, "Sbt_VIP"]) #
          clin_sample[i, "Dbt_MO"] <- 13.680 + (0.847 * clin_sample[i, "Dbt_VIP"]) #
        }
        if (clin_sample[i, "Age"] > 45 & clin_sample[i, "Age"] <= 55) {
          clin_sample[i, "Sbt_MO"] <- 12.723 + (0.906 * clin_sample[i, "Sbt_VIP"]) #
          clin_sample[i, "Dbt_MO"] <- 17.675 + (0.800 * clin_sample[i, "Dbt_VIP"]) #
        }
        if (clin_sample[i, "Age"] > 55) {
          clin_sample[i, "Sbt_MO"] <- 13.817 + (0.900 * clin_sample[i, "Sbt_VIP"]) #
          clin_sample[i, "Dbt_MO"] <- 15.084 + (0.836 * clin_sample[i, "Dbt_VIP"]) #
        }
      }
      if (clin_sample[i, "Sex"] == "male") {
        if (clin_sample[i, "Age"] <= 45) {
          clin_sample[i, "Sbt_MO"] <- 21.612 + (0.835 * clin_sample[i, "Sbt_VIP"]) #
          clin_sample[i, "Dbt_MO"] <- 14.463 + (0.848 * clin_sample[i, "Dbt_VIP"]) #
        }
        if (clin_sample[i, "Age"] > 45 & clin_sample[i, "Age"] <= 55) {
          clin_sample[i, "Sbt_MO"] <- 19.748 + (0.861 * clin_sample[i, "Sbt_VIP"]) #
          clin_sample[i, "Dbt_MO"] <- 13.390 + (0.878 * clin_sample[i, "Dbt_VIP"]) #
        }
        if (clin_sample[i, "Age"] > 55) {
          clin_sample[i, "Sbt_MO"] <- 20.246 + (0.853 * clin_sample[i, "Sbt_VIP"]) #
          clin_sample[i, "Dbt_MO"] <- 16.308 + (0.833 * clin_sample[i, "Dbt_VIP"]) #
        }
      }
    }
  }

  # serum cholesterol and triglycerides
  # before 2009-09-01 VIP triglycerides and cholesterol were measured on reflotron instead of clinical chemistry lab
  # all MONICA samples were analysed by clinical chemistry
  # therefore, values have to be recalculated
  m <- nrow(clin_sample)
  for (i in 1:m) {
    if (clin_sample[i, "Cohort"] == "MO") {
      clin_sample[i, "Tg"] <- 0.177 + (0.932 * clin_sample[i, "stg_mo"])
      clin_sample[i, "Chol_tot"] <- 0.170 + (0.939 * clin_sample[i, "skol_mo"])
    }
  }

  # Check results
  return_list$check_var_cleaning <-
    clin_sample[clin_sample$Cohort == "MO", c(
      "Tg",
      "stg_mo",
      "Chol_tot",
      "skol_mo"
    )]


  # Sample date and SCD date. Parse as date variables.
  clin_sample$SCD_date <- as.Date(
    clin_sample$SCD_date,
    format = "%y-%m-%d"
  )
  clin_sample$Sample_date <- as.Date(
    clin_sample$Sample_date,
    format = "%y-%m-%d"
  )
  clin_sample$years_TO_scd <- clin_sample$SCD_date - clin_sample$Sample_date
  clin_sample$years_TO_scd <- clin_sample$years_TO_scd / 365.25
  clin_sample$years_AT_scd <- clin_sample$Age + (clin_sample$years_TO_scd)
  clin_sample$years_TO_scd <- round(as.numeric(clin_sample$years_TO_scd), 1)
  clin_sample$years_AT_scd <- round(as.numeric(clin_sample$years_AT_scd), 1)

  # Education_2fct. Changing name from utbildning to Education and
  # 1 = elementary, 2 = elementary + real, 3 = gymnasium, 4 = higher
  # to
  # 0 = elementary, 1 = real, gymnasium or higher (secondary)
  clin_sample$Education_2fct <- as.numeric(clin_sample$Education)
  for (i in seq(1, nrow(clin_sample))) {
    if (is.na(clin_sample$Education_2fct[i]) == FALSE) {
      # Run recode
      if (clin_sample$Education[i] == 1) {
        clin_sample$Education_2fct[i] <- 0
      }
      if (clin_sample$Education[i] == 2) {
        clin_sample$Education_2fct[i] <- 1
      }
      if (clin_sample$Education[i] == 3) {
        clin_sample$Education_2fct[i] <- 1
      }
      if (clin_sample$Education[i] == 4) {
        clin_sample$Education_2fct[i] <- 1
      }
    }
  }
  clin_sample$Education_2fct <- factor(clin_sample$Education_2fct)
  levels(clin_sample$Education_2fct) <- c("no", "yes")

  # Education_4fct. Changing name from utbildning to Education where
  # 1 = elementary, 2 = elementary + real, 3 = gymnasium, 4 = higher
  clin_sample$Education_4fct <- factor(clin_sample$Education, ordered = TRUE)
  levels(clin_sample$Education_4fct) <- c("elem", "real", "gymn", "high")

  # Make lab_diabetes, using Glc_0h and Glc_2h to complement questionnaire diabetes
  clin_sample$lab_diabetes_hi0 <- ifelse(
    clin_sample$Glc_0h >= 7.0 & !is.na(clin_sample$Glc_0h),
    1, 0
  )
  clin_sample$lab_diabetes_hi2 <- ifelse(
    clin_sample$Glc_2h >= 11.1 & !is.na(clin_sample$Glc_2h),
    1, 0
  )
  clin_sample$lab_diabetes_fp <- ifelse(
    clin_sample$lab_diabetes_hi0 == 1 & clin_sample$lab_diabetes_hi2 == 0 &
      !is.na(clin_sample$Glc_0h) & !is.na(clin_sample$Glc_2h),
    1, 0
  )
  clin_sample$lab_diabetes_tp <- ifelse(
    (clin_sample$lab_diabetes_hi0 == 1 | clin_sample$lab_diabetes_hi2 == 1) &
      clin_sample$lab_diabetes_fp == 0,
    1, 0
  )
  clin_sample$lab_diabetes_tpq <- ifelse(
    clin_sample$lab_diabetes_tp == 1 | (clin_sample$Diabetes_Q == "yes" &
      !is.na(clin_sample$Diabetes_Q)),
    1, 0
  )
  clin_sample$lab_diabetes_tpq <- as.factor(clin_sample$lab_diabetes_tpq)
  levels(clin_sample$lab_diabetes_tpq) <- c("no", "yes")

  # SCD_timetodeath
  # Setting Ctrl NA's to control
  clin_sample$SCD_timetodeath <- as.character(clin_sample$SCD_timetodeath)
  clin_sample$SCD_timetodeath[is.na(clin_sample$SCD_timetodeath) & clin_sample$Case_control == "Ctrl"] <- "Ctrl"
  clin_sample$SCD_timetodeath[is.na(clin_sample$SCD_timetodeath) & clin_sample$Case_control == "Case"] <- "Unknown"
  clin_sample$SCD_timetodeath <- as.factor(clin_sample$SCD_timetodeath)
  levels(clin_sample$SCD_timetodeath) <- c("<1h", "1-24h", "ctrl", "unknown")

  #----------------------------------------------------------------------------
  # Pharmaceutical data
  #----------------------------------------------------------------------------
  # Because of different codings in the study database, questionnaire data about
  # medicines have to be extracted differently.

  med_data <- read.csv2(path_sample,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = rep("character", 113),
    na.strings = c("")
  )

  vars1 <- c(
    "pat_code",
    "fallkontroll_fia3",
    "med_C5a",
    "med_C5b",
    "med_C5c",
    "med_C5d",
    "med_C5e",
    "med_C5f",
    "smartmed",
    "med_asahjar",
    "med_asacvs",
    "diabetesbehandling_a",
    "diabetesbehandling_b",
    "diabetesbehandling_c",
    "diabetesbehandling_d"
  )

  med_data <- med_data[vars1]

  table(med_data[, "med_C5a"])
  table(med_data[, "med_C5b"])
  table(med_data[, "med_C5c"])
  table(med_data[, "med_C5d"])
  table(med_data[, "med_C5e"])
  table(med_data[, "med_C5f"])
  table(med_data[, "smartmed"])
  table(med_data[, "med_asahjar"])
  table(med_data[, "med_asacvs"])

  # For these questionnaires, "yes" answers were set as a number, eg "1",
  # and you could only answer yes, or else the variable is set as NA.
  # The NA coding is however the same as before. Therefore, they have to be
  # preprocessed separately.
  
  m <- nrow(med_data)
  n <- ncol(med_data)

  for (j in 1:n) {
    for (i in 1:m) {
      if (is.na(med_data[i, j]) == TRUE) {
        med_data[i, j] <- "0"
      }
      if (med_data[i, j] %in% c(5555, 6666, 7777, 8888, 9999)) {
        med_data[i, j] <- NA
      }
    }
  }

  # Blood pressure medication
  table(med_data[, "med_C5a"])

  med_data$med_C5a <- as.numeric(med_data$med_C5a)
  table(med_data[, "med_C5a"])

  # In VIP you can only answer 1=yes to this question. 
  # However, in monica you can also answer 2 = no, and 3 = unsure.
  # However, no one in MONICA was unsure about taking BP lowering drugs
  # Therefore, 2 is set to 1
  med_data$med_C5a[med_data$med_C5a == 2 & !is.na(med_data$med_C5a)] <- 0
  table(med_data[, "med_C5a"])

  med_data$med_C5a <- as.factor(med_data$med_C5a)
  table(med_data[, "med_C5a"])
  levels(med_data$med_C5a) <- c("no", "yes")
  table(med_data[, "med_C5a"])

  # Heart/angina medication (only VIP)
  table(med_data[, "med_C5b"])
  med_data$med_C5b <- as.factor(med_data$med_C5b)
  table(med_data[, "med_C5b"])
  levels(med_data$med_C5b) <- c("no", "yes")
  table(med_data[, "med_C5b"])

  # BZ/antihist medication (only VIP)
  table(med_data[, "med_C5c"])
  med_data$med_C5c <- as.factor(med_data$med_C5c)
  table(med_data[, "med_C5c"])
  levels(med_data$med_C5c) <- c("no", "yes")
  table(med_data[, "med_C5c"])

  # PPI-like (only VIP)
  table(med_data[, "med_C5d"])
  med_data$med_C5d <- as.factor(med_data$med_C5d)
  table(med_data[, "med_C5d"])
  levels(med_data$med_C5d) <- c("no", "yes")
  table(med_data[, "med_C5d"])

  # Lipid lowering
  table(med_data[, "med_C5e"])
  med_data$med_C5e <- as.numeric(med_data$med_C5e)
  table(med_data[, "med_C5e"])

  # In VIP you can only answer 1 = yes to this question. 
  # However, in monica you can also answer 2=no, and 3=unsure.
  # Therefore, 2 is set to 0
  med_data$med_C5e[med_data$med_C5e == 2 & !is.na(med_data$med_C5e)] <- 0
  # Two people were unsure about taking lipid lowering drugs
  # They are set to NA
  med_data$med_C5e[med_data$med_C5e == 3 & !is.na(med_data$med_C5e)] <- NA

  table(med_data[, "med_C5e"])
  med_data$med_C5e <- as.factor(med_data$med_C5e)
  table(med_data[, "med_C5e"])

  levels(med_data$med_C5e) <- c("no", "yes")
  table(med_data[, "med_C5e"])

  # pain medicine (relatively new variable in VIP, many obs missing)
  table(med_data[, "smartmed"])
  med_data$smartmed <- as.factor(med_data$smartmed)
  table(med_data[, "smartmed"])
  levels(med_data$smartmed) <- c("no", "yes")

  # MONICA ASA for myocardial infarction (only MO), 
  table(med_data[, "med_asahjar"])
  med_data$med_asahjar <- as.factor(med_data$med_asahjar)
  levels(med_data$med_asahjar) <- c("yes", "yother", "no")
  

  # MO ASA for stroke (only MO)
  table(med_data[, "med_asacvs"])
  med_data$med_asacvs <- as.factor(med_data$med_asacvs)
  levels(med_data$med_asacvs) <- c("yes", "yother", "no")
  
  # Add asa to VIP Heart_drug
  
  for (i in seq(med_data$med_C5b)) {
    if (!is.na(med_data$med_asahjar[i])) {
      if (med_data$med_asahjar[i] == "yes") {
        med_data$med_C5b[i] <- "yes"
      }
      if (med_data$med_asahjar[i] == "yother") {
        med_data$med_C5b[i] <- "yes"
      }
      if (med_data$med_asahjar[i] == "no") {
        med_data$med_C5b[i] <- "no"
      }
    }
  }
  
  for (i in seq(med_data$med_C5b)) {
    if (!is.na(med_data$med_asacvs[i])) {
      if (med_data$med_asacvs[i] == "yes" && med_data$med_C5b[i] != "yes") {
        med_data$med_C5b[i] <- "yes"
      }
      if (med_data$med_asacvs[i] == "yother" && med_data$med_C5b[i] != "yes") {
        med_data$med_C5b[i] <- "yes"
      }
      if (med_data$med_asacvs[i] == "no" && med_data$med_C5b[i] != "yes") {
        med_data$med_C5b[i] <- "no"
      }
    }
  }

  # diet and physical exercise treatment of diabetes
  table(med_data[, "diabetesbehandling_a"])
  med_data$diabetesbehandling_a <- as.factor(med_data$diabetesbehandling_a)
  levels(med_data$diabetesbehandling_a) <- c("no", "yes")

  # pill treatment (that were used from the 80's - 20's) of diabetes
  table(med_data[, "diabetesbehandling_b"])
  med_data$diabetesbehandling_b <- as.factor(med_data$diabetesbehandling_b)
  levels(med_data$diabetesbehandling_b) <- c("no", "yes")

  # insulin treatment
  table(med_data[, "diabetesbehandling_c"])
  med_data$diabetesbehandling_c <- as.factor(med_data$diabetesbehandling_c)
  levels(med_data$diabetesbehandling_c) <- c("no", "yes")

  # no current treatment of the diabetes, 1=yes (follow-up question if you answered yes on Diabetes_Q)
  table(med_data[, "diabetesbehandling_d"])
  med_data$diabetesbehandling_d <- as.factor(med_data$diabetesbehandling_d)
  levels(med_data$diabetesbehandling_d) <- c("no", "yes")
  
  # Rename and merge with clin_sample
  names(med_data) <- clin_meta$name2[clin_meta$name1 %in% names(med_data)]
  
  # Remove duplicates for med_data as well
  med_data <- subset(
    med_data, 
    (Patient_code %in% check_table_names & Case_control == "F") |
      !Patient_code %in% check_table_names
  )
  
  med_data <- med_data[!names(med_data) %in% "Case_control"]
  
  med_inc <- c("BP_drug",
               "Heart_drug",
               "BZoHist_drug",
               "PPI_drug",
               "Lipid_drug",
               "No_drug",
               "smartmed",
               "med_asahjar",
               "med_asacvs",
               "Diabet_diet",
               "Diabet_pill",
               "Diabet_insulin",
               "diabetesbehandling_d")
  
  clin_sample <- clin_sample[!names(clin_sample) %in% med_inc] #rm old med vars
  # from clin_sample
  
  clin_sample <- merge(clin_sample, med_data, by = "Patient_code")

  #----------------------------------------------------------------------------
  #  Filter elevant variables
  #----------------------------------------------------------------------------
  
  # Out varcleaned data before subsetting
  return_list$varcleaned_data <- clin_sample
  
  #Subset
  clin_sample <- clin_sample[c(
    "Patient_code",
    "Set_FIA3",
    "FIAnum_FIA3",
    "SCD_date",
    "Sample_date",
    "years_TO_scd",
    "years_AT_scd",
    "Sex",
    "Chol_tot",
    "Age",
    "Case_control",
    "ApoB100",
    "ApoA1",
    "BMI",
    "Sbt_VIP",
    "Dbt_VIP",
    "Smoker_2fct",
    "Smoker_3fct",
    "Smoker_5fct",
    "Glc_0h",
    "Glc_2h",
    "Diabetes_Q",
    "lab_diabetes_tpq",
    "Fast_sample",
    "case_freeze_thawed",
    "Education_2fct",
    "Education_4fct",
    "SCD_type",
    "SCD_timetodeath",
    "BP_drug",
    "Heart_drug",
    "BZoHist_drug",
    "PPI_drug",
    "Lipid_drug",
    "Diabet_diet",
    "Diabet_pill",
    "Diabet_insulin",
    "Lpa",
    "CRP"
  )]
  return_list$data_out <- clin_sample

  #----------------------------------------------------------------------------
  # FIN RETURN
  #----------------------------------------------------------------------------
  return(return_list)
}
