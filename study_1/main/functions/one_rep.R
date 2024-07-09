one_rep <- function(repNr, condition) {
  
  source("functions/generate_model.R", local = TRUE)
  source("functions/generate_data.R", local = TRUE)
  
  # generate random seed based on replication number
  set.seed(repNr + 2023)
  
  # target covariance matrices for x=0 and x=1
  Sigma_nT <- generate_model(beta_z = condition$beta_z_nT, load_z = condition$load_z)
  Sigma_wT <- generate_model(beta_z = condition$beta_z_wT, load_z = condition$load_z)
  
  num_z_vars <- sum(grepl("^z", colnames(Sigma_nT))) 
  
  cop_data_nT <- generate_data(n = condition$n, 
                               mean_y = 0, 
                               mean_z = condition$mean_z_nT,
                               Sigma = Sigma_nT, 
                               margin_dist_z = condition$marginsZ_nT,
                               margin_dist_y = condition$marginsY_nT,
                               copula_type_ZY = condition$copulaZY_nT,
                               Nmax = 10^5,
                               num_z_vars = num_z_vars)
  cop_data_wT <- generate_data(n = condition$n, 
                               mean_y = condition$treat_eff,
                               mean_z = condition$mean_z_wT,
                               Sigma = Sigma_wT, 
                               margin_dist_z = condition$marginsZ_wT, 
                               margin_dist_y = condition$marginsY_wT,
                               copula_type_ZY = condition$copulaZY_wT,
                               Nmax = 10^5,
                               num_z_vars = num_z_vars)
  
  samples_nT <- cop_data_nT[[1]]
  samples_wT <- cop_data_wT[[1]]
  
  # covariance bias estimate
  cov_matrix_est_bias_nT <- Sigma_nT - cov(samples_nT)
  cov_matrix_est_bias_wT <- Sigma_wT - cov(samples_wT)
  
  # add treatment variable
  samples_nT <- mutate(samples_nT, x = 0)
  samples_wT <- mutate(samples_wT, x = 1)
  
  discretize_data <- function(data, thresholds) {
    # This function discretizes continuous data into ordered categories.
    # Each category represents an interval of values determined by 'thresholds'. 
    disc_data <- cut(data, thresholds, labels = FALSE)
    return(disc_data)
  }
  
  discretize_columns <- function(sample, pattern, thresholds) {
    # This function finds columns in a data frame that match a specified pattern,
    # and applies the 'discretize_data' function to them, using a given set of thresholds. 
    cols_to_discretize <- grepl(pattern, colnames(sample))
    sample[, cols_to_discretize] <- lapply(sample[, cols_to_discretize], discretize_data, thresholds)
    return(sample)
  }
  
  sample_compl_disc_nT <- discretize_columns(samples_nT, "z", condition$thresholds[[1]])
  sample_compl_disc_wT <- discretize_columns(samples_wT, "z", condition$thresholds[[1]])
  
  # combine continuous and categorical tibbles with and without treatment
  sample_compl <- rbind(samples_nT, samples_wT)
  sample_compl_disc <- rbind(sample_compl_disc_nT, sample_compl_disc_wT)
  
  # make treatment variable a factor
  sample_compl$x <- as_factor(sample_compl$x) 
  sample_compl_disc$x <- as_factor(sample_compl_disc$x)
  
  # get column names that start with 'y'
  z_cols <- grep("^z", names(sample_compl_disc), value = TRUE)
  # make y columns ordered
  for (col in z_cols) {
    sample_compl_disc[[col]] <- as.ordered(sample_compl_disc[[col]])
  }
  # create a copy of the discretized data, but with y columns as numeric
  sample_compl_disc_cont <- sample_compl_disc
  for (col in z_cols) {
    sample_compl_disc_cont[[col]] <- as.numeric(sample_compl_disc_cont[[col]])
  }
  
  sample_compl_disc_cont <- sample_compl_disc_cont |> 
    mutate(across(c(z1, z2, z3), ~ . - mean(.)))
  
  mm_ML <- generateMeasurementModel(names = "xi",
                                    indicators = list("xi" = paste0("z", 1:num_z_vars)),
                                    ncells = 2,
                                    model = c("tau-cong"),
                                    data = as.data.frame(sample_compl_disc_cont))
  mm_ML <- gsub("c\\(1,1\\)\\*z1", 
                paste0("c\\(", condition$load_z, ",", condition$load_z, "\\)*z1"), 
                mm_ML)
  mm_ML <- gsub("z1 ~ c\\(0,0\\)\\*1", 
                paste0("z1 ~ c\\(", condition$mean_z_nT, ",", condition$mean_z_wT, "\\)*1"), 
                mm_ML)
  
  model_elr_cont <- effectLite(y = "y", x = "x", z = "xi",
                               data = as.data.frame(sample_compl_disc_cont),
                               method = "sem",
                               measurement = mm_ML,
                               interactions = "all",
                               estimator = "ML")
  model_elr_cont_robust <- effectLite(y = "y", x = "x", z = "xi",
                                      data = as.data.frame(sample_compl_disc_cont),
                                      method = "sem",
                                      measurement = mm_ML,
                                      interactions = "all",
                                      estimator = "MLR")
  
  mm_disc <- generateMeasurementModel(names = "xi",
                                      indicators = list("xi" = paste0("z", 1:num_z_vars)),
                                      ncells = 2,
                                      model = c("tau-cong-categorical"),
                                      data = as.data.frame(sample_compl_disc))
  mm_disc <- gsub("xi =~ c\\(1,1\\)\\*z1", 
                  paste0("xi =~ c(", condition$load_z, ",", condition$load_z, ")*z1"), 
                  mm_disc)
  for (i in 1:num_z_vars) {
    mm_disc <- gsub(paste0("z", i, " ~ c\\(0,0\\)\\*1"), 
                    paste0("z", i, " ~ c(", condition$mean_z_nT, ",", condition$mean_z_wT, ")*1"), 
                    mm_disc)
  }
  mm_disc <- gsub("z1 \\| c\\(0,0\\)\\*t1", 
                  paste0("z1 | c(", condition$thresholds[[1]][2], ",", 
                         condition$thresholds[[1]][2], ")*t1"), 
                  mm_disc)
  model_elr_disc <- effectLite(y = "y", x = "x", z = "xi",
                               data = as.data.frame(sample_compl_disc),
                               method = "sem",
                               measurement = mm_disc,
                               interactions = "all",
                               estimator = "DWLS")

  res <- list(ML = model_elr_cont,
              MLR = model_elr_cont_robust,
              DWLS = model_elr_disc,
              cov_matrix_pop = list("Sigma_nT" = Sigma_nT, 
                                    "Sigma_wT" = Sigma_wT),
              cop_data = list("vine_cop_nT" = cop_data_nT[[2]], 
                              "vine_cop_wT" = cop_data_wT[[2]]),
              Method = c("model_ml", "model_mlr", "model_dwls"),
              replication = repNr,
              Seed = .Random.seed)
  
  return(res)
  
}
