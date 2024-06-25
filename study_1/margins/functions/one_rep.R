one_rep <- function(repNr, condition) {
  source("functions/generate_model.R", local = TRUE)
  source("functions/generate_data.R", local = TRUE)
  
  set.seed(repNr + 2023)
  
  Sigma_nT <- generate_model(beta_z = condition$beta_z_nT, load_z = condition$load_z)
  Sigma_wT <- generate_model(beta_z = condition$beta_z_wT, load_z = condition$load_z)
  
  num_z_vars <- sum(grepl("^z", colnames(Sigma_nT)))
  
  generate_cop_data <- function(treatment) {
    generate_data(
      n = condition$n, 
      mean_y = ifelse(treatment, condition$treat_eff, 0),
      mean_z = ifelse(treatment, condition$mean_z_wT, condition$mean_z_nT),
      Sigma = ifelse(treatment, Sigma_wT, Sigma_nT),
      margin_dist_z = ifelse(treatment, condition$marginsZ_wT, condition$marginsZ_nT),
      marginsZ_rate = ifelse(treatment, condition$marginsZ_theta_wT, condition$marginsZ_theta_nT),
      margin_dist_y = ifelse(treatment, condition$marginsY_wT, condition$marginsY_nT),
      copula_type_ZY = ifelse(treatment, condition$copulaZY_wT, condition$copulaZY_nT),
      Nmax = 10^5,
      num_z_vars = num_z_vars
    )
  }
  
  cop_data_nT <- generate_cop_data(FALSE)
  cop_data_wT <- generate_cop_data(TRUE)
  
  samples_nT <- cop_data_nT[[1]] %>% mutate(x = 0)
  samples_wT <- cop_data_wT[[1]] %>% mutate(x = 1)
  
  discretize_data <- function(data, thresholds) {
    cut(data, thresholds, labels = FALSE)
  }
  
  discretize_columns <- function(sample, pattern, thresholds) {
    cols_to_discretize <- grep(pattern, names(sample))
    sample[cols_to_discretize] <- lapply(sample[cols_to_discretize], discretize_data, thresholds)
    sample
  }
  
  sample_compl_disc_nT <- discretize_columns(samples_nT, "^z", condition$thresholds[[1]])
  sample_compl_disc_wT <- discretize_columns(samples_wT, "^z", condition$thresholds[[1]])
  
  sample_compl <- rbind(samples_nT, samples_wT)
  sample_compl_disc <- rbind(sample_compl_disc_nT, sample_compl_disc_wT)
  
  sample_compl$x <- factor(sample_compl$x)
  sample_compl_disc$x <- factor(sample_compl_disc$x)
  
  z_cols <- grep("^z", names(sample_compl_disc), value = TRUE)
  sample_compl_disc[z_cols] <- lapply(sample_compl_disc[z_cols], as.ordered)
  
  sample_compl_disc_cont <- sample_compl_disc %>%
    mutate(across(all_of(z_cols), ~ as.numeric(.) - mean(as.numeric(.))))
  
  mm_ML <- generateMeasurementModel(
    names = "xi",
    indicators = list("xi" = paste0("z", 1:num_z_vars)),
    ncells = 2,
    model = c("tau-cong"),
    data = as.data.frame(sample_compl_disc_cont)
  ) %>%
    str_replace("c\\(1,1\\)\\*z1", sprintf("c(%s,%s)*z1", condition$load_z, condition$load_z)) %>%
    str_replace("z1 ~ c\\(0,0\\)\\*1", sprintf("z1 ~ c(%s,%s)*1", condition$mean_z_nT, condition$mean_z_wT))
  
  run_effectLite <- function(data, mm, estimator) {
    effectLite(y = "y", x = "x", z = "xi",
               data = as.data.frame(data),
               method = "sem",
               measurement = mm,
               interactions = "all",
               estimator = estimator)
  }
  
  model_elr_cont <- run_effectLite(sample_compl_disc_cont, mm_ML, "ML")
  model_elr_cont_robust <- run_effectLite(sample_compl_disc_cont, mm_ML, "MLR")
  
  mm_disc <- generateMeasurementModel(
    names = "xi",
    indicators = list("xi" = paste0("z", 1:num_z_vars)),
    ncells = 2,
    model = c("tau-cong-categorical"),
    data = as.data.frame(sample_compl_disc)
  ) %>%
    str_replace("xi =~ c\\(1,1\\)\\*z1", sprintf("xi =~ c(%s,%s)*z1", condition$load_z, condition$load_z)) %>%
    str_replace_all("z\\d+ ~ c\\(0,0\\)\\*1", sprintf("\\0 ~ c(%s,%s)*1", condition$mean_z_nT, condition$mean_z_wT)) %>%
    str_replace("z1 \\| c\\(0,0\\)\\*t1", sprintf("z1 | c(%s,%s)*t1", condition$thresholds[[1]][2], condition$thresholds[[1]][2]))
  
  model_elr_disc <- run_effectLite(sample_compl_disc, mm_disc, "DWLS")
  
  list(
    ML = model_elr_cont,
    MLR = model_elr_cont_robust,
    DWLS = model_elr_disc,
    cov_matrix_pop = list("Sigma_nT" = Sigma_nT, "Sigma_wT" = Sigma_wT),
    cop_data = list("vine_cop_nT" = cop_data_nT[[2]], "vine_cop_wT" = cop_data_wT[[2]]),
    Method = c("model_ml", "model_mlr", "model_dwls"),
    replication = repNr,
    Seed = .Random.seed
  )
}