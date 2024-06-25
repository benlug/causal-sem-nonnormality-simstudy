one_rep <- function(repNr, condition) {
  source("functions/generate_model.R", local = TRUE)
  source("functions/generate_data.R", local = TRUE)
  
  set.seed(repNr + 2023)
  
  generate_sigma <- function(group) {
    generate_model(
      sMZ = condition[[paste0("sMZG", group)]],
      sMW = condition[[paste0("sMWG", group)]],
      sYW = condition[[paste0("sYWG", group)]],
      sYM = condition[[paste0("sYMG", group)]],
      sYZ = condition[[paste0("sYZG", group)]],
      MpreYpre_corr = condition$MpreYpre_corr,
      loading = condition$loading
    )
  }
  
  Sigma_nT <- generate_sigma(1)
  Sigma_wT <- generate_sigma(2)
  
  num_m_vars <- sum(grepl("^m", colnames(Sigma_nT)))
  
  generate_cop_data <- function(group) {
    condition$mean_M_G2_latent <- ifelse(group == 2, condition$loading * condition$mean_M_G2, condition$mean_M_G1)
    generate_data(
      n = condition$n,
      mean_Ypre = condition[[paste0("mean_Ypre_G", group)]],
      mean_Y = condition[[paste0("mean_Y_G", group)]],
      mean_Mpre = condition[[paste0("mean_Mpre_G", group)]],
      mean_M = condition$mean_M_G2_latent,
      Sigma = ifelse(group == 1, Sigma_nT, Sigma_wT),
      margin_dist_M = condition[[paste0("marginsM_", ifelse(group == 1, "nT", "wT"))]],
      margin_dist_Ypre = condition[[paste0("marginsYpre_", ifelse(group == 1, "nT", "wT"))]],
      copula = condition[[paste0("copula_", ifelse(group == 1, "nT", "wT"))]],
      Nmax = 10^5,
      num_m_vars = num_m_vars
    )
  }
  
  cop_data_nT <- generate_cop_data(1)
  cop_data_wT <- generate_cop_data(2)
  
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
  
  sample_compl_disc_nT <- discretize_columns(samples_nT, "^m", condition$thresholds[[1]])
  sample_compl_disc_wT <- discretize_columns(samples_wT, "^m", condition$thresholds[[1]])
  
  sample_compl <- rbind(samples_nT, samples_wT)
  sample_compl_disc <- rbind(sample_compl_disc_nT, sample_compl_disc_wT) %>%
    mutate(x = factor(x))
  
  m_cols <- grep("^m", names(sample_compl_disc), value = TRUE)
  sample_compl_disc[m_cols] <- lapply(sample_compl_disc[m_cols], as.ordered)
  
  sample_compl_disc_cont <- sample_compl_disc %>%
    mutate(across(all_of(m_cols), as.numeric))
  
  model <- '
    M ~ c(sMZG1,sMZG2)*Mpre + c(sMWG1,sMWG2)*Ypre
    Y ~ c(sYWG1,sYWG2)*Ypre + c(sYMG1,sYMG2)*M + c(sYZG1,sYZG2)*Mpre
    Mpre ~~ Ypre
    
    Mpre ~ c(mZG1,mZG2)*1
    Ypre ~ c(mWG1,mWG2)*1
    M ~ c(iMG1,iMG2)*1
    Y ~ c(iYG1,iYG2)*1
    
    # Calculation of New Parameters
    mMG1 := iMG1 + sMZG1*mZG1 + sMWG1*mWG1
    mMG2 := iMG2 + sMZG2*mZG2 + sMWG2*mWG2              
    mZ := (1*mZG1 + 1*mZG2)/2
    mW := (1*mWG1 + 1*mWG2)/2
    mM := (1*mMG1 + 1*mMG2)/2
    mMadjG1 := iMG1 + sMZG1*mZ + sMWG1*mW
    mMadjG2 := iMG2 + sMZG2*mZ + sMWG2*mW              
    mYadjG1 := iYG1 + sYWG1*mW + sYZG1*mZ + sYMG1*mM
    mYadjG2 := iYG2 + sYWG2*mW + sYZG2*mZ + sYMG2*mM
    iYrG1 := iYG1 + sYMG1*iMG1
    iYrG2 := iYG2 + sYMG2*iMG2
    sYZrG1 := sYZG1 + sYMG1*sMZG1
    sYZrG2 := sYZG2 + sYMG2*sMZG2
    sYWrG1 := sYWG1 + sYMG1*sMWG1
    sYWrG2 := sYWG2 + sYMG2*sMWG2
    mYadtG1 := iYrG1 + sYWrG1*mW + sYZrG1*mZ
    mYadtG2 := iYrG2 + sYWrG2*mW + sYZrG2*mZ
    direct := mYadjG2 - mYadjG1             
    total := mYadtG2 - mYadtG1
    indirect := total - direct
    ga10 := iYG2 - iYG1
    ga11 := sYWG2 - sYWG1
    ga12 := sYZG2 - sYZG1
    ga13 := sYMG2 - sYMG1
    nde0 := ga10 + ga11*mW + ga12*mZ + ga13*mMadjG1
    nde1 := ga10 + ga11*mW + ga12*mZ + ga13*mMadjG2
    adet := ga10 + ga11*mWG2 + ga12*mZG2 + ga13*mMG2
  '
  
  mm_ML <- generateMeasurementModel(
    names = "M",
    indicators = list("M" = paste0("m", 1:num_m_vars)),
    ncells = 2,
    model = c("tau-cong"),
    data = as.data.frame(sample_compl_disc_cont)
  ) %>%
    str_replace("c\\(1,1\\)\\*m1", sprintf("c(%s,%s)*m1", condition$loading, condition$loading)) %>%
    str_replace("m1 ~ c\\(0,0\\)\\*1", sprintf("m1 ~ c(%s,%s)*1", condition$mean_Mpre_G1, condition$mean_Mpre_G2))
  
  mm_ML <- paste(mm_ML, model)
  
  run_sem <- function(mm, data, estimator) {
    sem(mm, data = data, group = "x", estimator = estimator)
  }
  
  model_ML <- run_sem(mm_ML, sample_compl_disc_cont, "ML")
  model_MLR <- run_sem(mm_ML, sample_compl_disc_cont, "MLR")
  
  mm_disc <- generateMeasurementModel(
    names = "M",
    indicators = list("M" = paste0("m", 1:num_m_vars)),
    ncells = 2,
    model = c("tau-cong-categorical"),
    data = as.data.frame(sample_compl_disc)
  ) %>%
    str_replace("M =~ c\\(1,1\\)\\*m1", sprintf("M =~ c(%s,%s)*m1", condition$loading, condition$loading)) %>%
    str_replace_all("m\\d+ ~ c\\(0,0\\)\\*1", sprintf("\\0 ~ c(%s,%s)*1", condition$mean_Mpre_G1, condition$mean_Mpre_G2)) %>%
    str_replace("m1 \\| c\\(0,0\\)\\*t1", sprintf("m1 | c(%s,%s)*t1", condition$thresholds[[1]][2], condition$thresholds[[1]][2]))
  
  mm_disc <- paste(mm_disc, model)
  
  model_DWLS <- run_sem(mm_disc, sample_compl_disc, "DWLS")
  
  list(
    ML = model_ML,
    MLR = model_MLR,
    DWLS = model_DWLS,
    cov_matrix_pop = list("Sigma_nT" = Sigma_nT, "Sigma_wT" = Sigma_wT),
    cop_data = list("vine_cop_nT" = cop_data_nT[[2]], "vine_cop_wT" = cop_data_wT[[2]]),
    Method = c("model_ml", "model_mlr", "model_dwls"),
    replication = repNr,
    Seed = .Random.seed
  )
}