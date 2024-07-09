one_rep <- function(repNr, condition) {
  
  source("functions/generate_model.R", local = TRUE)
  source("functions/generate_data.R", local = TRUE)
  
  # generate random seed based on replication number
  set.seed(repNr + 2023)
  
  # target covariance matrices for x=0 and x=1
  Sigma_nT <- generate_model(sMZ = condition$sMZG1,
                             sMW = condition$sMWG1,
                             sYW = condition$sYWG1,
                             sYM = condition$sYMG1, 
                             sYZ = condition$sYZG1,
                             MpreYpre_corr = condition$MpreYpre_corr,
                             loading = condition$loading)
  Sigma_wT <- generate_model(sMZ = condition$sMZG2,
                             sMW = condition$sMWG2,
                             sYW = condition$sYWG2,
                             sYM = condition$sYMG2,
                             sYZ = condition$sYZG2,
                             MpreYpre_corr = condition$MpreYpre_corr,
                             loading = condition$loading)
  
  num_m_vars <- sum(grepl("^m", colnames(Sigma_nT)))
  
  cop_data_nT <- generate_data(n = condition$n, 
                               mean_Ypre = condition$mean_Ypre_G1,
                               mean_Y = condition$mean_Y_G1,
                               mean_Mpre = condition$mean_Mpre_G1,
                               mean_M = condition$mean_M_G1,
                               Sigma = Sigma_nT, 
                               margin_dist_M = condition$marginsM_nT,
                               margin_dist_Ypre = condition$marginsYpre_nT,
                               copula = condition$copula_nT,
                               Nmax = 10^5,
                               num_m_vars = num_m_vars)
  cop_data_wT <- generate_data(n = condition$n, 
                               mean_Ypre = condition$mean_Ypre_G2,
                               mean_Y = condition$mean_Y_G2,
                               mean_Mpre = condition$mean_Mpre_G2,
                               mean_M = condition$mean_M_G2,
                               Sigma = Sigma_wT, 
                               margin_dist_M = condition$marginsM_wT,
                               margin_dist_Ypre = condition$marginsYpre_wT,
                               copula = condition$copula_wT,
                               Nmax = 10^5,
                               num_m_vars = num_m_vars)
  
  samples_nT <- cop_data_nT[[1]] |> 
    mutate(x = 0) 
    
  samples_wT <- cop_data_wT[[1]] |> 
    mutate(x = 1)
  
  # covariance bias estimate
  # cov_matrix_est_bias_nT <- Sigma_nT - cov(samples_nT)
  # cov_matrix_est_bias_wT <- Sigma_wT - cov(samples_wT)
  
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
  
  sample_compl_disc_nT <- discretize_columns(samples_nT, "m", condition$thresholds[[1]])
  sample_compl_disc_wT <- discretize_columns(samples_wT, "m", condition$thresholds[[1]])
  
  sample_compl <- rbind(samples_nT, samples_wT)
  sample_compl_disc <- rbind(sample_compl_disc_nT, sample_compl_disc_wT)
  sample_compl$x <- as_factor(sample_compl$x) 
  sample_compl_disc$x <- as_factor(sample_compl_disc$x)
  
  m_cols <- grep("^m", names(sample_compl_disc), value = TRUE)
  for (col in m_cols) {
    sample_compl_disc[[col]] <- as.ordered(sample_compl_disc[[col]])
  }
  sample_compl_disc_cont <- sample_compl_disc
  for (col in m_cols) {
    sample_compl_disc_cont[[col]] <- as.numeric(sample_compl_disc_cont[[col]])
  }
  
  model <- '
  
  M ~ c(sMZG1,sMZG2)*Mpre + c(sMWG1,sMWG2)*Ypre
  Y ~ c(sYWG1,sYWG2)*Ypre + c(sYMG1,sYMG2)*M + c(sYZG1,sYZG2)*Mpre
  Mpre ~~ Ypre
  
  Mpre ~ c(mZG1,mZG2)*1
  Ypre ~ c(mWG1,mWG2)*1
  M ~ c(iMG1,iMG2)*1
  Y ~ c(iYG1,iYG2)*1
  
  ## Calculation of New Parameters
  
  # Conditional means of M
  mMG1 := iMG1 + sMZG1*mZG1 + sMWG1*mWG1
  mMG2 := iMG2 + sMZG2*mZG2 + sMWG2*mWG2              
  
  # Unconditional means of Mpre, Ypre and M              
  mZ := (1*mZG1 + 1*mZG2)/2
  mW := (1*mWG1 + 1*mWG2)/2
  mM := (1*mMG1 + 1*mMG2)/2
  
  # Adjusted means of M
  mMadjG1 := iMG1 + sMZG1*mZ + sMWG1*mW
  mMadjG2 := iMG2 + sMZG2*mZ + sMWG2*mW              
  
  # Adjusted means of the outcome variable Y in both groups
  mYadjG1 := iYG1 + sYWG1*mW + sYZG1*mZ + sYMG1*mM
  mYadjG2 := iYG2 + sYWG2*mW + sYZG2*mZ + sYMG2*mM
  
  # Group specific regressions Y on Mpre and Ypre  
  iYrG1 := iYG1 + sYMG1*iMG1
  iYrG2 := iYG2 + sYMG2*iMG2
  
  sYZrG1 := sYZG1 + sYMG1*sMZG1
  sYZrG2 := sYZG2 + sYMG2*sMZG2
  
  sYWrG1 := sYWG1 + sYMG1*sMWG1
  sYWrG2 := sYWG2 + sYMG2*sMWG2
  
  # Adjusted means for the computation of the total effect
  # based on group specific regressions Y on Mpre and Ypre
  mYadtG1 := iYrG1 + sYWrG1*mW + sYZrG1*mZ
  mYadtG2 := iYrG2 + sYWrG2*mW + sYZrG2*mZ
  
  # Effects
  direct := mYadjG2 - mYadjG1             
  total := mYadtG2 - mYadtG1
  indirect := total - direct
  
  ga10 := iYG2 - iYG1
  ga11 := sYWG2 - sYWG1
  ga12 := sYZG2 - sYZG1
  ga13 := sYMG2 - sYMG1
  
  nde0 := ga10 + ga11*mW + ga12*mZ + ga13*mMadjG1 # natural direct effect 0
  nde1 := ga10 + ga11*mW + ga12*mZ + ga13*mMadjG2 # natural direct effect 1
  adet := ga10 + ga11*mWG2 + ga12*mZG2 + ga13*mMG2 # average direct effect on the treated
  
  '
  
  mm_ML <- generateMeasurementModel(names = "M",
                                    indicators = list("M" = paste0("m", 1:num_m_vars)),
                                    ncells = 2,
                                    model = c("tau-cong"),
                                    data = as.data.frame(sample_compl_disc_cont))
  mm_ML <- gsub("c\\(1,1\\)\\*m1",
                paste0("c\\(", condition$loading, ",", condition$loading, "\\)*m1"),
                mm_ML)
  mm_ML <- gsub("m1 ~ c\\(0,0\\)\\*1",
                paste0("m1 ~ c\\(", condition$mean_Mpre_G1, ",", condition$mean_Mpre_G2, "\\)*1"),
                mm_ML)
  
  mm_ML <- paste(mm_ML, model)
  
  model_ML <- sem(mm_ML, data = sample_compl_disc_cont, 
                  group = "x",
                  estimator = "ML")
  summary(model_ML)

  
  model_MLR <- sem(mm_ML, data = sample_compl_disc_cont, 
                  group = "x",
                  estimator = "MLR")
  summary(model_MLR)
  
  mm_disc <- generateMeasurementModel(names = "M",
                                      indicators = list("M" = paste0("m", 1:num_m_vars)),
                                      ncells = 2,
                                      model = c("tau-cong-categorical"),
                                      data = as.data.frame(sample_compl_disc))
  mm_disc <- gsub("M =~ c\\(1,1\\)\\*m1", 
                  paste0("M =~ c(", condition$loading, ",", condition$loading, ")*m1"), 
                  mm_disc)
  for (i in 1:num_m_vars) {
    mm_disc <- gsub(paste0("m", i, " ~ c\\(0,0\\)\\*1"),
                    paste0("m", i, " ~ c(", condition$mean_Mpre_G1, ",", condition$mean_Mpre_G2, ")*1"),
                    mm_disc)
  }
  mm_disc <- gsub("m1 \\| c\\(0,0\\)\\*t1", 
                  paste0("m1 | c(", condition$thresholds[[1]][2], ",", 
                         condition$thresholds[[1]][2], ")*t1"), 
                  mm_disc)
  mm_disc <- paste(mm_disc, model)
  
  model_DWLS <- sem(mm_disc, data = sample_compl_disc, 
                    group = "x", 
                    estimator = "DWLS")
  summary(model_DWLS)
  
  
  res <- list(ML = model_ML,
              MLR = model_MLR,
              DWLS = model_DWLS,
              cov_matrix_pop = list("Sigma_nT" = Sigma_nT, 
                                    "Sigma_wT" = Sigma_wT),
              cop_data = list("vine_cop_nT" = cop_data_nT[[2]], 
                              "vine_cop_wT" = cop_data_wT[[2]]),
              Method = c("model_ml", "model_mlr", "model_dwls"),
              replication = repNr,
              Seed = .Random.seed)
  return(res)
}
