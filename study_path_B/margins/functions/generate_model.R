generate_model <- function(sMZ,
                           sMW,
                           sYW,
                           sYM,
                           sYZ,
                           MpreYpre_corr,
                           loading) {
  
  # Define the SEM model specification
  model <- '
  Mpre =~ start(loading)*m1 + start(loading)*m2 + start(loading)*m3
  m1 ~~ start(latvar)*m1
  m2 ~~ start(latvar)*m2
  m3 ~~ start(latvar)*m3
  
  M ~ sMZ*Mpre 
  Y ~ sYW*Ypre + sYM*M 
  
  # Covariances
  Mpre ~~ MpreYpre_corr*Ypre
  '
  
  # Substitute placeholders with actual parameter values
  model <- gsub("sMZ", as.character(sMZ), model, fixed = TRUE)
  model <- gsub("sMW", as.character(sMW), model, fixed = TRUE)
  model <- gsub("sYW", as.character(sYW), model, fixed = TRUE)
  model <- gsub("sYM", as.character(sYM), model, fixed = TRUE)
  model <- gsub("sYZ", as.character(sYZ), model, fixed = TRUE)
  model <- gsub("MpreYpre_corr", as.character(MpreYpre_corr), model, fixed = TRUE)
  model <- gsub("loading", as.character(loading), model, fixed = TRUE)
  
  latvar <- 1 - loading^2
  model <- gsub("latvar", as.character(latvar), model, fixed = TRUE)
  
  # Fit SEM and capture potential errors
  fit <- try(sem(model = model, data = NULL), silent = TRUE)
  if (class(fit) == "try-error") {
    stop("Model fitting failed.")
  }
  
  # Extract and return the fitted covariance matrix
  Sigma <-  fitted(fit)$cov
  
  return(Sigma)
}
