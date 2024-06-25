generate_model <- function(beta_z, load_z) {
  
  # Define the SEM model specification
  model <- '
  z =~ start(load_z)*z1 + start(load_z)*z2 + start(load_z)*z3
  y ~ start(beta_z)*z
  y ~~ start(resvar)*y
  z1 ~~ start(latvar)*z1
  z2 ~~ start(latvar)*z2
  z3 ~~ start(latvar)*z3
  '
  
  # Compute latent and residual variances
  latvar <- 1 - load_z^2
  resvar <- 1 - beta_z^2
  
  # Substitute placeholders with actual parameter values
  model <- gsub("beta_z", as.character(beta_z), model, fixed = TRUE)
  model <- gsub("load_z", as.character(load_z), model, fixed = TRUE)
  model <- gsub("resvar", as.character(resvar), model, fixed = TRUE)
  model <- gsub("latvar", as.character(latvar), model, fixed = TRUE)
  
  fit <- sem(model = model, data = NULL)
  # Extract and return the fitted covariance matrix
  Sigma <-  fitted(fit)$cov
  
  return(Sigma)
}