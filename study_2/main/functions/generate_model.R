generate_model <- function(beta_z, load_z) {

  if (!is.numeric(beta_z) || !is.numeric(load_z)) {
    stop("Input parameters must be numeric!")
  }
  
  # Define the SEM model specification
  model <- '
  y =~ start(load_z)*y1 + start(load_z)*y2 + start(load_z)*y3
  y ~ start(beta_z)*z
  y ~~ start(resvar)*y
  y1 ~~ start(latvar)*y1
  y2 ~~ start(latvar)*y2
  y3 ~~ start(latvar)*y3
  '
  
  # Compute latent and residual variances
  latvar <- 1 - load_z^2
  resvar <- 1 - beta_z^2
  
  # Substitute placeholders with actual parameter values
  model <- gsub("beta_z", as.character(beta_z), model, fixed = TRUE)
  model <- gsub("load_z", as.character(load_z), model, fixed = TRUE)
  model <- gsub("resvar", as.character(resvar), model, fixed = TRUE)
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
