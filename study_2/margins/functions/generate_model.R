generate_model <- function(beta_z, load_z) {
  
  latvar <- 1 - load_z^2
  zvar <- 1 - load_z^2
  resvar <- 1 - beta_z^2
  
  model <- paste0(
    "y =~ start(", load_z, ")*y1 + start(", load_z, ")*y2 + start(", load_z, ")*y3\n",
    "y ~ start(", beta_z, ")*z\n",
    "y ~~ start(", resvar, ")*y\n",
    "y1 ~~ start(", latvar, ")*y1\n",
    "y2 ~~ start(", latvar, ")*y2\n",
    "y3 ~~ start(", latvar, ")*y3\n",
    "z ~~ start(", zvar, ")*z\n"
  )
  
  fit <- tryCatch({
    sem(model = model, data = NULL)
  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })
  
  Sigma <- fitted(fit)$cov
  
  return(Sigma)
}