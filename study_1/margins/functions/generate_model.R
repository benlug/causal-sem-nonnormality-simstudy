generate_model <- function(beta_z, load_z) {
  if (!is.numeric(beta_z) || !is.numeric(load_z)) {
    stop("Input parameters must be numeric!")
  }
  
  model <- paste0(
    "z =~ start(", load_z, ")*z1 + start(", load_z, ")*z2 + start(", load_z, ")*z3\n",
    "y ~ start(", beta_z, ")*z\n",
    "y ~~ start(", 1 - beta_z^2, ")*y\n",
    "z1 ~~ start(", 1 - load_z^2, ")*z1\n",
    "z2 ~~ start(", 1 - load_z^2, ")*z2\n",
    "z3 ~~ start(", 1 - load_z^2, ")*z3\n"
  )
  
  fit <- sem(model = model, data = NULL)
  Sigma <- fitted(fit)$cov
  
  return(Sigma)
}