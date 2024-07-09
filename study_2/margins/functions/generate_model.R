generate_model <- function(beta_z, load_z) {

  if (!is.numeric(beta_z) || !is.numeric(load_z)) {
    stop("Input parameters must be numeric!")
  }
  
  model <- '
  y =~ start(load_z)*y1 + start(load_z)*y2 + start(load_z)*y3
  y ~ start(beta_z)*z
  y ~~ start(resvar)*y
  y1 ~~ start(latvar)*y1
  y2 ~~ start(latvar)*y2
  y3 ~~ start(latvar)*y3
  '
  
  latvar <- 1 - load_z^2
  zvar <- 1 - load_z^2
  resvar <- 1 - beta_z^2
  
  model <- gsub("beta_z", as.character(beta_z), model, fixed = TRUE)
  model <- gsub("load_z", as.character(load_z), model, fixed = TRUE)
  model <- gsub("resvar", as.character(resvar), model, fixed = TRUE)
  model <- gsub("latvar", as.character(latvar), model, fixed = TRUE)
  model <- gsub("zvar", as.character(zvar), model, fixed = TRUE)
  
  fit <- sem(model = model, data = NULL)
  Sigma <-  fitted(fit)$cov
  
  return(Sigma)
}