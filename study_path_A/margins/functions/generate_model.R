generate_model <- function(sMZ, sMW, sYW, sYM, sYZ, MpreYpre_corr, loading) {
  latvar <- 1 - loading^2
  resvar <- 1 - sMZ^2
  
  model <- paste0(
    "M =~ start(", loading, ")*m1 + start(", loading, ")*m2 + start(", loading, ")*m3\n",
    "m1 ~~ start(", latvar, ")*m1\n",
    "m2 ~~ start(", latvar, ")*m2\n",
    "m3 ~~ start(", latvar, ")*m3\n",
    "M ~~ start(", resvar, ")*M\n",
    "M ~ ", sMZ, "*Mpre\n",
    "Y ~ ", sYW, "*Ypre + ", sYM, "*M\n",
    "Mpre ~~ ", MpreYpre_corr, "*Ypre\n"
  )
  
  fit <- tryCatch({
    sem(model = model, data = NULL)
  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })
  
  fitted(fit)$cov
}