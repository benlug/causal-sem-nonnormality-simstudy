generate_model <- function(sMZ, sMW, sYW, sYM, sYZ, MpreYpre_corr, loading) {
  model <- paste0(
    "M =~ start(", loading, ")*m1 + start(", loading, ")*m2 + start(", loading, ")*m3\n",
    "m1 ~~ start(", 1 - loading^2, ")*m1\n",
    "m2 ~~ start(", 1 - loading^2, ")*m2\n",
    "m3 ~~ start(", 1 - loading^2, ")*m3\n",
    "M ~~ start(", 1 - sMZ^2, ")*M\n",
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