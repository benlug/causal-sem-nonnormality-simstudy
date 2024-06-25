generate_data <- function(n, mean_y, mean_z, Sigma, 
                          margin_dist_y, margin_dist_z, 
                          marginsY_rate, 
                          copula_type_ZY, Nmax, 
                          num_y_vars) {
  
  marginal_distrs <- lapply(sqrt(diag(Sigma)), function(sd) list(distr = "norm", mean = 0, sd = sd))
  
  for (i in 1:num_y_vars) {
    y_key <- paste0("y", i)
    if (margin_dist_y == "exp") {
      marginal_distrs[[y_key]] <- list(distr = "exp", rate = marginsY_rate)
    } else {
      marginal_distrs[[y_key]]$mean <- mean_y
    }
  }
  marginal_distrs$z$mean <- mean_z
  
  pcs <- list(
    replicate(2, bicop_dist(copula_type_ZY), simplify = FALSE),
    bicop_dist("gaussian"),
    list(bicop_dist(copula_type_ZY), bicop_dist("gaussian")),
    list(bicop_dist("gaussian"))
  )
  
  vine_structure <- dvine_structure(1:4)
  vine_cop <- vinecop_dist(pcs, structure = vine_structure)
  
  calib <- vita(margins = marginal_distrs, 
                sigma.target = Sigma, 
                vc = vine_cop,
                cores = 1,
                Nmax = Nmax)
  
  sample <- rvine(n, vine = calib) %>%
    as_tibble() %>%
    setNames(names(marginal_distrs))
  
  list(sample, calib)
}