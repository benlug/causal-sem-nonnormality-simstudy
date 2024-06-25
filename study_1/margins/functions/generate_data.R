generate_data <- function(n, mean_y, mean_z, Sigma, 
                          margin_dist_y, margin_dist_z, 
                          marginsZ_rate, 
                          copula_type_ZY, Nmax, 
                          num_z_vars) {
  
  marginal_distrs <- lapply(sqrt(diag(Sigma)), function(sd) list(distr = "norm", mean = 0, sd = sd))
  
  for (i in 1:num_z_vars) {
    z_key <- paste0("z", i)
    if (margin_dist_z == "exp") {
      marginal_distrs[[z_key]] <- list(distr = "exp", rate = marginsZ_rate)
    } else {
      marginal_distrs[[z_key]]$mean <- mean_z
    }
  }
  marginal_distrs$y$mean <- mean_y
  
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