generate_data <- function(n, mean_y, mean_z, Sigma, 
                          margin_dist_y, margin_dist_z, 
                          marginsY_rate, 
                          copula_type_ZY, Nmax, 
                          num_y_vars) {
  
  # Define marginal distributions based on Sigma diagonal
  marginal_distrs <- lapply(X = sqrt(diag(Sigma)), 
                            function(X) list(distr = "norm", mean = 0, sd = X))
  
  for (i in 1:num_y_vars) {
    marginal_distrs[[paste0("y", i)]]$mean <- mean_y
  }
  marginal_distrs$z$mean <- mean_z
  
  if (margin_dist_y == "exp") {
    for (i in 1:num_y_vars) {
      y_key <- paste("y", i, sep = "")
      marginal_distrs[[y_key]]$distr <- margin_dist_y
      marginal_distrs[[y_key]]$rate <- marginsY_rate    
      
      marginal_distrs[[y_key]]$mean <- NULL
      marginal_distrs[[y_key]]$sd <- NULL 
    }
  }
  
  # Pair Copula Specifications for a D-vine structure
  pcs <- list(
    list(bicop_dist(as.character(copula_type_ZY)), bicop_dist(as.character(copula_type_ZY)), bicop_dist("gaussian")),
    list(bicop_dist(as.character(copula_type_ZY)), bicop_dist("gaussian")),
    list(bicop_dist("gaussian"))
  )
  
  # Vine Structure
  # Assuming a D-vine structure for simplicity
  vine_structure <- dvine_structure(1:4) # For 4 variables (y1, y2, y3, z)
  
  # Create Vine Copula Distribution
  vine_cop <- vinecop_dist(pcs, structure = vine_structure)
  
  # Generate the copula calibration with error handling
  calib <- vita(margins = marginal_distrs, 
                sigma.target = Sigma, 
                vc = vine_cop,
                # family_set = copula_type_ZY,
                cores = 1,
                Nmax = Nmax)

  sample <- rvine(n, vine = calib)
  colnames(sample) <- names(marginal_distrs)
  sample <- as_tibble(sample)
  
  return(list(sample, calib))
}



