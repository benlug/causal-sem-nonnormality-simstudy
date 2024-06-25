generate_data <- function(n, mean_y, mean_z, Sigma, 
                          margin_dist_y, margin_dist_z, 
                          copula_type_ZY, Nmax, 
                          num_z_vars) {
  
  # Define marginal distributions based on Sigma diagonal
  marginal_distrs <- lapply(X = sqrt(diag(Sigma)), 
                            function(X) list(distr = "norm", mean = 0, sd = X))
  
  # Adjust means for y variables and z variable
  for (i in 1:num_z_vars) {
    marginal_distrs[[paste0("z", i)]]$mean <- mean_z
  }
  marginal_distrs$y$mean <- mean_y
  
  # Total number of variables
  total_vars <- num_z_vars + 1
  
  # Pair Copula Specifications for a D-vine structure
  pcs <- vector("list", total_vars - 1)
  for (i in 1:(total_vars - 1)) {
    pcs[[i]] <- vector("list", total_vars - i)
    for (j in 1:(total_vars - i)) {
      if (i == 1 && j <= num_z_vars) {
        pcs[[i]][[j]] <- bicop_dist(as.character(copula_type_ZY))
      } else {
        pcs[[i]][[j]] <- bicop_dist("gaussian")
      }
    }
  }
  
  # Vine Structure
  vine_structure <- dvine_structure(total_vars)
  
  # Create Vine Copula Distribution
  vine_cop <- vinecop_dist(pcs, structure = vine_structure)
  
  # Generate the copula calibration with error handling
  calib <- vita(margins = marginal_distrs, 
                sigma.target = Sigma, 
                vc = vine_cop,
                cores = 1,
                Nmax = Nmax)
  sample <- rvine(n, vine = calib)
  colnames(sample) <- names(marginal_distrs)
  sample <- as_tibble(sample)
  
  return(list(sample, calib))
}