generate_data <- function(n, 
                          mean_Ypre, mean_Y, 
                          mean_Mpre, mean_M,
                          Sigma, 
                          margin_dist_M, margin_dist_Ypre, 
                          marginsM_rate,
                          copula, Nmax, 
                          num_m_vars) {
  
  marginal_distrs <- lapply(X = sqrt(diag(Sigma)), 
                            function(X) list(distr = "norm", mean = 0, sd = X))
  
  marginal_distrs$Y$mean <- mean_Y
  marginal_distrs$Ypre$mean <- mean_Ypre
  marginal_distrs$Mpre$mean <- mean_Mpre
  
  for (i in 1:num_m_vars) {
    marginal_distrs[[paste0("m", i)]]$mean <- mean_M
  }
  
  if (margin_dist_M == "exp") {
    for (i in 1:num_m_vars) {
      m_key <- paste("m", i, sep = "")
      marginal_distrs[[m_key]]$distr <- margin_dist_M
      marginal_distrs[[m_key]]$rate <- marginsM_rate    
      
      marginal_distrs[[m_key]]$mean <- NULL
      marginal_distrs[[m_key]]$sd <- NULL 
    }
  }
  
  pcs <- list(
    # First tree: 5 pair copulas (for 6 variables)
    list(
      bicop_dist(as.character(copula)), 
      bicop_dist(as.character(copula)),
      bicop_dist("gaussian"),
      bicop_dist("gaussian"),
      bicop_dist("gaussian")
    ),
    # Second tree: 4 pair copulas
    list(
      bicop_dist(as.character(copula)), 
      bicop_dist("gaussian"),
      bicop_dist("gaussian"),
      bicop_dist("gaussian")
    ),
    # Third tree: 3 pair copulas
    list(
      bicop_dist("gaussian"),
      bicop_dist("gaussian"),
      bicop_dist("gaussian")
    ),
    # Fourth tree: 2 pair copulas
    list(
      bicop_dist("gaussian"),
      bicop_dist("gaussian")
    ),
    # Fifth tree: 1 pair copula
    list(bicop_dist("gaussian"))
  )
  
  # Vine Structure for 6 variables (m1, m2, m3, M, Y, Ypre)
  vine_structure <- dvine_structure(1:6)
  
  # Create Vine Copula Distribution
  vine_cop <- vinecop_dist(pcs, structure = vine_structure)
  
  calib <- vita(margins = marginal_distrs, 
                sigma.target = Sigma, 
                # family_set = copula,
                vc = vine_cop,
                cores = 1,
                Nmax = Nmax)
  
  sample <- rvine(n, vine = calib)
  colnames(sample) <- names(marginal_distrs)
  sample <- as_tibble(sample)
  
  return(list(sample, calib))
  
}

