generate_data <- function(n, 
                          mean_Ypre, mean_Y, 
                          mean_Mpre, mean_M,
                          Sigma, 
                          margin_dist_M, margin_dist_Ypre, 
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
      bicop_dist(as.character(copula)),
      bicop_dist("gaussian"),
      bicop_dist("gaussian")
    ),
    # Third tree: 3 pair copulas
    list(
      bicop_dist(as.character(copula)),
      bicop_dist(as.character(copula)),
      bicop_dist("gaussian")
    ),
    # Fourth tree: 2 pair copulas
    list(
      bicop_dist(as.character(copula)),
      bicop_dist(as.character(copula))
    ),
    # Fifth tree: 1 pair copula
    list(bicop_dist(as.character(copula)))
  )
  
  # Vine Structure for 6 variables (m1, m2, m3, M, Y, Ypre)
  vine_structure <- dvine_structure(1:6)
  
  # Create Vine Copula Distribution
  vine_cop <- vinecop_dist(pcs, structure = vine_structure)
  
  # Assuming 'Sigma' is your matrix
  # variable_names <- colnames(Sigma)
  # 
  # for (tree in 1:length(pcs)) {
  #   cat("Tree", tree, ":\n")
  #   for (pair in 1:(length(variable_names) - tree)) {
  #     # Map indices to variable names
  #     var1_name <- variable_names[pair]
  #     var2_name <- variable_names[pair + 1]
  # 
  #     # Get copula type
  #     copula_type <- pcs[[tree]][[pair]]
  # 
  #     # Extract the family of the copula
  #     copula_family <- if(is.list(copula_type) && !is.null(copula_type$family)) {
  #       copula_type$family
  #     } else {
  #       "<Unknown Copula Family>"
  #     }
  # 
  #     # Print the pair and their copula family
  #     cat(var1_name, "and", var2_name, ":", copula_family, "\n")
  #   }
  #   cat("\n")
  # }
  
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

