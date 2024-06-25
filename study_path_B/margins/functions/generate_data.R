generate_data <- function(n, mean_Ypre, mean_Y, mean_Mpre, mean_M,
                          Sigma, margin_dist_M, margin_dist_Ypre, 
                          marginsM_rate, copula, Nmax, num_m_vars) {
  
  marginal_distrs <- lapply(sqrt(diag(Sigma)), function(sd) list(distr = "norm", mean = 0, sd = sd))
  
  marginal_distrs$Y$mean <- mean_Y
  marginal_distrs$Ypre$mean <- mean_Ypre
  marginal_distrs$M$mean <- mean_M
  
  if (margin_dist_M == "exp") {
    for (i in 1:num_m_vars) {
      m_key <- paste0("m", i)
      marginal_distrs[[m_key]] <- list(distr = "exp", rate = marginsM_rate)
    }
  } else {
    for (i in 1:num_m_vars) {
      marginal_distrs[[paste0("m", i)]]$mean <- mean_Mpre
    }
  }
  
  pcs <- list(
    c(replicate(2, bicop_dist(copula), simplify = FALSE), replicate(3, bicop_dist("gaussian"), simplify = FALSE)),
    c(replicate(2, bicop_dist(copula), simplify = FALSE), replicate(2, bicop_dist("gaussian"), simplify = FALSE)),
    c(replicate(2, bicop_dist(copula), simplify = FALSE), list(bicop_dist("gaussian"))),
    replicate(2, bicop_dist(copula), simplify = FALSE),
    list(bicop_dist(copula))
  )
  
  vine_structure <- dvine_structure(1:6)
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