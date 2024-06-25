# Libraries and Setup
library(pacman)
p_load(copula, rvinecopulib, covsim, mvtnorm, lavaan, EffectLiteR, 
       tidyverse, foreach, progressr, parallel, doParallel)

# Set working directory and configure progress handling
setwd(this.path::this.dir())
handlers(global = TRUE)
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
  width    = 80,
  complete = "+"
))

# Parallel processing setup
n_cores <- parallel::detectCores() - 1
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = if (Sys.info()['sysname'] == "Linux") "FORK" else "PSOCK"
)
doParallel::registerDoParallel(cl = my_cluster)

# Helper function
calculate_bins <- function(n_bins) {
  max_mean <- 1 / 0.9
  upper_bound <- max_mean * 3
  bin_edges <- seq(0, upper_bound, length.out = n_bins + 1)
  bin_edges[length(bin_edges)] <- Inf
  bin_edges
}

# Simulation parameters
sim_conditions <- expand_grid(
  n = c(100, 250),
  sim_rep = 200,
  sMZG1 = 0.8, sMZG2 = 0.8,
  sMWG1 = 0, sMWG2 = 0, 
  sYWG1 = 0.75, sYWG2 = 0.75,
  sYMG1 = 0.5, sYMG2 = 0.5,
  sYZG1 = 0, sYZG2 = 0,
  mean_Ypre_G1 = 0, mean_Ypre_G2 = 0,
  mean_Y_G1 = 0, mean_Y_G2 = 1,
  mean_Mpre_G1 = 0, mean_Mpre_G2 = 0,
  mean_M_G1 = 0, mean_M_G2 = 0.3,
  MpreYpre_corr = 0.3,
  loading = c(0.8, 0.9),
  thresholds = list(
    sym_3 = calculate_bins(3),
    sym_4 = calculate_bins(4),
    sym_5 = calculate_bins(5)
  ),
  marginsM_nT = "exp", marginsM_wT = "exp",
  marginsM_theta_nT = 1.2,
  marginsM_theta_wT = c(1.1, 1.3),
  marginsYpre_nT = "norm", marginsYpre_wT = "norm",
  copula_nT = c("gauss", "clayton", "joe"),
  copula_wT = c("gauss", "clayton", "joe")
)

# Simulation function
run_simulation_condition <- function(cond_nr, sim_fnc) {
  condition <- sim_conditions[cond_nr, ]
  
  results <- foreach(rep = seq_len(condition$sim_rep),
                     .packages = names(sessionInfo()$otherPkgs),
                     .errorhandling = "pass") %dopar% {
                       sim_fnc(repNr = rep, condition = condition)
                     }
  
  save(results, file = paste0("simres_", cond_nr, ".rda"))
}

# Load simulation function
source("functions/one_rep.R")

# Run simulation
start_time <- Sys.time()
with_progress({
  p <- progressor(along = seq_len(nrow(sim_conditions)))
  for (cond_nr in seq_len(nrow(sim_conditions))) {
    run_simulation_condition(cond_nr = cond_nr, sim_fnc = one_rep)
    p()
  }
})

# Save simulation info
simtime <- difftime(Sys.time(), start_time, units = "auto")
sesinfo <- sessionInfo()
save(simtime, sesinfo, file = "siminfo.rda")
save(sim_conditions, file = "sim_conditions.rda")

# Stop cluster
parallel::stopCluster(cl = my_cluster)