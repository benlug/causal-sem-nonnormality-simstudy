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
calculate_bins <- function(n_bins, n_std_devs = 3) {
  bin_edges <- c(0, seq(0, n_std_devs, length.out = n_bins)[-1])
  bin_edges[length(bin_edges)] <- Inf
  bin_edges
}

# Simulation parameters
sim_conditions <- expand_grid(
  n = c(100, 250),
  sim_rep = 200,
  beta_z_nT = 0.3,
  beta_z_wT = c(0.5, 0.7),
  load_z = c(0.8, 0.9),
  mean_z_nT = 0,
  mean_z_wT = 0,
  treat_eff = 0.3,
  thresholds = list(
    sym_3 = calculate_bins(3),
    sym_4 = calculate_bins(4),
    sym_5 = calculate_bins(5)
  ),
  marginsY_nT = "exp",
  marginsY_wT = "exp",
  marginsZ_nT = "norm",
  marginsZ_wT = "norm",
  marginsY_theta_nT = 1.2,
  marginsY_theta_wT = c(1.1, 1.3),
  copulaZY_nT = c("gauss", "clayton", "joe"),
  copulaZY_wT = c("gauss", "clayton", "joe")
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