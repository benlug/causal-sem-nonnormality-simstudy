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
my_cluster <- makeCluster(
  n_cores, 
  type = if (Sys.info()['sysname'] == "Linux") "FORK" else "PSOCK"
)
registerDoParallel(cl = my_cluster)

# Helper functions
calculate_bins <- function(n_bins, n_std_devs) {
  bin_edges <- seq(-n_std_devs, n_std_devs, length.out = n_bins + 1)
  bin_edges[c(1, length(bin_edges))] <- c(-Inf, Inf)
  bin_edges
}

calculate_skewed_bins <- function(n_bins, skew = "left") {
  linear_probs <- seq(0, 1, length.out = n_bins + 1)
  adjusted_probs <- linear_probs ^ 1.5 / max(linear_probs ^ 1.5)
  bin_edges <- qnorm(adjusted_probs)
  
  if (skew != "left") {
    bin_edges[2:(length(bin_edges)-1)] <- rev(bin_edges[2:(length(bin_edges)-1)] * -1)
  }
  
  bin_edges
}

# Simulation parameters
sim_conditions <- expand_grid(
  n = c(100, 250),
  sim_rep = 200,
  beta_z_nT = 0.3,
  beta_z_wT = c(0.5, 0.7),
  load_z = c(0.9, 0.8),
  mean_z_nT = 0,
  mean_z_wT = 0,
  treat_eff = 1,
  thresholds = list(
    sym_3 = calculate_bins(3, 3),
    sym_4 = calculate_bins(4, 3),
    sym_5 = calculate_bins(5, 3),
    left_4 = calculate_skewed_bins(4, skew = "left"),
    left_5 = calculate_skewed_bins(5, skew = "left"),
    right_4 = calculate_skewed_bins(4, skew = "right"),
    right_5 = calculate_skewed_bins(5, skew = "right")
  ),
  marginsY_nT = "norm",
  marginsY_wT = "norm",
  marginsZ_nT = "norm",
  marginsZ_wT = "norm",
  copulaZY_nT = c("gauss", "clayton", "joe"),
  copulaZY_wT = c("gauss", "clayton", "joe")
)

# Simulation function
mc_run <- function(cond_nr, sim_fnc) {
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
  p <- progressor(steps = nrow(sim_conditions))
  for (cond_nr in seq_len(nrow(sim_conditions))) {
    mc_run(cond_nr = cond_nr, sim_fnc = one_rep)
    p(sprintf("Completed condition %d of %d", cond_nr, nrow(sim_conditions)))
  }
})

# Save simulation info
simtime <- difftime(Sys.time(), start_time, units = "auto")
sesinfo <- sessionInfo()
save(simtime, sesinfo, file = "siminfo.rda")
save(sim_conditions, file = "sim_conditions.rda")

# Stop cluster
stopCluster(cl = my_cluster)
