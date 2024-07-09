# libraries + working dir + parallel backend ------------------------------

packages <- c("copula", "rvinecopulib", "covsim", "mvtnorm", 
              "lavaan", "EffectLiteR", 
              "tidyverse", 
              "foreach", "progressr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

handlers(global = TRUE)
handlers(
  handler_progress(
    format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
    width    = 80,
    complete = "+")
)

# set working directory
setwd(this.path::this.dir())

# create and register cluster
n_cores <- parallel::detectCores()-1
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = if_else(Sys.info()['sysname'] == "Linux", "FORK", "PSOCK")  # FORK=linux | PS0CK=windows
)
# print(my_cluster)  # check cluster and number of nodes
doParallel::registerDoParallel(cl = my_cluster)  # register cluster
# foreach::getDoParRegistered() # check if cluster is registered

calculate_bins <- function(n_bins, n_std_devs) {
  bin_edges <- seq(-n_std_devs, n_std_devs, length.out = n_bins + 1)
  bin_edges[1] <- -Inf
  bin_edges[length(bin_edges)] <- Inf
  return(bin_edges)
}
calculate_skewed_bins <- function(n_bins, mean = 0, variance = 1, skew = "left") {
  if (n_bins < 3) {
    stop("Number of bins should be at least 3")
  }
  
  linear_probs <- seq(0, 1, length.out = n_bins + 1)
  
  # Adjust probabilities based on skewness
  if (skew == "left") {
    adjusted_probs <- linear_probs ^ 1.5
  } else {  # For right skew
    adjusted_probs <- (1 - (1 - linear_probs) ^ 1.5)
  }
  
  # Normalize the probabilities to ensure they range from 0 to 1
  adjusted_probs <- adjusted_probs / max(adjusted_probs)
  
  # Adjust the first and last probabilities to avoid -Inf and Inf
  adjusted_probs[1] <- 1e-10
  adjusted_probs[length(adjusted_probs)] <- 1 - 1e-10
  
  std_dev <- sqrt(variance)
  bin_edges <- qnorm(adjusted_probs, mean, std_dev)
  bin_edges[1] <- -Inf
  bin_edges[length(bin_edges)] <- Inf
  
  return(bin_edges)
}

# simulation parameters ---------------------------------------------------

# n:                  sample size
# sim_rep:            number of simulation replications
# beta_z_nT:          regression coefficient for ordinal covarate in control group
# beta_z_wT:          regression coefficient for ordinal covarate in treatment group
# load_z:             factor loading for y 
# mean_z_nT:          mean of covariate Z in control group
# mean_z_wT:          mean of covariate Z in treatment group
# treat_eff:          effect of treatment `x` on dependent variable
# thresholds:         thresholds for discretizing the latent variable
# marginsY_nT:        marginal distribution of response variable in control group
# marginsY_wT:        marginal distribution of response variable in treatment group
# marginsZ_nT:        marginal distribution of ordinal covariate in control group
# marginsZ_wT:        marginal distribution of ordinal covariate in treatment group
# copulaZY_nT:        copula for ordinal covariate and response in control group
# copulaZY_wT:        copula for ordinal covariate and response in treatment group  

sim_conditions <-  expand_grid(n = c(100, 250), sim_rep = 200,
                               sMZG1 = c(0.8), sMZG2 = c(0.8),
                               sMWG1 = c(0), sMWG2 = c(0), 
                               sYWG1 = c(0.75), sYWG2 = c(0.75),
                               sYMG1 = c(0.5), sYMG2 = c(0.5),
                               sYZG1 = c(0), sYZG2 = c(0),
                               mean_Ypre_G1 = 0, mean_Ypre_G2 = 0,
                               mean_Y_G1 = 0, mean_Y_G2 = 1,
                               mean_Mpre_G1 = 0, mean_Mpre_G2 = 0,
                               mean_M_G1 = 0, mean_M_G2 = 1,
                               MpreYpre_corr = c(0.3),
                               loading = c(0.8),  # 0.9
                               thresholds = list(sym_3 = calculate_bins(3, 1.8),
                                                 sym_4 = calculate_bins(4, 1.8),
                                                 sym_5 = calculate_bins(5, 1.8),
                                                 left_4 = calculate_skewed_bins(4, skew = "left"),
                                                 left_5 = calculate_skewed_bins(5, skew = "left"),
                                                 right_4 = calculate_skewed_bins(4, skew = "right"),
                                                 right_5 = calculate_skewed_bins(5, skew = "right")),
                               marginsM_nT = c("norm"), marginsM_wT = c("norm"),
                               marginsYpre_nT = c("norm"), marginsYpre_wT = c("norm"),
                               copula_nT = c("gauss", "clayton", "joe"),
                               copula_wT = c("gauss", "clayton", "joe"))

# set specific condition (for testing purposes) ---------------------------

repNr <- 22
cond_nr <- 21
condition <- sim_conditions[cond_nr, ]

# parallelize simulation replications -------------------------------------

mc_run <- function(cond_nr, sim_fnc) {
  
  condition <- sim_conditions[cond_nr, ]
  
  results <- foreach(rep = 1L:condition$sim_rep,
                     .packages = names(sessionInfo()$otherPkgs),
                     .errorhandling = "pass") %dopar% {
                       sim_fnc(repNr = rep, condition = condition)
                     }
  # save results
  save(results, file = paste0("simres_", cond_nr, ".rda"))
}

# choose simulation function ----------------------------------------------
source("functions/one_rep.R")

# run simulation ----------------------------------------------------------
tmp <- Sys.time()
with_progress({
  p <- progressor(along = seq_len(nrow(sim_conditions)))  # seq_len 
  for (cond_nr in seq_len(nrow(sim_conditions))) {
    sim_dat <- mc_run(cond_nr = cond_nr, sim_fnc = one_rep)
    p()
  }
})
simtime <- Sys.time() - tmp
sesinfo <- sessionInfo()
save(simtime, sesinfo, file = paste0("siminfo.rda"))
save(sim_conditions, file = paste0("sim_conditions.rda"))

# stop cluster ------------------------------------------------------------

parallel::stopCluster(cl = my_cluster)
