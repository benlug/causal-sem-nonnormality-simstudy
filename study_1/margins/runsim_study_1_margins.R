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

setwd(this.path::this.dir())

# create and register cluster
n_cores <- parallel::detectCores()-2
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = if_else(Sys.info()['sysname'] == "Linux", "FORK", "PSOCK")  # FORK=linux | PS0CK=windows
)
# print(my_cluster)  # check cluster and number of nodes
doParallel::registerDoParallel(cl = my_cluster)  # register cluster
# foreach::getDoParRegistered() # check if cluster is registered

calculate_bins <- function(n_bins, n_std_devs = 3) {
  bin_edges <- seq(0, n_std_devs, length.out = n_bins + 1)
  bin_edges[1] <- 0
  bin_edges[length(bin_edges)] <- Inf
  return(bin_edges)
}

# simulation parameters ---------------------------------------------------

# n:                  sample size
# sim_rep:            number of simulation replications
# beta_z_nT:          regression coefficient for ordinal covarate in control group
# beta_z_wT:          regression coefficient for ordinal covarate in treatment group
# load_z:             factor loading for z
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
                               beta_z_nT = c(0.3), beta_z_wT = c(0.5, 0.7),
                               load_z = c(0.9, 0.8), 
                               mean_z_nT = c(0), mean_z_wT = c(0),
                               treat_eff = c(1),
                               thresholds = list(sym_3 = calculate_bins(3),
                                                 sym_4 = calculate_bins(4),
                                                 sym_5 = calculate_bins(5)),
                               marginsY_nT = c("norm"), marginsY_wT = c("norm"),
                               marginsZ_nT = c("exp"), marginsZ_wT = c("exp"),
                               marginsZ_theta_nT = c(1),
                               marginsZ_theta_wT = c(0.9, 1, 1.1),
                               copulaZY_nT = c("gauss", "clayton", "joe"),
                               copulaZY_wT = c("gauss", "clayton", "joe"))

# set specific condition (for testing purposes) ---------------------------

repNr <- 1
cond_nr <- 1
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
  p <- progressor(along = 515:(nrow(sim_conditions)))  # seq_len 
  for (cond_nr in 515:(nrow(sim_conditions))) {        # seq_len 
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
