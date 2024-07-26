## ########################## RUN THE AGE-SPECIFIC MODEL ##########################

# the here() package is used to provide relative file locations

here::i_am("r/run_ad.R")

library(mgcv)
library(rstan)
library(here)

# rstan options
rstan_options(auto_write = FALSE)
options(mc.cores = 3)

# detect if the program is being run on the cluster
cluster <- length(grep("cpu", Sys.info()["nodename"], fixed = TRUE)) > 0
## ### Specify some hard-wired defaults
model.ind <- 4
## Specify the type of exclusions made from the dataset
data.id <- 6
## Has an appropriate spline already been coded up?
spl.file <- NULL
# spl.file.new <- here("bug/ptens_10_8.bug")
spl.file.new <- here("bug/tps_80.bug")
## Set the number of iterations
n.iter <- 1000

cat("Model ID = ", model.ind, "\n")

### Model indices
## 1) Ptens - qt, no age diag (standard)
## 2) Ptens - qt, no age diag, but linear trend
## 3) Ptens - qt, age diag 4: intercept - age specific
## 4) Ptens - qt, age diag 5: intercept - age and time specific
## 5) Ptens - qt, age diag: time trend - age specific
## 6) Ptens - qt, age diag 1: time trend - age and state specific
## 7) Ptens - qt, age diag 2, no urep - indept rw for each age class
## 8) Ptens - qt, age diag 3, no urep - splines for dx probs

## ############################################################ INITIAL STUFF ###################################################

## ### Loading Data
## ### loading initial prevalence
load(here("data/init_ad.RData")) ### Age prog, age indpt diag, HIV (same as 3 states)

### Quarterly data
load(here("data/test_data_ad.RData")) ## Will load a data object, `data.95'.

if (!model.ind %in% c(1:10)) stop("model.ind incorrectly specified")

if (model.ind == 1) model.txt <- "ptens_quar_nourep.stan"
if (model.ind == 2) model.txt <- "ptens_quar_nourep_aidiag1.stan"
if (model.ind == 3) model.txt <- "ptens_quar_agediag4_nourep.stan"
# if(model.ind==4) model.txt <- "ptens_quar_agediag5_nourep.stan"
if (model.ind == 4) model.txt <- "tps_quar_agediag5_nourep.stan"
if (model.ind == 5) model.txt <- "ptens_quar_agediag_nourep.stan"
if (model.ind == 6) model.txt <- "ptens_quar_agediag1_nourep.stan"
if (model.ind == 7) model.txt <- "ptens_quar_agediag2_nourep.stan"
if (model.ind == 8) model.txt <- "ptens_quar_agediag3_nourep.stan"

################################################################## USEFUL FUNCTIONS ###############################################################

## ### mgcv splines specifications
## dummy data
tmp <- runif(data.95$num.quarters * data.95$num.ages)

yrs <- 1:data.95$num.quarters ## yrs is actually quarters, but convenient to call it like his
ages <- 1:data.95$num.ages

## Creating a dataset for jagam
jags.data <- list(
  age = rep(ages, each = data.95$num.quarters),
  yrs = rep(yrs, times = data.95$num.ages),
  D = tmp
)

m.list <- list(c(2, 1), c(2, 1))
if (is.null(spl.file)) {
  jagam.out <- jagam(D ~ te(yrs, age, bs = c("ps", "ps"), k = c(10, 8), m = m.list),
    family = gaussian, data = jags.data,
    file = spl.file.new, diagonalize = TRUE
  )
}

x <- jagam.out$jags.data$X
s1 <- jagam.out$jags.data$S1
rm(jagam.out, jags.data, tmp, yrs, ages)

## Final dataset
stan.data <- list(
  "HIV" = data.95$HIV.diagnoses, "AIDS" = data.95$AIDS.diagnoses, "CD4" = data.95$CD4.cell.proportions,
  "nage" = data.95$num.ages, "q" = q.list.4states, "X" = x, "init_prev" = prev.mat, "ninfpars" = ncol(x),
  "S1" = s1, "nmissing" = 0, "nquar" = data.95$num.quarters, "sigma_inform" = 0
)

if (model.ind %in% c(2:6)) p.save <- c("beta", "lambda", "vardelta", "d", "alpha") ## saving age-specific intercept dx parameters
if (model.ind %in% c(8)) p.save <- c("beta", "lambda", "lambda_d1", "lambda_d2", "lambda_d3", "lambda_d4", "d", "delta1", "delta2", "delta3", "delta4")

model <- stan_model(file = here(paste0("stan/", model.txt)))

fit <- sampling(model,
  data = stan.data,
  iter = n.iter,
  chains = 3,
  pars = p.save,
  control = list(adapt_delta = 0.8)
)

save(fit, stan.data, model.ind, data.id, file = here("data/postsim_ad.RData"))
