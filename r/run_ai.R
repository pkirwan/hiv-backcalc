################## Age Independent CD4 Back Calculation ##################

here::i_am("r/run_ai.R")

library(here)
library(rstan)

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(here("r/ai_functions.R"))

# quarterly data
load(here("data/test_data_ai.RData"))

# Incidence models
## 1) Random walk starting in 1978
## 2) Random walk starting in 1995
## 3) Thin plate regression spline with shrinkage (Wood 2003 + Marra and Wood 2011) - 10 knots
## 4) Cubic B-spline with a first order difference penalty (Eilers and Marx 1996)  - 10 knots
## 5) GP with 0 mean and quadratic exponential kernel

# Diagnosis models
## 1) Random walk
## 2) Random walk with added (state-specific) linear term giving a temporal trend
## 3) RW ?

# build_ai() selects the correct model code, generates a jagam design matrix, and builds the model
model <- build_ai(data = test_data_ai, incidence_model = 3, diagnosis_model = 1)

# sample from the backcalculation model
fit <- sampling(
  object = model$stan_model,
  data = model$stan_data,
  iter = 2000,
  chains = 4,
  seed = sample.int(.Machine$integer.max, 1),
  pars = model$pars_save,
  control = list(adapt_delta = 0.95),
  save_warmup = FALSE
)

model_fit_ai <- as.matrix(fit)

# post_process_backcalc() processes the model samples and returns a list of derived outputs
postproc_ai <- post_process_backcalc(model_fit_ai)

save(list = c(model_fit_ai, postproc_ai), file = here("data/postproc_ai.RData"))
