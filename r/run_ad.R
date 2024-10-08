## ########################## AGE-SPECIFIC MODEL ##########################

here::i_am("r/run_ad.R")

library(here)
library(rstan)

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(here("r/ad_functions.R"))

# quarterly data
load(here("data/model_ad_rita_2011_2023_data6.RData"))

### Model ids
## 1) Ptens - qt, no age diag (standard)
## 2) Ptens - qt, no age diag, but linear trend
## 3) Ptens - qt, age diag 4: intercept - age specific
## 4) Ptens - qt, age diag 5: intercept - age and state specific
## 5) Ptens - qt, age diag: time trend - age specific
## 6) Ptens - qt, age diag 1: time trend - age and state specific
## 7) Ptens - qt, age diag 2, no urep - indept rw for each age class
## 8) Ptens - qt, age diag 3, no urep - splines for dx probs
## 9) Ptens - qt, age diag and rita: time trend - age and state specific, with rita

model <- build_ad(data = stan_data, model_id = 9, data_id = 6)

fit <- sampling(
  object = model$stan_model,
  data = model$stan_data,
  iter = 1000,
  chains = 3,
  pars = model$pars_save,
  control = list(adapt_delta = 0.8)
)

save(model, fit, file = here("data/postsim_ad_rita_2011_2023.RData"))
