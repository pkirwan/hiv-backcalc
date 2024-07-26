################## Age Independent CD4 Back Calculation ##################

here::i_am("r/run_ai.R")

library(here)
library(future.callr)
plan(callr)

source(here("r/run_ai_functions.R"))
load(here("data/test_data_ai.RData"))

# Incidence models
## 1) Random walk starting in 1978
## 2) Random walk starting in 1995
## 3) Thin plate regression spline with shrinkage (Wood 2003 + Marra and Wood 2011) - 10 knots
## 4) Cubic B-spline with a first order difference penalty (Eilers and Marx 1996)  - 10 knots
## 5) GP with 0 mean and quadratic exponential kernel

postsim_ai <- main_ai(stan_data, incidence_model = 3)

save(postsim_ai, file = here("data/postsim_ai.RData"))
