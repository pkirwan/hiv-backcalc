################## Age Independent CD4 Back Calculation ##################
################### Process results of the MCMC model ####################

# setup
here::i_am("r/process_ai.R")

# libraries
library(here)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

# source rcpp and r functions
Rcpp::sourceCpp(here("cpp/backcalcfns.cpp"))
source(here("r/main_ai_functions.R"))

# 2022 results
load(here("data/postsim_ai.RData"))

# outputs
hiv_ai <- hivlist(postsim_ai)

save(hiv_ai, file = here("data/postproc_ai.RData"))
