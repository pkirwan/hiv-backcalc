main_ai <- function(stan_data, incidence_model = 1, cf = "", diagnosis_model = 1) {
    # load the required libraries
    library(Rcpp)
    library(RcppArmadillo)
    library(rstan)
    library(posterior)
    library(mgcv)

    # rstan options
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)

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

    # Counterfactual models
    ## 1) CF1
    ## 2) CF2

    ## ## Specify a random seed
    rseed <- sample.int(.Machine$integer.max, 1)

    if (incidence_model == 1) {
        model_code <- "RW1978"
        pars_save <- c("gamma", "var1", "var2", "vardelta", "d", "exp_HIV_dx", "exp_AIDS_dx", "exp_p", "under_rep", "c_78", "c_84", "prev_mat")
    }
    if (incidence_model == 2) {
        model_code <- "RW1995"
        pars_save <- c("gamma", "var1", "vardelta", "d", "log_lik", "exp_HIV_dx", "exp_AIDS_dx", "exp_p", "prev_mat")
    }
    if (incidence_model == 3) {
        b <- "ts"
        k <- 10 ## Setting a default number of knots
        m <- 2 ## Setting the default order of the penalty
        model_code <- "spl"
        pars_save <- c("beta", "gamma", "lambda", "sigma", "vardelta", "d", "log_lik", "exp_HIV_dx", "exp_AIDS_dx", "exp_p", "prev_mat")
    }
    if (incidence_model == 4) {
        b <- "bs"
        k <- 10
        m <- c(2, 1)
        model_code <- "spl"
        pars_save <- c("beta", "gamma", "lambda", "sigma", "vardelta", "d", "log_lik", "exp_HIV_dx", "exp_AIDS_dx", "exp_p", "prev_mat")
    }
    if (incidence_model == 5) {
        stan_data$sigma_sq <- 0.001 ## fixed nugget for gp
        range01 <- function(x) {
            (x - min(x)) / (max(x) - min(x))
        } ## scales inputs of GP to 0-1.
        stan_data$x <- range01(1:stan_data$nquar)
        model_code <- "GP"
        pars_save <- c("gamma", "inv_rho", "eta", "vardelta", "d", "log_lik", "exp_HIV_dx", "exp_AIDS_dx", "exp_p", "prev_mat")
    }
    if (incidence_model == 6) {
        b <- "ts"
        k <- 10 ## Setting a default number of knots
        m <- 2 ## Setting the default order of the penalty
        model_code <- "rita"
        pars_save <- c("beta", "gamma", "lambda", "sigma", "vardelta", "d", "log_lik", "exp_HIV_dx", "exp_AIDS_dx", "exp_p", "prev_mat")
    }

    # if(diagnosis_model == 1) nothing to do

    if (diagnosis_model == 2) {
        model_code <- paste0(model_code, "_linDelta")
        pars_save <- c(pars_save, "alpha")
    } else if (diagnosis_model == 3) {
        model_code <- paste0(model_code, "_incDelta")
    }

    # Spline models require a design matrix (x), generated using the jagam function from the mgcv package
    if (incidence_model %in% c(3, 4, 6)) {
        ## fake data
        jags_data <- list(
            qts = 1:stan_data$nquar,
            D = runif(stan_data$nquar)
        )

        # This creates the spline object and using a bug file with some stan code
        # diagonalize = TRUE makes a reparameterization so that the priors are iid normal
        jagam_out <- jagam(D ~ s(qts, bs = b, k = k, m = m),
            family = gaussian, data = jags_data,
            file = here("stan/AI_spl.bug"), diagonalize = TRUE
        )

        ## Grab the design matrix and augment stan_data
        stan_data$X <- jagam_out$jags.data$X
        stan_data$ninfpars <- ncol(stan_data$X)
    }

    ## Derive the correct STAN model file and build the model
    model <- stan_model(file = here("stan", paste0("AI_", model_code, cf, ".stan")))

    ## Sample from model
    ## Don't worry about warnings: "The following numerical problems occurred the indicated number of times on chain...." and "Rejecting initial value..."
    ## Only worry when warning: "There are a number of divergent transitions" (that should not happen though)

    # equivalent for rstan

    fit <- sampling(model,
        data = stan_data,
        iter = 2000,
        chains = 4,
        seed = rseed,
        pars = pars_save,
        control = list(adapt_delta = 0.95),
        save_warmup = FALSE
    )

    fit_mat <- as.matrix(fit)

    # equivalent cmdstanr code
    # model <- cmdstan_model(stan_file = here("stan", paste0("AI_", model_code, cf, ".stan")))
    #
    # fit <- model$sample(
    #     data = stan_data,
    #     seed = rseed,
    #     chains = 4,
    #     parallel_chains = 4,
    #     adapt_delta = 0.95,
    #     refresh = 500
    # )

    # fit_mat <- as_draws_matrix(
    #     fit$draws(variables = pars_save)
    # )

    return(list(fit_mat, stan_data))
}

hivlist <- function(fit_list, dist_fns = "r/process_ai_functions.R", start_yr = 1995, seed = NA) {
    ## Setup for running in parallel
    # plan(multicore) # may not be able to allocate enough memory so only use 4 workers
    # https://cran.r-project.org/web/packages/future/vignettes/future-4-non-exportable-objects.html
    # plan(multisession) and plan(callr) are no good with RCpp functions

    library(future.apply)
    plan(tweak(multicore, workers = 4))

    fit_mat <- fit_list[[1]]
    stan_data <- fit_list[[2]]

    source(here(dist_fns))

    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    
    inf_dist <- future_sapply(1:nrow(fit_mat), function(i) infectDistByYearDiagnosis(stan_data, fit_mat, i, seed = seed), future.seed = TRUE)
    prev_dist <- future_sapply(1:nrow(fit_mat), function(i) prevByYearInfection(stan_data, fit_mat, i, seed = seed), future.seed = TRUE)

    incidence <- incidence.df(stan_data, prev_dist, start.yr = start_yr, annual = FALSE)
    incidence_annual <- incidence.df(stan_data, prev_dist, start.yr = start_yr, annual = TRUE)
    undiag_prev <- undiag.prev.df(stan_data, prev_dist, inf_dist, start.yr = start_yr)
    diag_prob <- diag.prob.df(stan_data, fit_mat, start.yr = start_yr)
    diagnoses <- diagnoses.df(stan_data, prev_dist, start.yr = start_yr)

    hiv_list <- list(
        incidence = incidence,
        incidence_annual = incidence_annual,
        undiag_prev = undiag_prev,
        diag_prob = diag_prob,
        diagnoses = diagnoses
    )

    return(hiv_list)
}
