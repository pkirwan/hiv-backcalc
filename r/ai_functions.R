build_ai <- function(
    data,
    incidence_model = 1,
    diagnosis_model = 1) {
  # libraries
  library(mgcv)
  library(rstan)

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
    data$sigma_sq <- 0.001 ## fixed nugget for gp
    range01 <- function(x) {
      (x - min(x)) / (max(x) - min(x))
    } ## scales inputs of GP to 0-1.
    data$x <- range01(1:data$nquar)
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
      qts = 1:data$nquar,
      D = runif(data$nquar)
    )

    # This creates the spline object and using a bug file with some stan code
    # diagonalize = TRUE makes a reparameterization so that the priors are iid normal
    jagam_out <- jagam(D ~ s(qts, bs = b, k = k, m = m),
      family = gaussian, data = jags_data,
      file = here("bug/AI_spl.bug"), diagonalize = TRUE
    )

    ## Grab the design matrix and augment data
    data$X <- jagam_out$jags.data$X
    data$ninfpars <- ncol(data$X)
  }


  ## build the model
  stan_model <- stan_model(file = here("stan", paste0("AI_", model_code, ".stan")))

  return(list(
    "stan_data" = data,
    "stan_model" = stan_model,
    "pars_save" = pars_save
  ))
}

post_process_backcalc <- function(model_fit, data, start_yr = 1995, seed = NA) {
  # libraries
  library(tidyverse)
  library(Rcpp)
  library(RcppArmadillo)

  # Setup for running in parallel using future_sapply
  # use plan(multicore) as plan(multisession) and plan(callr) do not work well with Rcpp
  library(future)
  library(future.apply)
  plan(tweak(multicore, workers = 4))

  # source rcpp and r functions
  Rcpp::sourceCpp(here("cpp/backcalcfns.cpp"))

  if (is.na(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }

  inf_dist <- future_sapply(1:nrow(model_fit), function(i) get_incidence(data, model_fit, i, seed = seed), future.seed = TRUE)
  prev_dist <- future_sapply(1:nrow(model_fit), function(i) get_prevalence(data, model_fit, i, seed = seed), future.seed = TRUE)

  incidence <- incidence.df(data, prev_dist, start.yr = start_yr, annual = FALSE)
  incidence_annual <- incidence.df(data, prev_dist, start.yr = start_yr, annual = TRUE)
  undiag_prev <- undiag.prev.df(data, prev_dist, inf_dist, start.yr = start_yr)
  diag_prob <- diag.prob.df(data, model_fit, start.yr = start_yr)
  diagnoses <- diagnoses.df(data, prev_dist, start.yr = start_yr)

  out <- list(
    incidence = incidence,
    incidence_annual = incidence_annual,
    undiag_prev = undiag_prev,
    diag_prob = diag_prob,
    diagnoses = diagnoses
  )

  return(out)
}

############################# INFECTION AND DIAGNOSED PREVALENCE DISTRIBUTION #############################

### this function generates the HIV.Dx.Infect.Dist and AIDS.Dx.Infect.Dist, by year of diagnosis
get_incidence <- function(stan_data, fit_mat, i, seed = 1707) {
  output <- list()

  n.diag.st <- ncol(stan_data$CD4)
  n.CD4 <- nrow(stan_data$CD4)
  nquar <- stan_data$nquar
  n.iters <- nrow(fit_mat)
  N.sample <- rowSums(stan_data$CD4)
  # diags0 <- c(stan_data$CD4[1, ], stan_data$AIDS[1])

  if (!exists("init_prev", where = stan_data)) {
    init.prev <- rep(0, 4) ### This allows modification for initial prevalence
  } else {
    init.prev <- stan_data$init_prev
  }

  cat("Iteration Number = ", i, "\n")

  d1.ind <- grep("d\\[1,", colnames(fit_mat))
  d2.ind <- grep("d\\[2,", colnames(fit_mat))
  d3.ind <- grep("d\\[3,", colnames(fit_mat))
  d4.ind <- grep("d\\[4,", colnames(fit_mat))

  # these are the diagnosis probabilities for this iteration of the model, for each CD4 state
  d1 <- fit_mat[i, d1.ind]
  d2 <- fit_mat[i, d2.ind]
  d3 <- fit_mat[i, d3.ind]
  d4 <- fit_mat[i, d4.ind]

  # this is the number of new infections for this iteration of the model
  infs.ind <- grep("gamma", colnames(fit_mat))
  temp.h <- exp(fit_mat[i, infs.ind])

  diagnoses <- forward_calculate(init.prev, temp.h, d1 = d1, d2 = d2, d3 = d3, d4 = d4, q = stan_data$q)

  if ("under_rep" %in% colnames(fit_mat)) {
    under.report <- 1 - fit_mat[i, "under_rep"]
    under.report <- c(rep(0, 88), rep(under.report, nquar - 88))
  } else {
    under.report <- rep(0, nquar)
  }

  set.seed(seed)
  unreported.AIDS <- rpois(nquar, under.report * diagnoses$AIDS[1, ])
  for (j in 2:nquar)
  {
    set.seed(seed)
    unreported.AIDS <- rbind(unreported.AIDS, rpois(nquar, under.report * diagnoses$AIDS[j, ]))
  }

  Dx.no.CD4 <- stan_data$HIV - c(rep(0, nquar - n.CD4), N.sample)

  temp.HIV.dx <- temp_diags_fun(diagnoses)
  d.to.i.HIV <- diag_to_infec_wrap(temp.HIV.dx)
  d.to.i.AIDS <- normalize_mat(t(diagnoses$AIDS))

  set.seed(seed)
  temp <- sapply(2:nquar, function(x) rmultinom(1, Dx.no.CD4[x], d.to.i.HIV[, , x]))
  temp <- array(c(matrix(0, n.diag.st, nquar), temp), dim = c(n.diag.st, nquar, nquar))

  output$HIV.Dx.Infect.Dist <- temp

  ### Simulating AIDS diagnoses
  set.seed(seed)
  tempvec <- sapply(5:nquar, function(x) rmultinom(1, stan_data$AIDS[x], d.to.i.AIDS[x, ]))

  output$AIDS.Dx.Infect.Dist <- cbind(0, 0, 0, 0, tempvec) + unreported.AIDS

  for (state in 1:n.diag.st)
  {
    d.to.i.temp <- diag_to_infec_wrap1(diagnoses$prevalence, state)

    HIV.diags <- c(rep(0, nquar - nrow(stan_data$CD4)), stan_data$CD4[, state])

    set.seed(seed)
    tempvec1 <- sapply((state + 1):nquar, function(x) rmultinom(1, HIV.diags[x], d.to.i.temp[x - state, ]))
    output$HIV.Dx.Infect.Dist[state, , ] <- output$HIV.Dx.Infect.Dist[state, , ] + cbind(matrix(0, nquar, nquar - ncol(tempvec1)), tempvec1)
  }

  return(output)
}

# this function generates the unreported.AIDS, undiagnosed, HIV.Dx.no.CD4, HIV.Dx.CD4 and reported.AIDS
# prevalence distributions, by year of infection and as estimated in the final quarter of the model

get_prevalence <- function(stan_data, fit_mat, i, seed = 1707) {
  output <- list()

  n.diag.st <- ncol(stan_data$CD4) ## number of diagnoses states
  n.CD4 <- nrow(stan_data$CD4)
  nquar <- stan_data$nquar
  n.iters <- nrow(fit_mat)
  N.sample <- rowSums(stan_data$CD4) ## number of diagnoses linked to a CD4

  if (!exists("init_prev", where = stan_data)) {
    init.prev <- rep(0, 4) ### This allows modification for initial prevalence
  } else {
    init.prev <- stan_data$init_prev
  }

  cat("Iteration Number = ", i, "\n")

  d1.ind <- grep("d\\[1,", colnames(fit_mat))
  d2.ind <- grep("d\\[2,", colnames(fit_mat))
  d3.ind <- grep("d\\[3,", colnames(fit_mat))
  d4.ind <- grep("d\\[4,", colnames(fit_mat))

  d1 <- fit_mat[i, d1.ind]
  d2 <- fit_mat[i, d2.ind]
  d3 <- fit_mat[i, d3.ind]
  d4 <- fit_mat[i, d4.ind]

  infs.ind <- grep("gamma", colnames(fit_mat))
  temp.h <- exp(fit_mat[i, infs.ind]) ## picking the infection parameters from chain and iteration of interest

  diagnoses <- forward_calculate(init.prev, temp.h, d1 = d1, d2 = d2, d3 = d3, d4 = d4, q = stan_data$q)

  output$AIDS.diagnoses <- matrix(rpois(length(diagnoses$AIDS.diagnoses), diagnoses$AIDS.diagnoses), nquar)
  output$HIV.diagnoses <- matrix(rpois(length(diagnoses$HIV.diagnoses), diagnoses$HIV.diagnoses), nquar)
  # output$expected.p <- matrix(rpois(length(diagnoses$expected.p),diagnoses$expected.p),4,nquar,nquar)

  if ("under_rep" %in% colnames(fit_mat)) {
    under.report <- 1 - fit_mat[i, "under_rep"]
    under.report <- c(rep(0, 88), rep(under.report, nquar - 88))
  } else {
    under.report <- rep(0, nquar)
  }

  set.seed(seed)
  unreported.AIDS <- rpois(nquar, under.report * diagnoses$AIDS[1, ])
  for (j in 2:nquar)
  {
    set.seed(seed)
    unreported.AIDS <- rbind(unreported.AIDS, rpois(nquar, under.report * diagnoses$AIDS[j, ]))
  }

  output$unreported.AIDS <- rowSums(unreported.AIDS)

  # it appears that the undiagnosed component will only hold data for the most recent year
  set.seed(seed)
  output$undiagnosed <- matrix(rpois(n = prod(dim(diagnoses$prevalence[nquar, 1:n.diag.st, ])), lambda = diagnoses$prevalence[nquar, 1:n.diag.st, ]), nr = n.diag.st)

  Dx.no.CD4 <- stan_data$HIV - c(rep(0, nquar - n.CD4), N.sample) ## Number of CD4 count for every quarter

  temp.HIV.dx <- temp_diags_fun(diagnoses)
  d.to.i.HIV <- diag_to_infec_wrap(temp.HIV.dx)
  d.to.i.AIDS <- normalize_mat(t(diagnoses$AIDS))

  set.seed(seed)
  temp <- sapply(2:nquar, function(x) rmultinom(1, Dx.no.CD4[x], d.to.i.HIV[, , x])) ## do not do 1 otherwise troubles (cos not possible to be dx at time 1)
  temp <- array(c(matrix(0, n.diag.st, nquar), temp), dim = c(n.diag.st, nquar, nquar))

  output$HIV.Dx.no.CD4 <- apply(temp, 1:2, sum)
  output$HIV.Dx.CD4 <- matrix(NA, n.diag.st, nquar) ## Initialization needed later

  set.seed(seed)
  tempvec <- sapply(5:nquar, function(x) rmultinom(1, stan_data$AIDS[x], d.to.i.AIDS[x, ]))

  output$reported.AIDS <- rowSums(tempvec)
  for (state in 1:n.diag.st)
  {
    d.to.i.temp <- diag_to_infec_wrap1(diagnoses$prevalence, state)
    HIV.diags <- c(rep(0, nquar - nrow(stan_data$CD4)), stan_data$CD4[, state])
    set.seed(seed)
    tempvec1 <- sapply((state + 1):nquar, function(x) rmultinom(1, HIV.diags[x], d.to.i.temp[x - state, ]))
    output$HIV.Dx.CD4[state, ] <- rowSums(tempvec1)
  }
  return(output)
}


############################# FUNCTIONS FOR GENERATING INCIDENCE, UNDX PREVALENCE AND DIAGNOSIS PROB DATAFRAMES #############################

### Incidence

incidence.df <- function(stan_data, PrevbyYearInf, start.yr = 1978, annual = FALSE) {
  nyr <- stan_data$nquar / 4

  ### Total infections
  ## - Those undx at end of surveillance period (in any states)
  ## - Those Dx with HIV and no CD4 (PrevbyYearInf["HIV.Dx.no.CD4",]) - actually the CD4 is "imputed"
  ## - Those Dx with HIV and CD4 (PrevbyYearInf["HIV.Dx.CD4",])
  ## - Those Dx with AIDS and reported (PrevbyYearInf["unreported.AIDS",])
  ## - Those Dx with AIDS and not reported (PrevbyYearInf["reported.AIDS",])

  infections.total <- sapply(PrevbyYearInf["undiagnosed", ], function(x) colSums(x)) +
    sapply(PrevbyYearInf["HIV.Dx.no.CD4", ], function(x) apply(x, 2, sum)) +
    sapply(PrevbyYearInf["HIV.Dx.CD4", ], function(x) apply(x, 2, sum)) +
    simplify2array(PrevbyYearInf["unreported.AIDS", ]) +
    simplify2array(PrevbyYearInf["reported.AIDS", ])

  if (annual == FALSE) {
    q.infections.total <- apply(infections.total, 1, quantile, probs = c(0.025, 0.5, 0.975))
    df <- as.data.frame(t(q.infections.total))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> pivot_longer(-year, names_to = "quartile")
  } else {
    infections.years <- apply(infections.total, 2, function(x) tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum))
    q.infections.years <- apply(infections.years, 1, quantile, probs = c(0.025, 0.5, 0.975))
    df <- as.data.frame(t(q.infections.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> pivot_longer(-year, names_to = "quartile")
  }
  return(df)
}

incidence.by.state.df <- function(stan_data, PrevbyYearInf, state = 1, start.yr = 1978, annual = FALSE) {
  nyr <- stan_data$nquar / 4

  ### Total infections
  ## - Those undx at end of surveillance period (in any states)
  ## - Those Dx with HIV and no CD4 (PrevbyYearInf["HIV.Dx.no.CD4",]) - actually the CD4 is "imputed"
  ## - Those Dx with HIV and CD4 (PrevbyYearInf["HIV.Dx.CD4",])
  ## - Those Dx with AIDS and reported (PrevbyYearInf["unreported.AIDS",])
  ## - Those Dx with AIDS and not reported (PrevbyYearInf["reported.AIDS",])

  if (state < 5) {
    infections.total <- sapply(PrevbyYearInf["undiagnosed", ], function(x) tail(head(x, n = state), n = 1)) +
      sapply(PrevbyYearInf["HIV.Dx.no.CD4", ], function(x) tail(head(x, n = state), n = 1)) +
      sapply(PrevbyYearInf["HIV.Dx.CD4", ], function(x) tail(head(x, n = state), n = 1))
  } else {
    infections.total <- simplify2array(PrevbyYearInf["unreported.AIDS", ]) +
      simplify2array(PrevbyYearInf["reported.AIDS", ])
  }

  if (annual == FALSE) {
    q.infections.total <- apply(infections.total, 1, quantile, probs = c(0.025, 0.5, 0.975))
    df <- as.data.frame(t(q.infections.total))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> pivot_longer(-year, names_to = "quartile")
  } else {
    infections.years <- apply(infections.total, 2, function(x) tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum))
    q.infections.years <- apply(infections.years, 1, quantile, probs = c(0.025, 0.5, 0.975))
    df <- as.data.frame(t(q.infections.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> pivot_longer(-year, names_to = "quartile")
  }
  return(df)
}

### Undiagnosed prevalence

undiag.prev.df <- function(stan_data, PrevbyYearInf, InfDistnbyYearDx, start.yr = 1978, annual = FALSE) {
  nyr <- stan_data$nquar / 4
  # start.yr <- 1995
  ## Undiagnosed prevalence is equal to the number of people
  ## that are undiagnosed by the end of the surveillance period
  undiagnosed.prevalence <- sapply(PrevbyYearInf["undiagnosed", ], function(x) cumsum(colSums(x)))

  ### Plus all those that were diagnosed (AIDS or HIV) # why??
  ### between the infection time and dx time
  ### these are counted by the prev_fun(1)
  undiagnosed.prevalence.AIDS <- sapply(InfDistnbyYearDx["AIDS.Dx.Infect.Dist", ], function(x) t(prev_fun(x)))
  undiagnosed.prevalence.HIV <- sapply(InfDistnbyYearDx["HIV.Dx.Infect.Dist", ], function(x) t(prev_fun1(x)))
  undiagnosed.prevalence <- undiagnosed.prevalence + undiagnosed.prevalence.AIDS + undiagnosed.prevalence.HIV

  q.undiagnosed.prevalence <- apply(undiagnosed.prevalence, 1, quantile, probs = c(0.025, 0.5, 0.975))[, 1:(stan_data$nquar)]
  df <- as.data.frame(t(q.undiagnosed.prevalence))
  df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
  df <- df |> pivot_longer(-year, names_to = "quartile")

  if (annual == TRUE) {
    df <- df |>
      filter(year %% 1 == 0.75) |>
      mutate(year = year - 0.75)
  }

  return(df)
}

### Diagnosis probabilities
diag.prob.df <- function(stan_data, fit_mat, start.yr = 1978) {
  nyr <- stan_data$nquar / 4

  d1.ind <- grep("d\\[1,", colnames(fit_mat))
  d2.ind <- grep("d\\[2,", colnames(fit_mat))
  d3.ind <- grep("d\\[3,", colnames(fit_mat))
  d4.ind <- grep("d\\[4,", colnames(fit_mat))

  ### Quantiles of diagnosis
  d1.q <- apply(fit_mat[, d1.ind], 2, quantile, probs = c(0.025, 0.5, 0.975))
  d2.q <- apply(fit_mat[, d2.ind], 2, quantile, probs = c(0.025, 0.5, 0.975))
  d3.q <- apply(fit_mat[, d3.ind], 2, quantile, probs = c(0.025, 0.5, 0.975))
  d4.q <- apply(fit_mat[, d4.ind], 2, quantile, probs = c(0.025, 0.5, 0.975))

  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  d1.q <- t(d1.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4>500") |>
    pivot_longer(c(-year, -strata), names_to = "quartile")
  d2.q <- t(d2.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4 350-500") |>
    pivot_longer(c(-year, -strata), names_to = "quartile")
  d3.q <- t(d3.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4 200-350") |>
    pivot_longer(c(-year, -strata), names_to = "quartile")
  d4.q <- t(d4.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4<200") |>
    pivot_longer(c(-year, -strata), names_to = "quartile")

  df <- rbind(d1.q, d2.q, d3.q, d4.q)
  return(df)
}

diagnoses.df <- function(stan_data, PrevbyYearInf, start.yr = 1978, annual = FALSE) {
  nyr <- (stan_data$nquar) / 4

  ### Plus all those that were diagnosed (AIDS or HIV)
  ### between the infection time and dx time
  ### these are counted by the prev_fun(1)

  diagnosed.AIDS <- sapply(PrevbyYearInf["AIDS.diagnoses", ], function(x) colSums(x))

  diagnosed.HIV <- sapply(PrevbyYearInf["HIV.diagnoses", ], function(x) colSums(x))

  q.diagnosed <- apply(diagnosed.HIV + diagnosed.AIDS, 1, quantile, probs = c(0.025, 0.5, 0.975))[, 1:(stan_data$nquar)]
  df <- as.data.frame(t(q.diagnosed))
  df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
  df <- df |> pivot_longer(-year, names_to = "quartile")

  if (annual == TRUE) {
    df <- df |>
      filter(year %% 1 == 0.75) |>
      mutate(year = year - 0.75)
  }

  return(df)
}
