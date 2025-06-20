## ## ################ PROCESSING THE AGE-DEPENDENT POSTERIOR OUTPUT FROM STAN

### This is a general code to clean the back-calculation output produced by STAN
### It can deal with different lengths of data (e.g. time periods)
### And quarterly and yearly models
### with and without age-dependent diagnoses probabilities

here::i_am("r/process_ad.R")

library(here)
library(rstan)
library(Rcpp)
library(RcppArmadillo)

# load the functions to process the output
source(here("r/ad_functions.R"))

## Load results from stan model
# load(here("data/postsim_ad_rita_2011_2023.RData"))
load(here("data/postsim_ad_2011_2023.RData"))

############################################################# DATA PROCESSING ####################################################################

out <- list()

## Useful constants

## nt is the the number of period considered
## i.e. nquar if model quarterly
## or nyr if model yearly
## is.null allows distinction between quarter and years
nt <- ifelse(is.null(model$stan_data$nquar) == FALSE, model$stan_data$nquar, model$stan_data$nyr)
start_year <- 2011
nage <- model$stan_data$nage

### Specifications
# out$data <- 1 # tentatively set data.rm <- 0
out$nt <- nt
# out$data.id <- data.id

### Fit summaries
fit.mat <- as.matrix(fit)
fit.summ <- summary(fit)[[1]]

### Flags for the model considered
## Is this quarterly? Yearly?
## AgeDx Probs? or not
quar.flag <- ifelse(is.null(model$stan_data$nquar) == FALSE, TRUE, FALSE)
ageprobs.flag <- ifelse(length(grep("d\\[", colnames(fit.mat))) == nt * 4, FALSE, TRUE)
ageprobs.spl.flag <- ifelse(length(grep("d\\[", colnames(fit.mat))) == nt * 4 * nage, TRUE, FALSE)
delta.flag <- ifelse(length(grep("delta1", colnames(fit.mat))) > 0, TRUE, FALSE)
rita.flag <- ifelse(length(grep("d\\[1,5", colnames(fit.mat))) > 0, TRUE, FALSE)
urep.flag <- ifelse(("under_rep" %in% colnames(fit.mat)) == TRUE, TRUE, FALSE)

out$quar.flag <- quar.flag
out$ageprobs.flag <- ageprobs.flag
out$ageprobs.spl.flag <- ageprobs.spl.flag
out$rita.flag <- rita.flag
out$urep.flag <- urep.flag

### General information
## Lack of convergence
out$conv.chk <- length(fit.summ[, 10][fit.summ[, 10] > 1.05]) == 0

## Divergent stuff
sampler.params <- get_sampler_params(fit, inc_warmup = FALSE)
out$n.divg <- sum(sapply(sampler.params, function(x) x[, 5]))

## Time taken
out$time <- max(apply(get_elapsed_time(fit), 1, sum))

## Infections (yearly)
## Splines - beta: spline pars, gamma = X*beta
## infs matrix with log-infs for every simulation
## has same structure for GP and splines

beta.ind <- grep("beta", colnames(fit.mat))
lambda.ind <- grep("lambda", colnames(fit.mat))
x <- model$stan_data$X
infs <- t(apply(fit.mat[, beta.ind], 1, function(y) x %*% y))

if (length(lambda.ind) > 1) {
  lambda <- apply(fit.mat[, lambda.ind], 2, summary.fct)
} else { ## if length==1 not a matrix thus does not work
  lambda <- summary.fct(fit.mat[, lambda.ind])
}


### Posterior mean of smoothing parameter of interest
out$lambda <- lambda
### Posterior mean of alpha and sigma coefficients of interest
alpha.ind <- grep("alpha", colnames(fit.mat))
sigma2.ind <- grep("vardelta", colnames(fit.mat))

if (length(alpha.ind) != 0) out$alpha <- fit.mat[, alpha.ind]
if (length(sigma2.ind) != 0) out$sigma2 <- fit.mat[, sigma2.ind]

## Posterior mean for all time and ages
out$infs.all <- apply(infs, 2, summary.fct)

# Time profile and age-specific time profile
# yearly
out$infs <- infs.fct(infs, 15, 66, nt, type = "time", qt.flag = quar.flag)
out$infs.1524 <- infs.fct(infs, 15, 24, nt, type = "time", qt.flag = quar.flag)
out$infs.2534 <- infs.fct(infs, 25, 34, nt, type = "time", qt.flag = quar.flag)
out$infs.3544 <- infs.fct(infs, 35, 44, nt, type = "time", qt.flag = quar.flag)
out$infs.45 <- infs.fct(infs, 45, 66, nt, type = "time", qt.flag = quar.flag)

# quarterly
out$infs.qt <- infs.fct.qt(infs, 15, 66, nt, qt.flag = quar.flag)
out$infs.qt.1524 <- infs.fct.qt(infs, 15, 24, nt, qt.flag = quar.flag)
out$infs.qt.2534 <- infs.fct.qt(infs, 25, 34, nt, qt.flag = quar.flag)
out$infs.qt.3544 <- infs.fct.qt(infs, 35, 44, nt, qt.flag = quar.flag)
out$infs.qt.45 <- infs.fct.qt(infs, 45, 66, nt, qt.flag = quar.flag)

## ## incidence
## environment(inc.fct) <- environment()

## out$inc <- inc.fct(infs,"Tot", l.msmdata+1, nt, qt.flag = quar.flag)
## out$inc.1534 <- inc.fct(infs,"15-34", l.msmdata+1, nt, qt.flag = quar.flag)
## out$inc.3544 <- inc.fct(infs,"35-44", l.msmdata+1, nt, qt.flag = quar.flag)
## out$inc.45 <- inc.fct(infs,"45+", l.msmdata+1, nt, qt.flag = quar.flag)


## peaks in infections...
out$peak <- peak.fct(infs, 15, 66, nt, last.yr = 10, qt.flag = quar.flag)
out$peak.1524 <- peak.fct(infs, 15, 24, nt, last.yr = 10, qt.flag = quar.flag)
out$peak.2534 <- peak.fct(infs, 25, 34, nt, last.yr = 10, qt.flag = quar.flag)
out$peak.3544 <- peak.fct(infs, 35, 44, nt, last.yr = 10, qt.flag = quar.flag)
out$peak.45 <- peak.fct(infs, 45, 66, nt, last.yr = 10, qt.flag = quar.flag)

### trends in incidence in the last two years
out$trend.2 <- trend.fct(infs, 15, 66, nt, last.yr = 2, qt.flag = quar.flag)
out$trend.2.1524 <- trend.fct(infs, 15, 24, nt, last.yr = 2, qt.flag = quar.flag)
out$trend.2.2534 <- trend.fct(infs, 25, 34, nt, last.yr = 2, qt.flag = quar.flag)
out$trend.2.3544 <- trend.fct(infs, 35, 44, nt, last.yr = 2, qt.flag = quar.flag)
out$trend.2.45 <- trend.fct(infs, 45, 66, nt, last.yr = 2, qt.flag = quar.flag)

### trends in incidence in the last three years
out$trend.3 <- trend.fct(infs, 15, 66, nt, last.yr = 3, qt.flag = quar.flag)
out$trend.3.1524 <- trend.fct(infs, 15, 24, nt, last.yr = 3, qt.flag = quar.flag)
out$trend.3.2534 <- trend.fct(infs, 25, 34, nt, last.yr = 3, qt.flag = quar.flag)
out$trend.3.3544 <- trend.fct(infs, 35, 44, nt, last.yr = 3, qt.flag = quar.flag)
out$trend.3.45 <- trend.fct(infs, 45, 66, nt, last.yr = 3, qt.flag = quar.flag)

### Diagnoses probabilities
if (!ageprobs.flag) {
  d1.ind <- grep("d\\[1,", colnames(fit.mat))
  d2.ind <- grep("d\\[2,", colnames(fit.mat))
  d3.ind <- grep("d\\[3,", colnames(fit.mat))
  d4.ind <- grep("d\\[4,", colnames(fit.mat))

  out$d1.q <- apply(fit.mat[, d1.ind], 2, summary.fct)
  out$d2.q <- apply(fit.mat[, d2.ind], 2, summary.fct)
  out$d3.q <- apply(fit.mat[, d3.ind], 2, summary.fct)
  out$d4.q <- apply(fit.mat[, d4.ind], 2, summary.fct)
} else {
  if (!ageprobs.spl.flag) {
    d1.ind.y <- grep("d\\[1,1", colnames(fit.mat))
    d2.ind.y <- grep("d\\[1,2", colnames(fit.mat))
    d3.ind.y <- grep("d\\[1,3", colnames(fit.mat))
    d4.ind.y <- grep("d\\[1,4", colnames(fit.mat))
    d1.ind.my <- grep("d\\[2,1", colnames(fit.mat))
    d2.ind.my <- grep("d\\[2,2", colnames(fit.mat))
    d3.ind.my <- grep("d\\[2,3", colnames(fit.mat))
    d4.ind.my <- grep("d\\[2,4", colnames(fit.mat))
    d1.ind.mo <- grep("d\\[3,1", colnames(fit.mat))
    d2.ind.mo <- grep("d\\[3,2", colnames(fit.mat))
    d3.ind.mo <- grep("d\\[3,3", colnames(fit.mat))
    d4.ind.mo <- grep("d\\[3,4", colnames(fit.mat))
    d1.ind.o <- grep("d\\[4,1", colnames(fit.mat))
    d2.ind.o <- grep("d\\[4,2", colnames(fit.mat))
    d3.ind.o <- grep("d\\[4,3", colnames(fit.mat))
    d4.ind.o <- grep("d\\[4,4", colnames(fit.mat))

    out$d1.q.y <- apply(fit.mat[, d1.ind.y], 2, summary.fct)
    out$d2.q.y <- apply(fit.mat[, d2.ind.y], 2, summary.fct)
    out$d3.q.y <- apply(fit.mat[, d3.ind.y], 2, summary.fct)
    out$d4.q.y <- apply(fit.mat[, d4.ind.y], 2, summary.fct)
    out$d1.q.my <- apply(fit.mat[, d1.ind.my], 2, summary.fct)
    out$d2.q.my <- apply(fit.mat[, d2.ind.my], 2, summary.fct)
    out$d3.q.my <- apply(fit.mat[, d3.ind.my], 2, summary.fct)
    out$d4.q.my <- apply(fit.mat[, d4.ind.my], 2, summary.fct)
    out$d1.q.mo <- apply(fit.mat[, d1.ind.mo], 2, summary.fct)
    out$d2.q.mo <- apply(fit.mat[, d2.ind.mo], 2, summary.fct)
    out$d3.q.mo <- apply(fit.mat[, d3.ind.mo], 2, summary.fct)
    out$d4.q.mo <- apply(fit.mat[, d4.ind.mo], 2, summary.fct)
    out$d1.q.o <- apply(fit.mat[, d1.ind.o], 2, summary.fct)
    out$d2.q.o <- apply(fit.mat[, d2.ind.o], 2, summary.fct)
    out$d3.q.o <- apply(fit.mat[, d3.ind.o], 2, summary.fct)
    out$d4.q.o <- apply(fit.mat[, d4.ind.o], 2, summary.fct)

    if (rita.flag) {
      d5.ind.y <- grep("d\\[1,5", colnames(fit.mat))
      d5.ind.my <- grep("d\\[2,5", colnames(fit.mat))
      d5.ind.mo <- grep("d\\[3,5", colnames(fit.mat))
      d5.ind.o <- grep("d\\[4,5", colnames(fit.mat))

      out$d5.q.y <- apply(fit.mat[, d5.ind.y], 2, summary.fct)
      out$d5.q.my <- apply(fit.mat[, d5.ind.my], 2, summary.fct)
      out$d5.q.mo <- apply(fit.mat[, d5.ind.mo], 2, summary.fct)
      out$d5.q.o <- apply(fit.mat[, d5.ind.o], 2, summary.fct)
    }
  } else {
    d1.ind <- grep("d\\[1,", colnames(fit.mat))
    d2.ind <- grep("d\\[2,", colnames(fit.mat))
    d3.ind <- grep("d\\[3,", colnames(fit.mat))
    d4.ind <- grep("d\\[4,", colnames(fit.mat))
    d1.q <- apply(fit.mat[, d1.ind], 2, summary.fct)
    d2.q <- apply(fit.mat[, d2.ind], 2, summary.fct)
    d3.q <- apply(fit.mat[, d3.ind], 2, summary.fct)
    d4.q <- apply(fit.mat[, d4.ind], 2, summary.fct)
    out$d1.q <- array(d1.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
    out$d2.q <- array(d2.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
    out$d3.q <- array(d3.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
    out$d4.q <- array(d4.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
  }
}

if (delta.flag) {
  delta1.ind <- grep("delta1", colnames(fit.mat))
  delta2.ind <- grep("delta2", colnames(fit.mat))
  delta3.ind <- grep("delta3", colnames(fit.mat))
  delta4.ind <- grep("delta4", colnames(fit.mat))

  out$delta1.q <- apply(fit.mat[, delta1.ind], 2, summary.fct)
  out$delta2.q <- apply(fit.mat[, delta2.ind], 2, summary.fct)
  out$delta3.q <- apply(fit.mat[, delta3.ind], 2, summary.fct)
  out$delta4.q <- apply(fit.mat[, delta4.ind], 2, summary.fct)

  if (rita.flag) {
    delta5.ind <- grep("delta5", colnames(fit.mat))
    out$delta5.q <- apply(fit.mat[, delta5.ind], 2, summary.fct)
  }
} else {
  # deltas are the same for every age group so only derive once
  out$delta1.q <- apply(log((fit.mat[, d1.ind.y]) / (1 - fit.mat[, d1.ind.y])) - fit.mat[, alpha.ind][, 1], 2, summary.fct)
  out$delta2.q <- apply(log((fit.mat[, d2.ind.y]) / (1 - fit.mat[, d2.ind.y])) - fit.mat[, alpha.ind][, 2], 2, summary.fct)
  out$delta3.q <- apply(log((fit.mat[, d3.ind.y]) / (1 - fit.mat[, d3.ind.y])) - fit.mat[, alpha.ind][, 3], 2, summary.fct)
  out$delta4.q <- apply(log((fit.mat[, d4.ind.y]) / (1 - fit.mat[, d4.ind.y])) - fit.mat[, alpha.ind][, 4], 2, summary.fct)
}

## Loading Rcpp functions
## for post-simulating the epidemics
if (rita.flag) {
  Rcpp::sourceCpp(here("cpp/stancppfcts_quarter_rita.cpp"))
} else if (!ageprobs.flag & !quar.flag & !ageprobs.spl.flag) {
  Rcpp::sourceCpp(here("cpp/stancppfcts_simpleprev_yrinfsqt.cpp"))
} else if (ageprobs.flag & !quar.flag & !ageprobs.spl.flag) {
  Rcpp::sourceCpp(here("cpp/stancppfcts_simpleperv_agediag_qtinfs.cpp"))
} else if (!ageprobs.flag & quar.flag & !ageprobs.spl.flag) {
  Rcpp::sourceCpp(here("cpp/stancppfcts_quarter_simpleprev.cpp"))
} else if (ageprobs.flag & quar.flag & !ageprobs.spl.flag) {
  Rcpp::sourceCpp(here("cpp/stancppfcts_quarter_simpleprev_agediag.cpp"))
} else if (ageprobs.flag & quar.flag & ageprobs.spl.flag) {
  Rcpp::sourceCpp(here("cpp/stancppfcts_quarter_simpleprev_agediag_spl.cpp"))
} else {
  stop("Flags somehow uncorrectly specified")
}

## Need to set environment of d.mat.iter
## to be global environment (otherwise too many arguments to take care of ...)
environment(d.mat.iter) <- environment()

## GOODNESS OF FIT
exp.diags <- lapply(1:nrow(fit.mat), function(i) gof_iter_fct(exp(infs[i, ]), d.mat.iter(fit.mat, i, age.dx.flag = ageprobs.flag, age.dx.spl.flag = ageprobs.spl.flag, rita.flag, nt), model$stan_data$q, model$stan_data$init_prev, nt, nage))
exp.prev <- lapply(1:nrow(fit.mat), function(i) prev_iter_fct(exp(infs[i, ]), d.mat.iter(fit.mat, i, age.dx.flag = ageprobs.flag, age.dx.spl.flag = ageprobs.spl.flag, rita.flag, nt), model$stan_data$q, model$stan_data$init_prev, nt, nage))

### List to array (easier to handle)
exp.diags.arr <- simplify2array(exp.diags)
exp.prev.arr <- simplify2array(exp.prev)

### Prev.mat

if (rita.flag) {
  prev.mat.2011 <- rbind(
    apply(exp.prev.arr[1, , 1, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 2, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 3, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 4, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 5, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 6, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 7, ], 1, summary.fct)[2, ]
  )
} else if (start_year == 1995) {
  prev.mat.1995 <- rbind(
    apply(exp.prev.arr[1, , 1, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 2, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 3, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 4, ], 1, summary.fct)[2, ]
  )
  prev.mat.2011 <- rbind(
    apply(exp.prev.arr[65, , 1, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[65, , 2, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[65, , 3, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[65, , 4, ], 1, summary.fct)[2, ]
  )
} else if (start_year == 2011) {
  prev.mat.2011 <- rbind(
    apply(exp.prev.arr[1, , 1, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 2, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 3, ], 1, summary.fct)[2, ],
    apply(exp.prev.arr[1, , 4, ], 1, summary.fct)[2, ]
  )
}

### Expected number of CD4
nCD4 <- apply(model$stan_data$CD4, c(1, 2), sum)
if (rita.flag) {
  CD4.st1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 1, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  CD4.st2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 2, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  CD4.st3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 3, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  CD4.st4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 4, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  CD4.st5 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 5, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
} else {
  CD4.st1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 1, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
  CD4.st2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 2, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
  CD4.st3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 3, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
  CD4.st4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) nCD4 * (exp.diags.arr[, , 4, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
}

# to here:
### Remember that because of quarterly
### dynamics no CD4 dx can occur
### at times 5,9,13,... for age 1
### Thus replace all na by 0
if (quar.flag) {
  CD4.st1[which(is.na(CD4.st1) == TRUE)] <- 0
  CD4.st2[which(is.na(CD4.st2) == TRUE)] <- 0
  CD4.st3[which(is.na(CD4.st3) == TRUE)] <- 0
  CD4.st4[which(is.na(CD4.st4) == TRUE)] <- 0
  if (rita.flag) {
    CD4.st5[which(is.na(CD4.st5) == TRUE)] <- 0
  }
}

### Expected HIV diagnoses
if (rita.flag) {
  HIV <- exp.diags.arr[, , 1, ] + exp.diags.arr[, , 2, ] + exp.diags.arr[, , 3, ] + exp.diags.arr[, , 4, ] + exp.diags.arr[, , 5, ]
} else {
  HIV <- exp.diags.arr[, , 1, ] + exp.diags.arr[, , 2, ] + exp.diags.arr[, , 3, ] + exp.diags.arr[, , 4, ]
}

### Expected AIDS diagnoses
### Depending on whether under-reporting present
### (With real-data always)
if (urep.flag) {
  urep.vec <- fit.mat[, "under_rep"]
  AIDS <- simplify2array(lapply(1:nrow(fit.mat), function(i) urep.fct(urep.vec[i], exp.diags.arr[, , 5, i], urep.flag = urep.flag, qt.flag = quar.flag, nt = nt)))
} else if (rita.flag) {
  AIDS <- exp.diags.arr[, , 6, ]
} else {
  AIDS <- exp.diags.arr[, , 5, ]
}

### Storing summary statistics
if (rita.flag) {
  out$HIV <- gof.fct.rita("HIV", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$HIV.1524 <- gof.fct.rita("HIV", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$HIV.2534 <- gof.fct.rita("HIV", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$HIV.3544 <- gof.fct.rita("HIV", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$HIV.45 <- gof.fct.rita("HIV", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)

  out$AIDS <- gof.fct.rita("AIDS", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$AIDS.1524 <- gof.fct.rita("AIDS", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$AIDS.2534 <- gof.fct.rita("AIDS", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$AIDS.3544 <- gof.fct.rita("AIDS", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$AIDS.45 <- gof.fct.rita("AIDS", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)

  out$CD4.1 <- gof.fct.rita("CD4.st1", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.2 <- gof.fct.rita("CD4.st2", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.3 <- gof.fct.rita("CD4.st3", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.4 <- gof.fct.rita("CD4.st4", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)

  out$CD4.1.1524 <- gof.fct.rita("CD4.st1", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.1.2534 <- gof.fct.rita("CD4.st1", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.1.3544 <- gof.fct.rita("CD4.st1", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.1.45 <- gof.fct.rita("CD4.st1", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)

  out$CD4.2.1524 <- gof.fct.rita("CD4.st2", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.2.2534 <- gof.fct.rita("CD4.st2", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.2.3544 <- gof.fct.rita("CD4.st2", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.2.45 <- gof.fct.rita("CD4.st2", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)

  out$CD4.3.1524 <- gof.fct.rita("CD4.st3", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.3.2534 <- gof.fct.rita("CD4.st3", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.3.3544 <- gof.fct.rita("CD4.st3", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.3.45 <- gof.fct.rita("CD4.st3", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)

  out$CD4.4.1524 <- gof.fct.rita("CD4.st4", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.4.2534 <- gof.fct.rita("CD4.st4", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.4.3544 <- gof.fct.rita("CD4.st4", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
  out$CD4.4.45 <- gof.fct.rita("CD4.st4", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, CD4.st5, nt, quar.flag)
} else {
  out$HIV <- gof.fct("HIV", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$HIV.1524 <- gof.fct("HIV", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$HIV.2534 <- gof.fct("HIV", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$HIV.3544 <- gof.fct("HIV", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$HIV.45 <- gof.fct("HIV", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)

  out$AIDS <- gof.fct("AIDS", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$AIDS.1524 <- gof.fct("AIDS", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$AIDS.2534 <- gof.fct("AIDS", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$AIDS.3544 <- gof.fct("AIDS", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$AIDS.45 <- gof.fct("AIDS", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)

  out$CD4.1 <- gof.fct("CD4.st1", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.2 <- gof.fct("CD4.st2", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.3 <- gof.fct("CD4.st3", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.4 <- gof.fct("CD4.st4", 15, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)

  out$CD4.1.1524 <- gof.fct("CD4.st1", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.1.2534 <- gof.fct("CD4.st1", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.1.3544 <- gof.fct("CD4.st1", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.1.45 <- gof.fct("CD4.st1", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)

  out$CD4.2.1524 <- gof.fct("CD4.st2", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.2.2534 <- gof.fct("CD4.st2", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.2.3544 <- gof.fct("CD4.st2", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.2.45 <- gof.fct("CD4.st2", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)

  out$CD4.3.1524 <- gof.fct("CD4.st3", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.3.2534 <- gof.fct("CD4.st3", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.3.3544 <- gof.fct("CD4.st3", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.3.45 <- gof.fct("CD4.st3", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)

  out$CD4.4.1524 <- gof.fct("CD4.st4", 15, 24, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.4.2534 <- gof.fct("CD4.st4", 25, 34, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.4.3544 <- gof.fct("CD4.st4", 35, 44, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
  out$CD4.4.45 <- gof.fct("CD4.st4", 45, 66, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, quar.flag)
}

out$prev <- prev.fct(exp.prev.arr, "all", 15, 66)

out$prev.st1 <- prev.fct(exp.prev.arr, 1, 15, 66)
out$prev.st2 <- prev.fct(exp.prev.arr, 2, 15, 66)
out$prev.st3 <- prev.fct(exp.prev.arr, 3, 15, 66)
out$prev.st4 <- prev.fct(exp.prev.arr, 4, 15, 66)

out$prev.1524 <- prev.fct(exp.prev.arr, "all", 15, 24)
out$prev.2534 <- prev.fct(exp.prev.arr, "all", 25, 34)
out$prev.3544 <- prev.fct(exp.prev.arr, "all", 35, 44)
out$prev.45 <- prev.fct(exp.prev.arr, "all", 45, 66)

out$prev.st1.1524 <- prev.fct(exp.prev.arr, 1, 15, 24)
out$prev.st1.2534 <- prev.fct(exp.prev.arr, 1, 25, 34)
out$prev.st1.3544 <- prev.fct(exp.prev.arr, 1, 35, 44)
out$prev.st1.45 <- prev.fct(exp.prev.arr, 1, 45, 66)

out$prev.st2.1524 <- prev.fct(exp.prev.arr, 2, 15, 24)
out$prev.st2.2534 <- prev.fct(exp.prev.arr, 2, 25, 34)
out$prev.st2.3544 <- prev.fct(exp.prev.arr, 2, 35, 44)
out$prev.st2.45 <- prev.fct(exp.prev.arr, 2, 45, 66)

out$prev.st3.1524 <- prev.fct(exp.prev.arr, 3, 15, 24)
out$prev.st3.2534 <- prev.fct(exp.prev.arr, 3, 25, 34)
out$prev.st3.3544 <- prev.fct(exp.prev.arr, 3, 35, 44)
out$prev.st3.45 <- prev.fct(exp.prev.arr, 3, 45, 66)

out$prev.st4.1524 <- prev.fct(exp.prev.arr, 4, 15, 24)
out$prev.st4.2534 <- prev.fct(exp.prev.arr, 4, 25, 34)
out$prev.st4.3544 <- prev.fct(exp.prev.arr, 4, 35, 44)
out$prev.st4.45 <- prev.fct(exp.prev.arr, 4, 45, 66)

if (rita.flag) {
  out$prev.st5 <- prev.fct(exp.prev.arr, 5, 15, 66)
  out$prev.st5.1524 <- prev.fct(exp.prev.arr, 5, 15, 24)
  out$prev.st5.2534 <- prev.fct(exp.prev.arr, 5, 25, 34)
  out$prev.st5.3544 <- prev.fct(exp.prev.arr, 5, 35, 44)
  out$prev.st5.45 <- prev.fct(exp.prev.arr, 5, 45, 66)

  out$prev.st6 <- prev.fct(exp.prev.arr, 6, 15, 66)
  out$prev.st6.1524 <- prev.fct(exp.prev.arr, 6, 15, 24)
  out$prev.st6.2534 <- prev.fct(exp.prev.arr, 6, 25, 34)
  out$prev.st6.3544 <- prev.fct(exp.prev.arr, 6, 35, 44)
  out$prev.st6.45 <- prev.fct(exp.prev.arr, 6, 45, 66)

  out$prev.st7 <- prev.fct(exp.prev.arr, 7, 15, 66)
  out$prev.st7.1524 <- prev.fct(exp.prev.arr, 7, 15, 24)
  out$prev.st7.2534 <- prev.fct(exp.prev.arr, 7, 25, 34)
  out$prev.st7.3544 <- prev.fct(exp.prev.arr, 7, 35, 44)
  out$prev.st7.45 <- prev.fct(exp.prev.arr, 7, 45, 66)
}
### Getting data posterior predictive distributions
## set.seed to ensure consistent results
N.Poiss <- 300

HIV.m <- apply(HIV, c(1, 2), mean)

set.seed(407)
HIV.p <- apply(HIV.m, c(1, 2), function(x) rpois(N.Poiss, lambda = x))

AIDS.m <- apply(AIDS, c(1, 2), mean)

set.seed(407)
AIDS.p <- apply(AIDS.m, c(1, 2), function(x) rpois(N.Poiss, lambda = x))

if (rita.flag) {
  p1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 1, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  p2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 2, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  p3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 3, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  p4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 4, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
  p5 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 5, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i] + exp.diags.arr[, , 5, i]))))
} else {
  p1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 1, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
  p2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 2, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
  p3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 3, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
  p4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) (exp.diags.arr[, , 4, i] / (exp.diags.arr[, , 1, i] + exp.diags.arr[, , 2, i] + exp.diags.arr[, , 3, i] + exp.diags.arr[, , 4, i]))))
}

if (quar.flag) {
  p1[which(is.na(p1) == TRUE)] <- 0
  p2[which(is.na(p2) == TRUE)] <- 0
  p3[which(is.na(p3) == TRUE)] <- 0
  p4[which(is.na(p4) == TRUE)] <- 0
  if (rita.flag) {
    p5[which(is.na(p5) == TRUE)] <- 0
  }
}

p1.m <- apply(p1, c(1, 2), mean)
p2.m <- apply(p2, c(1, 2), mean)
p3.m <- apply(p3, c(1, 2), mean)
p4.m <- apply(p4, c(1, 2), mean)
if (rita.flag) {
  p5.m <- apply(p4, c(1, 2), mean)
}

if (rita.flag) {
  CD4.Poiss <- array(0, dim = c(nt, nage, 5, N.Poiss))
} else {
  CD4.Poiss <- array(0, dim = c(nt, nage, 4, N.Poiss))
}

if (quar.flag && rita.flag) { ### if quarterly there cannot be CD4 for a=1, at t=5, t=9,...
  ### as anyway nCD4[,1] = 0 then ignore and put 0
  for (t in 1:nt) {
    for (a in 2:nage) {
      set.seed(407)
      CD4.Poiss[t, a, , ] <- rmultinom(N.Poiss, size = nCD4[t, a], prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a], p5.m[t, a]))
    }
  }
} else if (quar.flag) {
  for (t in 1:nt) {
    for (a in 2:nage) {
      set.seed(407)
      CD4.Poiss[t, a, , ] <- rmultinom(N.Poiss, size = nCD4[t, a], prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a]))
    }
  }
} else {
  for (t in 1:nt) {
    for (a in 1:nage) {
      set.seed(407)
      CD4.Poiss[t, a, , ] <- rmultinom(N.Poiss, size = nCD4[t, a], prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a]))
    }
  }
}

### Selecting number of people in specific states
CD4.1.Poss.q <- CD4.Poiss[, , 1, ]
CD4.2.Poss.q <- CD4.Poiss[, , 2, ]
CD4.3.Poss.q <- CD4.Poiss[, , 3, ]
CD4.4.Poss.q <- CD4.Poiss[, , 4, ]

if (rita.flag) {
  CD4.5.Poss.q <- CD4.Poiss[, , 5, ]
}

if (quar.flag) {
  qt.indx <- rep(1:(nt / 4), each = 4)
  call <- "function(x) tapply(rowSums(x), qt.indx, sum)"
} else {
  call <- "rowSums"
}

out$HIV.Poiss <- apply(apply(HIV.p, 1, eval(parse(text = call))), 1, summary.fct)
out$HIV.1524.Poiss <- apply(apply(HIV.p[, , 1:10], 1, eval(parse(text = call))), 1, summary.fct)
out$HIV.2534.Poiss <- apply(apply(HIV.p[, , 11:20], 1, eval(parse(text = call))), 1, summary.fct)
out$HIV.3544.Poiss <- apply(apply(HIV.p[, , 21:30], 1, eval(parse(text = call))), 1, summary.fct)
out$HIV.45.Poiss <- apply(apply(HIV.p[, , 31:52], 1, eval(parse(text = call))), 1, summary.fct)

out$AIDS.Poiss <- apply(apply(AIDS.p, 1, eval(parse(text = call))), 1, summary.fct)
out$AIDS.1524.Poiss <- apply(apply(AIDS.p[, , 1:10], 1, eval(parse(text = call))), 1, summary.fct)
out$AIDS.2534.Poiss <- apply(apply(AIDS.p[, , 11:20], 1, eval(parse(text = call))), 1, summary.fct)
out$AIDS.3544.Poiss <- apply(apply(AIDS.p[, , 21:30], 1, eval(parse(text = call))), 1, summary.fct)
out$AIDS.45.Poiss <- apply(apply(AIDS.p[, , 31:52], 1, eval(parse(text = call))), 1, summary.fct)

out$CD4.1.Mult <- apply(apply(CD4.1.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.2.Mult <- apply(apply(CD4.2.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.3.Mult <- apply(apply(CD4.3.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.4.Mult <- apply(apply(CD4.4.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)

out$CD4.1.1524.Mult <- apply(apply(CD4.1.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.1.2534.Mult <- apply(apply(CD4.1.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.1.3544.Mult <- apply(apply(CD4.1.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.1.45.Mult <- apply(apply(CD4.1.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

out$CD4.2.1524.Mult <- apply(apply(CD4.2.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.2.2534.Mult <- apply(apply(CD4.2.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.2.3544.Mult <- apply(apply(CD4.2.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.2.45.Mult <- apply(apply(CD4.2.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

out$CD4.3.1524.Mult <- apply(apply(CD4.3.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.3.2534.Mult <- apply(apply(CD4.3.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.3.3544.Mult <- apply(apply(CD4.3.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.3.45.Mult <- apply(apply(CD4.3.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

out$CD4.4.1524.Mult <- apply(apply(CD4.4.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.4.2534.Mult <- apply(apply(CD4.4.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.4.3544.Mult <- apply(apply(CD4.4.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
out$CD4.4.45.Mult <- apply(apply(CD4.4.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

if (quar.flag) {
  call <- "rowSums"

  out$HIV.Quar.Poiss <- apply(apply(HIV.p, 1, eval(parse(text = call))), 1, summary.fct)
  out$HIV.1524.Quar.Poiss <- apply(apply(HIV.p[, , 1:10], 1, eval(parse(text = call))), 1, summary.fct)
  out$HIV.2534.Quar.Poiss <- apply(apply(HIV.p[, , 11:20], 1, eval(parse(text = call))), 1, summary.fct)
  out$HIV.3544.Quar.Poiss <- apply(apply(HIV.p[, , 21:30], 1, eval(parse(text = call))), 1, summary.fct)
  out$HIV.45.Quar.Poiss <- apply(apply(HIV.p[, , 31:52], 1, eval(parse(text = call))), 1, summary.fct)

  out$AIDS.Quar.Poiss <- apply(apply(AIDS.p, 1, eval(parse(text = call))), 1, summary.fct)
  out$AIDS.1524.Quar.Poiss <- apply(apply(AIDS.p[, , 1:10], 1, eval(parse(text = call))), 1, summary.fct)
  out$AIDS.2534.Quar.Poiss <- apply(apply(AIDS.p[, , 11:20], 1, eval(parse(text = call))), 1, summary.fct)
  out$AIDS.3544.Quar.Poiss <- apply(apply(AIDS.p[, , 21:30], 1, eval(parse(text = call))), 1, summary.fct)
  out$AIDS.45.Quar.Poiss <- apply(apply(AIDS.p[, , 31:52], 1, eval(parse(text = call))), 1, summary.fct)

  out$CD4.1.Quar.Mult <- apply(apply(CD4.1.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.2.Quar.Mult <- apply(apply(CD4.2.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.3.Quar.Mult <- apply(apply(CD4.3.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.4.Quar.Mult <- apply(apply(CD4.4.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)

  out$CD4.1.1524.Quar.Mult <- apply(apply(CD4.1.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.1.2534.Quar.Mult <- apply(apply(CD4.1.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.1.3544.Quar.Mult <- apply(apply(CD4.1.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.1.45.Quar.Mult <- apply(apply(CD4.1.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

  out$CD4.2.1524.Quar.Mult <- apply(apply(CD4.2.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.2.2534.Quar.Mult <- apply(apply(CD4.2.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.2.3544.Quar.Mult <- apply(apply(CD4.2.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.2.45.Quar.Mult <- apply(apply(CD4.2.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

  out$CD4.3.1524.Quar.Mult <- apply(apply(CD4.3.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.3.2534.Quar.Mult <- apply(apply(CD4.3.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.3.3544.Quar.Mult <- apply(apply(CD4.3.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.3.45.Quar.Mult <- apply(apply(CD4.3.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)

  out$CD4.4.1524.Quar.Mult <- apply(apply(CD4.4.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.4.2534.Quar.Mult <- apply(apply(CD4.4.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.4.3544.Quar.Mult <- apply(apply(CD4.4.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.4.45.Quar.Mult <- apply(apply(CD4.4.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)
}

if (rita.flag) {
  out$CD4.5.Mult <- apply(apply(CD4.5.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.1524.Mult <- apply(apply(CD4.5.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.2534.Mult <- apply(apply(CD4.5.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.3544.Mult <- apply(apply(CD4.5.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.45.Mult <- apply(apply(CD4.5.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)
  
  out$CD4.5.Quar.Mult <- apply(apply(CD4.5.Poss.q, 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.1524.Quar.Mult <- apply(apply(CD4.5.Poss.q[, 1:10, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.2534.Quar.Mult <- apply(apply(CD4.5.Poss.q[, 11:20, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.3544.Quar.Mult <- apply(apply(CD4.5.Poss.q[, 21:30, ], 3, eval(parse(text = call))), 1, summary.fct)
  out$CD4.5.45.Quar.Mult <- apply(apply(CD4.5.Poss.q[, 31:52, ], 3, eval(parse(text = call))), 1, summary.fct)
}
# ### Time to diagnosis stuff
# if(!ageprobs.flag & !quar.flag & !ageprobs.spl.flag){
#   stop("Time to diagnosis not coded up yet")
# }else if(ageprobs.flag & !quar.flag & !ageprobs.spl.flag){
#   stop("Time to diagnosis not coded up yet")
# }else if(!ageprobs.flag & quar.flag & !ageprobs.spl.flag){
#   Rcpp::sourceCpp('HIV backcalculation/Age dependency/MSMdata16/TimeToDxFct_aidiags.cpp')
# }else if(ageprobs.flag & quar.flag & !ageprobs.spl.flag){
#   Rcpp::sourceCpp('HIV backcalculation/Age dependency/MSMdata16/TimeToDxFct_acdiags.cpp')
# }else if(ageprobs.flag & quar.flag & ageprobs.spl.flag){
#   Rcpp::sourceCpp('HIV backcalculation/Age dependency/MSMdata16/TimeToDxFct_adiags.cpp')
# }else{
#   stop("Flags somehow uncorrectly specified")
# }
# d.ind <- grep("d\\[",colnames(fit.mat))
# d.mat <- fit.mat[,d.ind]
#
# ### the call depends on model... maybe include nac in fct
# t.to.d <- t_to_d(nt, nage, 4, n_ind_sim=20, d.mat, model$stan_data$q)
# snaps <- snap(nt, nage, 4, n_ind_sim=20, t_max=150,  d.mat, model$stan_data$q)
#
# t.to.d.dist <- apply(t.to.d[[1]], c(1,2), quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
#
# t.to.d.dist[t.to.d.dist==(nt+1)] <- NA
# ### Scaling quarters to years
# out$t.to.d.dist  <- t.to.d.dist/4
#
# snaps.dist <- apply(snaps[[1]], c(1,2), quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
#
# ### Scaling quarters to years
# out$snaps.dist  <- snaps.dist/4
#

# save(out, file = here("data/postproc_ad_rita_2011_2023.RData"))
save(out, file = here("data/postproc_ad_2011_2023.RData"))
# save(prev.mat.2011, file = here("data/inits_2011.RData"))
