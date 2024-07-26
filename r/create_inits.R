here::i_am("r/create_inits.R")

library(here)

load(here("data/postsim_ad.RData"))

mcmc <- as.array(fit)

inds <- list(
    sigma = grep("lambda", dimnames(mcmc)[[3]], fixed = TRUE),
    beta = grep("beta", dimnames(mcmc)[[3]], fixed = TRUE),
    vardelta = grep("vardelta", dimnames(mcmc)[[3]], fixed = TRUE),
    alpha = grep("alpha", dimnames(mcmc)[[3]], fixed = TRUE)
)

mcmc <- mcmc[dim(mcmc)[1], , ]

logit <- function(p) log(p / (1 - p))

get.delta.raw <- function(x, delta0 = -3.2, sigma, sigma0 = 0.2, alpha) {
    delta <- logit(x)
    delta.raw <- (1 / sigma0) * (delta[1] - delta0 - alpha)
    for (t in 2:length(x)) {
        delta.raw[t] <- (1 / sigma) * (delta[t] - delta[t - 1])
    }
    delta.raw
}

inits1 <- list(
    sigma = 1 / sqrt(mcmc[1, inds$sigma]),
    beta = mcmc[1, inds$beta],
    vardelta = mcmc[1, inds$vardelta],
    alpha = mcmc[1, inds$alpha],
    delta_raw1 = get.delta.raw(mcmc[1, grep("d[1,1,", dimnames(mcmc)[[2]], fixed = TRUE)],
        sigma = mcmc[1, inds$vardelta[1]],
        alpha = mcmc[1, inds$alpha[1]]
    ),
    delta_raw2 = get.delta.raw(mcmc[1, grep("d[1,2,", dimnames(mcmc)[[2]], fixed = TRUE)],
        sigma = mcmc[1, inds$vardelta[2]],
        alpha = mcmc[1, inds$alpha[2]]
    ),
    delta_raw3 = get.delta.raw(mcmc[1, grep("d[1,3,", dimnames(mcmc)[[2]], fixed = TRUE)],
        delta0 = -3.0,
        sigma = mcmc[1, inds$vardelta[3]],
        alpha = mcmc[1, inds$alpha[3]]
    ),
    delta_raw4 = get.delta.raw(mcmc[1, grep("d[1,4,", dimnames(mcmc)[[2]], fixed = TRUE)],
        delta0 = -2.5,
        sigma = mcmc[1, inds$vardelta[4]],
        sigma0 = 0.3,
        alpha = mcmc[1, inds$alpha[4]]
    )
)

inits2 <- list(
    sigma = 1 / sqrt(mcmc[2, inds$sigma]),
    beta = mcmc[2, inds$beta],
    vardelta = mcmc[2, inds$vardelta],
    alpha = mcmc[2, inds$alpha],
    delta_raw1 = get.delta.raw(mcmc[2, grep("d[1,1,", dimnames(mcmc)[[2]], fixed = TRUE)],
        sigma = mcmc[2, inds$vardelta[1]],
        alpha = mcmc[2, inds$alpha[1]]
    ),
    delta_raw2 = get.delta.raw(mcmc[2, grep("d[1,2,", dimnames(mcmc)[[2]], fixed = TRUE)],
        sigma = mcmc[2, inds$vardelta[2]],
        alpha = mcmc[2, inds$alpha[2]]
    ),
    delta_raw3 = get.delta.raw(mcmc[2, grep("d[1,3,", dimnames(mcmc)[[2]], fixed = TRUE)],
        delta0 = -3.0,
        sigma = mcmc[2, inds$vardelta[3]],
        alpha = mcmc[2, inds$alpha[3]]
    ),
    delta_raw4 = get.delta.raw(mcmc[2, grep("d[1,4,", dimnames(mcmc)[[2]], fixed = TRUE)],
        delta0 = -2.5,
        sigma = mcmc[2, inds$vardelta[4]],
        sigma0 = 0.3,
        alpha = mcmc[2, inds$alpha[4]]
    )
)

inits3 <- list(
    sigma = 1 / sqrt(mcmc[3, inds$sigma]),
    beta = mcmc[3, inds$beta],
    vardelta = mcmc[3, inds$vardelta],
    alpha = mcmc[3, inds$alpha],
    delta_raw1 = get.delta.raw(mcmc[3, grep("d[1,1,", dimnames(mcmc)[[2]], fixed = TRUE)],
        sigma = mcmc[3, inds$vardelta[1]],
        alpha = mcmc[3, inds$alpha[1]]
    ),
    delta_raw2 = get.delta.raw(mcmc[3, grep("d[1,2,", dimnames(mcmc)[[2]], fixed = TRUE)],
        sigma = mcmc[3, inds$vardelta[2]],
        alpha = mcmc[3, inds$alpha[2]]
    ),
    delta_raw3 = get.delta.raw(mcmc[3, grep("d[1,3,", dimnames(mcmc)[[2]], fixed = TRUE)],
        delta0 = -3.0,
        sigma = mcmc[3, inds$vardelta[3]],
        alpha = mcmc[3, inds$alpha[3]]
    ),
    delta_raw4 = get.delta.raw(mcmc[3, grep("d[1,4,", dimnames(mcmc)[[2]], fixed = TRUE)],
        delta0 = -2.5,
        sigma = mcmc[3, inds$vardelta[4]],
        sigma0 = 0.3,
        alpha = mcmc[3, inds$alpha[4]]
    )
)

init.list <- list(inits1, inits2, inits3)

save(init.list, file = here("data/init_ad.RData"))
