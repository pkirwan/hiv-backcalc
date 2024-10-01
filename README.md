# HIV back-calculation model

This repository contains the necessary code to implement the CD4-staged HIV back-calculation model, a Bayesian, discrete-time, multi-state model used to estimate HIV incidence from diagnosis data in the United Kingdom.

## Background

In the United Kingdom, data concerning HIV and AIDS diagnoses are collected centrally via epidemiological case surveillance, resulting in a time-series of observed counts of quarterly HIV and AIDS diagnoses. Provided an estimate of the incubation period (i.e. the time from exposure to symptoms) is available, this can be combined with the observed data on diagnoses to 'back-calculate' the incidence of infection. The general idea of back-calculation is that the time of symptom onset for an individual is equal to the sum of the time of exposure and the incubation period (Egan and Hall 2015). In continuous time this can be expressed by the following convolution:

$$
  d(t) = \int_{t_0}^t{h(s)f(t-s) \mathrm{d}s}
$$

where $d(t)$ is the rate of diagnosis at time $t$, $t_0$ is the starting time of the epidemic, $h(s)$ is the rate of new infections at time $s$, and $f(t-s)$ is the incubation distribution (Brookmeyer, Gail et al. 1994).

The occurrence of infections over time is typically modelled as a non-homogeneous Poisson process, with the number of new infections in a particular time interval being i.i.d. Poisson random variables with mean $h(s)$ (Becker et al. 1991; Rosenberg and Gail 1991). As linear combinations of Poisson random variables are themselves Poisson-distributed, so the number of new diagnoses is Poisson distributed, with mean $d(t)$ (Chiang 1980).

The back-calculation method was first developed to estimate HIV prevalence from AIDS cases, with the incubation distribution representing the time from HIV exposure to onset of AIDS (Brookmeyer and Gail 1988). With the introduction of HIV testing, individuals could be diagnosed before AIDS symptoms became apparent and the back-calculation model was therefore extended to consider pre-AIDS diagnoses. This extension introduced a time-dependent incubation distribution, $f(t-s \mid s)$, to capture changes in testing (Aalen, Farewell et al. 1997). More recent extensions of the model have incorporated HIV progression estimates from CD4 count data to better characterise the time between infection and diagnosis, and to estimate trends in the probability of HIV diagnosis by disease stage (Birrell, Chadborn et al. 2012; van Sighem, Nakagawa et al. 2015; Sweeting, De Angelis and Aalen 2005), and age group (Brizzi, Birrell, Plummer et al. 2019).

## Model implementation

The back-calculation model is implemented in the probabilistic programming language Stan (Carpenter et al. 2017) using the R interface `rstan` (Stan Development Team 2018). The models can be broadly categorised as age-independent or age-dependent, depending on whether the diagnosis probabilities and HIV progression probabilities are assumed to vary with age. The age-independent model is a generalisation of the model proposed by Birrell, Chadborn et al. (2012), and the age-dependent model is a generalisation of the model proposed by Brizzi, Birrell, Plummer et al. (2019).

### Project structure:

- `bug` is an empty folder used to store the jagam spline object;
- `cpp` contains C++ code used for post-hoc processing of the model output;
- `data` contains example data files with the required data structure for the age-independent and age-dependent models and initial conditions;
- `r` contains the R scripts used to fit the age-independent and age-dependent models and generate plots of estimates;
- `stan` contains the Stan model code.

### Running the models

For each model there are several R scripts that need to be run in sequence to fit the model, process the output, and generate plots of the results, these are:

- `r/install_packages.R` installs the necessary R packages;
- `r/run_ai.R` and `r/run_ad.R` run the age-independent/age-dependent model, using `rstan` with three chains, see below for the different Stan models which can be used;
- the `ShinyStan` package allows for inspection of the resulting trace plots and model fit;
- `r/process_ad.R` processes the fitted age-dependent model, using post-hoc methods to generate the required posterior datasets;
- `r/plots_ai.R` and `r/plots_ad.R` contain `ggplot2` code to visualise several quantities of interest;
- `r/future_pred.R` contains code to generate projected estimates of HIV incidence by extrapolating the fitted incidence splines forward.

For 100 calendar quarters of data, in a high-performance computing environment the age-independent model can be run in around 30 minutes, while the age-dependent model can take around 60 hours to run.

## Stan models

Several Stan files are available which differ in how they specify the sub-model for the quarterly diagnosis probabilities or incidence. These can be modelled as evolving over time as random walks, or smoothed using a Gaussian process or a thin-plate spline.

### Age-independent models

- `stan/AI_RW1978.stan` is the Stan model code for the age-independent model with incidence evolving as a random walk, from a starting point in 1978;
- `stan/AI_RW1995.stan` is the Stan model code for the age-independent model with incidence evolving as a random walk, from a starting point in 1995;
- `stan/AI_spl.stan` is the Stan model code for the age-independent model with incidence smoothed using an approximation to an optimal thin-plate spline containing 10 knots;
- `stan/AI_GP.stan` is the Stan model code for the age-independent model with incidence smoothed using a Gaussian process with mean zero and a quadratic exponential kernel.

### Age-dependent models

- `stan/ptens_quar_agediag1_nourep.stan` is the Stan model code for the age-dependent model with diagnosis probabilities governed by a state-specific random-walk with age- around an age- and state-specific linear trend.

$$
\begin{align*}
  \delta_{a,k,1} &= \delta_{0,k} + \sigma_{0,k}\epsilon_{k,1}\\
  \delta_{a,k,t} &= \delta_{a,k,t-1} + \alpha_{a,k} + \sigma_k \epsilon_{k,t}\\
  \epsilon_{k,t} &\sim \mathrm{N}\left(0,1\right)\\
  \alpha_{a,k} &\sim \mathrm{N}\left(0,1\right)\\
  \sigma_k &\sim \Gamma\left(1,8\right)
\end{align*}
$$

- `stan/ptens_quar_agediag2_nourep.stan` is the Stan model code for the age-dependent model with diagnosis probabilities evolving as age- and state-specific random walks.

$$
\begin{align*}
  \delta_{a,k,1} &= \delta_{0,k} + \sigma_{0,k}\epsilon_{a,k,1}\\
  \delta_{a,k,t} &= \delta_{a,k,t-1} + \sigma_{a,k}\epsilon_{a,k,t}\\
  \epsilon_{a,k,t} &\sim \mathrm{N}\left(0,1\right)\\
  \sigma_{a,k} &\sim \Gamma\left(1,8\right)
\end{align*}
$$

- `stan/ptens_quar_agediag3_nourep.stan` is the Stan model code for the age-dependent model with diagnosis probabilities smoothed using an approximation to an optimal thin-plate spline containing the same number of knots as the same type of spline used to smooth the incidence.

- `stan/ptens_quar_agediag4_nourep.stan` is the Stan model code for the age-dependent model with age-dependent diagnosis probabilities, (state-specific) random-walk with age-specific offset.

$$
\begin{align*}
  \delta_{a,k,1} &= \delta_{0,k} + \alpha_a + \sigma_{0,k}\epsilon_{a,k,1}\\
  \delta_{a,k,t} &= \delta_{a,k,t-1} + \sigma_k \epsilon_{k,t}\\
  \alpha_a &\sim \mathrm{N}\left(0,1\right)
\end{align*}
$$

- `stan/ptens_quar_agediag5_nourep.stan` is the Stan model code for the age-dependent model with age-dependent diagnosis probabilities, (state-specific) random-walk with age- and state-specific offsets (intercepts).

$$
\begin{align*}
  \delta_{a,k,1} &= \delta_{0,k} + \alpha_{a,k} + \sigma_{0,k}\epsilon_{k,1}\\
  \delta_{a,k,t} &= \delta_{a,k,t-1} + \sigma_k \epsilon_{k,t}
\end{align*}
$$

- `stan/ptens_quar_nourep_aidiag1.stan` is the Stan model code for the age-dependent model with age-independent diagnosis probabilities, (state-specific) random-walk around a linear trend.

$$
\begin{align*}
  \delta_{a,k,1} &= \delta_{0,k} + \sigma_{0,k}\epsilon_{k,1}\\
  \delta_{a,k,t} &= \delta_{a,k,t-1} + \alpha_k + \sigma_k\epsilon_{k,t}\\
  \alpha_k &\sim \mathrm{N}\left(0,1\right)
\end{align*}
$$

## References

Aalen, O. O., Farewell, V. T., De Angelis, D., Day, N. E. et al. (1997). ‘A Markov model for HIV disease progression including the effect of HIV diagnosis and treatment: application to AIDS prediction in England and Wales’. Stat. Med. 16 (19), pp. 2191–2210. doi: 10.1002/(sici)1097-0258(19971015)16:19<2191::aid-sim645>3.0.co;2-5.

Becker, N. G., Watson, L. F. and Carlin, J. B. (1991). ‘A method of non-parametric back-projection and its application to AIDS data’. Stat. Med. 10 (10), pp. 1527–1542. doi: 10.1002/sim.4780101005.

Birrell, P. J., Chadborn, T. R., Gill, O. N., Delpech, V. C. et al. (2012). ‘Estimating trends in incidence, time-to-diagnosis and undiagnosed prevalence using a CD4-based Bayesian back-calculation’. Stat. Commun. Infect. Dis. 4 (1). doi: 10.1515/1948-4690.1055.

Brizzi, F., Birrell, P. J., Plummer, M. T., Kirwan, P. et al. (2019). ‘Extending Bayesian back-calculation to estimate age and time specific HIV incidence’. Lifetime Data Anal. 25 (4),pp. 757–780. doi: 10.1007/s10985-019-09465-1.

Brookmeyer, R. and Gail, M. H. (1988). ‘A Method for Obtaining Short-Term Projections and Lower Bounds on the Size of the AIDS Epidemic’. J. Am. Stat. Assoc. 83 (402), pp. 301–308. doi: 10.1080/01621459.1988.10478599.

Brookmeyer, R., Gail, M. H., Medical Statistician Epidemiology and Biostatistics Program Mitchell H Gail and . Gail, M. H. (1994). AIDS Epidemiology: A Quantitative Approach. Oxford University Press. 354 pp.

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D. et al. (2017). ‘Stan: A probabilistic programming language’. J. Stat. Softw. 76 (1), pp. 1–32. doi: 10.18637/jss.v076.i01.

Chiang, C. L. (1980). An introduction to stochastic processes and their applications. RE Krieger Publishing Company New York.

Egan, J. R. and Hall, I. M. (2015). ‘A review of back-calculation techniques and their potential to inform mitigation strategies with application to non-transmissible acute infectious diseases’. J. R. Soc. Interface. 12 (106). doi: 10.1098/rsif.2015.0096.

Rosenberg, P. S. and Gail, M. H. (1991). ‘Backcalculation of flexible linear models of the human immunodeficiency virus infection curve’. J. R. Stat. Soc. Ser. C. Appl. Stat. 40 (2), pp. 269–282. doi: 10.2307/2347592.

Sweeting, M. J., De Angelis, D. and Aalen, O. O. (2005). ‘Bayesian back-calculation using a multi-state model with application to HIV’. Stat. Med. 24 (24), pp. 3991–4007. doi: 10.1002/sim.2432.

van Sighem, A., Nakagawa, F., De Angelis, D., Quinten, C. et al. (2015). ‘Estimating HIV Incidence, Time to Diagnosis, and the Undiagnosed HIV Epidemic Using Routine Surveillance Data’. Epidemiology. 26 (5), pp. 653–660. doi: 10.1097/EDE.0000000000000324.

## License

&copy; 2024. This work is openly licensed via [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
