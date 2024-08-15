# Running this R file will install all of the packages required for the CD4 back-calculation project.
# The packages are installed using the pak package, which enables installing all the other packages in parallel, and installing any missing dependencies as well.
# The packages are installed in the user's home directory, so this script should be run as the user who will be running the back-calculation code.

install.packages("pak")

pak::pkg_install(
  c(
    "abind", "bayesplot", "BH", "data.table", "future.apply", "future.callr",
    "here", "httpgd", "jsonlite", "knitr", "languageserver", "lintr", "mgcv",
    "nanoparquet", "parallel", "patchwork", "posterior", "purrr", "Rcpp", "StanHeaders",
    "RcppArmadillo", "RcppEigen", "readstata13", "rmarkdown", "rstan", "scales", "styler",
    "tidyverse"
  )
)
