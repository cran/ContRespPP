## ---- include = FALSE---------------------------------------------------------
#knitr::opts_chunk$set(
#  collapse = TRUE,
#  comment = "#>"
#)
old <- options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(ContRespPP)
data("exData")

## ----experimentDataArgs-------------------------------------------------------
# Design Matrix
X <- exData[,c(2:14)]
head(X)
# Response Column
Y <- as.matrix(exData[,1], ncol=1) # could also be a vector e.g., Y <- exData[,1]
# Observations Seen
n.seen <- 75

## ----priorArgs----------------------------------------------------------------
# Means for betas
beta.mean <- matrix(
  c(400, 50, 50, -25, -50, 100, 100, rep(0, 6)),
  ncol = 1
) 
# could also be vector: beta.mean <- c(400, 50, 50, -25, -50, 100, 100, rep(0, 6))

# Precisions for betas
beta.precision <- matrix(
  c(1/10000, 1/10000, 1/10000, 1/2500, 1/2500, 1/10000, 1/10000, rep(1/10000, 6)), 
  ncol = 1
) 
# could also be vector: beta.precision <- c(1/10000, 1/10000, 1/10000, 1/2500, 1/2500, 1/10000, 1/10000, rep(1/10000, 6))

# Precision for tau
shape <- 0.0001
rate <- 0.0001

## ----simParamArgs-------------------------------------------------------------
# Conditional Draws
b.sim <- 20000
b.burnin <- 2000

# Non-Conditional Draws
n.sim <- 1000
y.burnin <- 100

## ----analysisParamArgs--------------------------------------------------------
# Threshold Value
theta.t <- 0.8

# Metric Threshold Value
phi.0 <- 400

## -----------------------------------------------------------------------------
num.factors <- 5
num.factor.levels <- c(2, 2, 3, 2, 2)
likelihood.encountering <- c(1/2, 1/2, 4/9, 5/9, 1/3, 1/3, 1/3, 1/2, 1/2, 1/2, 1/2)

prob <- prob.creator(num.factors, num.factor.levels, likelihood.encountering, print.result = TRUE)

## -----------------------------------------------------------------------------
# Mission Probabilities (likelihood of encountering a factor level)	
prob <- matrix(
  c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 1/2, 1/2, 4/9, 5/9, 1/3, 1/3, 1/3, 1/2, 1/2, 1/2, 1/2), 
  ncol = 2, 
  dimnames = list(NULL, c("factor", "probability"))
)

## ----optParamArg1-------------------------------------------------------------
factor.no.2way <- c(3)

## ----optParamArg2-------------------------------------------------------------
colnames.pick <- c(
  "eta", "alpha", "beta", "omega2", "omega3", 
  "theta", "gamma", "alphabeta", "alphatheta", 
  "alphagamma", "betatheta", "betagamma", "thetagamma", "tau"
)

## ----gibbsSampler, eval = FALSE-----------------------------------------------
#  results <- gibbs.sampler(
#    X = X,
#    Y = Y,
#    n.seen = n.seen,
#    beta.mean = beta.mean,
#    beta.precision = beta.precision,
#    shape = shape,
#    rate = rate,
#    n.sim = n.sim,
#    y.burnin = y.burnin,
#    b.sim = b.sim,
#    b.burnin = b.burnin,
#    phi.0 = phi.0,
#    theta.t = theta.t,
#    prob = prob,
#    factor.no.2way = factor.no.2way,
#    colnames.pick = colnames.pick,
#    seed = 512
#  )

## ---- echo = FALSE------------------------------------------------------------
cat("Running Burn-in 1 of 100")
cat("Running Simulation 1 of 1000")
cat("Simulation Complete")

## ----readResults, include=FALSE-----------------------------------------------
# Internal data read in quietly
# (not run on each build due to long runtime)
results <- ContRespPP:::exResults[[1]] # this is the predictive rjags example

## ----ppResults----------------------------------------------------------------
results
results$pp

## ----posteriorResults---------------------------------------------------------
head(results$posterior)

## ----indicatorResults---------------------------------------------------------
head(results$indicator)

## ---- include=FALSE-----------------------------------------------------------
options(old)

