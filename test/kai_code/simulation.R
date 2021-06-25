library(grf)
library(MASS)
library(dplyr)
library(data.table)
# library(survival)
# library(survminer)
library(doParallel)

registerDoParallel(1)
args = commandArgs(trailingOnly=TRUE)

# source('survival.imputation.R')  # load impute.survival function
# source('hte.validation.R')     # load validation methods for HTE
# source('simulation.scenarios.R') # load different simulation functions
# source('causal_inference_models.R')

generate.covariates <- function(n, p = 10){
    X = matrix(rnorm(n*p), nrow = n, ncol = p)
    X[, seq(2, p, by = 2)] = (X[, seq(2, p, by = 2)] > 0)
    colnames(X) <- paste0('v', seq(p))
    return(X)
}

randomized.tx <- function(n){
    rbinom(n, size = 1, p = 0.5) # random treatment assignment
}

biased.tx <- function(mu.x, tau.x){
    diff <- mu.x - tau.x/2
    probs <- exp(diff) / (1 + exp(diff))
    # tx <- sapply(probs, function(i) rbinom(1, size = 1, p = i))
    tx <- rbinom(length(mu.x), 1, probs)
    tx
}

generate.scenario <- function(n, p, mu.fun, tau.fun, noise.sd = 1, randomized = TRUE, event.rate = 0.1){
    x <- generate.covariates(n, p)

    mu.x <- do.call(mu.fun, list(x))
    tau.x <- do.call(tau.fun, list(x))
    # mu.x <- f8(x)
    # tau.x <- f1(x)

    if(randomized == TRUE){
        tx <-randomized.tx(n)
    }else{
        tx <- biased.tx(mu.x, tau.x)
    }

    y <- mu.x + (tx - 1/2) * tau.x + rnorm(n, sd = noise.sd)
    # log.y <- mu.x + tx  * tau.x + rnorm(n, sd = noise.sd)

    # y <- exp(log.y)
    # y <- rexp(n, rate = 1/y.bar)
    # res <- generate.censored.outcome(y, event.rate)

    # max.censored <- max(res[, 1][res[, 2] == 0])
    # res[, 2][res[, 1] == max.censored] <- 1

    # imputed.y <- impute.survival(res[, 1], res[, 2])

    # d <- cbind(y, tau.x, tx, imputed.y, res[, 2], x)
    # d <- cbind(y, tau.x, tx, log(y), rep(0, length(y)), x)
    # colnames(d)[1:5] <- c('y',  'tau', 'tx', 'y.imputed', 'censor')

    d <- cbind(x, y, tx, tau.x)
    colnames(d)[(ncol(d)-2):ncol(d)] <- c('y',  'w', 'tau_true')
    return(d)
}