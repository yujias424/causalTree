library(dplyr)
library(rlist)
library(foreach)

cl <- parallel::makeForkCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

source("./test/kai_code/simulation.R")
source("./test/kai_code/simulation.scenarios.R")

d1 <- generate.scenario(600, 400, 'f8', 'f1')

# create formulas for estimation
p <- 100
f <- ""
nextx <- ""
if (p>1) {
  for (ii in 1:(p-1)) {
    nextx <- paste("v",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
    f <- paste(f, nextx, "+", sep="")
  }
  f <- paste(f, "v", ii+1, sep="")
} else if (p==1) {
  f <- "v1"
}

# set training and testing data
dataTrain <- d16[1:600, ]
dataTest <- d16[601:1000, ]

ncov_sample <- floor(p/3)
ncov_sample <- p
ncolx <- p

cf <- causalForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                   split.Rule="CT", double.Sample = T, split.Honest=T,  split.Bucket=F, bucketNum = 5,
                   bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                   split.alpha = 0.5, cv.alpha = 0.5,
                   sample.size.total = floor(nrow(dataTrain) * 0.25), sample.size.train.frac = .5,
                   mtry = ceiling(ncol(dataTrain)/3), nodesize = 3, num.trees= 200,ncolx=ncolx,ncov_sample=ncov_sample
) 

cfpredtest <- predict(cf, newdata=dataTest, weight.type = NULL, type="vector")
plot(dataTest$tau_true, cfpredtest)
cor.test(dataTest$tau_true, cfpredtest)

source("./R/causalForest.R")

cf_const <- consistentcausalForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                                   split.Rule="CT", double.Sample = T, split.Honest=T,  split.Bucket=F, bucketNum = 5,
                                   bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                                   split.alpha = 0.5, cv.alpha = 0.5,
                                   sample.size.total = floor(nrow(dataTrain) * 0.5), sample.size.train.frac = .5,
                                   mtry = ceiling(ncol(dataTrain)/3), nodesize = 3, num.trees= 100,ncolx=ncolx, ncov_sample=ncov_sample
) 

cfpredtest <- predict(cf_const, newdata=dataTest, weight.type = 0, type="vector")
cor.test(dataTest$tau_true, cfpredtest)
plot(dataTest$tau_true, cfpredtest, ylim = c(0,3))

cfpredtest <- predict(cf_const, newdata=dataTest, weight.type = 2, type="vector")
cor.test(dataTest$tau_true, cfpredtest)
plot(dataTest$tau_true, cfpredtest, ylim = c(0,3))

cfpredtest <- predict(cf_const, newdata=dataTest, weight.type = NULL, type="vector")
cor.test(dataTest$tau_true, cfpredtest)
plot(dataTest$tau_true, cfpredtest, ylim = c(0,3))
