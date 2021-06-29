# set eight scenarios for estimating HTE, details refers to Powers et al. 2018

f1 = function(x) {
    rep(0, nrow(x))
}

f2 = function(x) {
    5 * (x[, 1] > 1) - 5 * pnorm(-1)
    # 5 * (x[, 1] > 1) - 5
} 
    

f3 = function(x) {
    5 * x[, 1]
}

f4 = function(x) {
    1 * x[,2] * x[,4] * x[,6] + 
    2 * x[,2] * x[,4] * (1-x[,6]) +
    3 * x[,2] * (1-x[,4]) * x[,6] + 
    4 * x[,2] * (1-x[,4]) * (1-x[,6]) +
    5 * (1-x[,2]) * x[,4] * x[,6] + 
    6 * (1-x[,2]) * x[,4] * (1-x[,6]) +
    7 * (1-x[,2]) * (1-x[,4]) * x[,6] + 
    8 * (1-x[,2]) * (1-x[,4]) * (1-x[,6]) - 
    5 - (-0.5)
}

f5 = function(x) {
    x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] # remove + x[, 9] - 0.5 
    # x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2
}

f6 = function(x) {
    4 * (x[,1]>1) * (x[,3]>0) + 4 * (x[,5]>1) * (x[,7]>0) + 2 * x[,8] * x[,9] - 1 - (4*pnorm(-1)-1)
}

f7 = function(x){
  ( x[, 1]^3 + 
    x[, 2] + 
    x[, 3]^3 + 
    x[, 4] + 
    x[, 5]^3 + 
    x[, 6] + 
    x[, 7]^3 +
    x[, 8] + 
    x[, 9]^3) / sqrt(2)  # remove - 5 - (7/sqrt(2) - 5)
}

# covered in other scenario
f8 = function(x){
    (f4(x) + f5(x)) / sqrt(2)
}    

# covered in other scenario
f9 = function(x) {
    (x[,1])^2 - 1
}

d1 <- as.data.frame(generate.scenario(600, 400, 'f8', 'f1'))
d9 <- as.data.frame(generate.scenario(600, 400, 'f8', 'f1', randomized = FALSE))

d2 <- as.data.frame(generate.scenario(600, 400, 'f5', 'f2', noise.sd = 0.25))
d10 <- as.data.frame(generate.scenario(600, 400, 'f5', 'f2', noise.sd = 0.25, randomized = FALSE))

d3 <- as.data.frame(generate.scenario(800, 300, 'f4', 'f3'))
d11 <- as.data.frame(generate.scenario(800, 300, 'f4', 'f3', randomized = FALSE))

d4 <- as.data.frame(generate.scenario(800, 300, 'f7', 'f4', noise.sd = 0.25))
d12 <- as.data.frame(generate.scenario(800, 300, 'f7', 'f4', noise.sd = 0.25, randomized = FALSE))

d5 <- as.data.frame(generate.scenario(900, 200, 'f3', 'f5'))
d13 <- as.data.frame(generate.scenario(900, 200, 'f3', 'f5', randomized = FALSE))

d6 <- as.data.frame(generate.scenario(900, 200, 'f1', 'f6'))
d14 <- as.data.frame(generate.scenario(900, 200, 'f1', 'f6', randomized = FALSE))

d7 <- as.data.frame(generate.scenario(1000, 100, 'f2', 'f7', noise.sd = 4))
d15 <- as.data.frame(generate.scenario(1000, 100, 'f2', 'f7', noise.sd = 4, randomized = FALSE))

d8 <- as.data.frame(generate.scenario(300, 100, 'f6', 'f8', noise.sd = 4))
d16 <- as.data.frame(generate.scenario(1000, 100, 'f6', 'f8', noise.sd = 4, randomized = FALSE))

generate_simulation_function <- function(no_covar, coef_min = 0.3, coef_max = 0.5, constant = 0){
    # set.seed(seed)
    rand_coef = runif(no_covar, min = coef_min, max = coef_max)

    sd <- sqrt(sum(rand_coef^2))
    args <- sapply(1:no_covar, function(i) paste0(sqrt(5)*rand_coef[i]/sd, "*x[,", i, "]"))
    body <- paste0(args, collapse = '+')

    eval(parse(text = paste0('f <- function(x) {(' , body , ')}')))
    return(get('f'))
}
