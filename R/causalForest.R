init.causalForest <- function(formula, data, treatment, weights=F, cost=F, num.trees, ncov_sample) { 
  num.obs <- nrow(data)
  trees <- vector("list", num.trees)
  inbag <- matrix(0, num.obs, num.trees) 
  cov_sample <- matrix(0, num.trees, ncov_sample)
  inbag.Est <- matrix(0, num.obs, num.trees)
  nameall_sample <- matrix(0, num.trees, ncov_sample + 2) #2 end cols for y,w,no tau_true
  fsample <- vector("list", num.trees)
  causalForestobj <- list(trees = trees, formula = formula, data = data, treatment = treatment, weights = weights, cost = cost, ntree = num.trees, inbag = inbag, cov_sample = cov_sample, fsample = fsample, nameall_sample = nameall_sample, inbag.Est = inbag.Est) 
  class(causalForestobj) <- "causalForest" 
  return(causalForestobj)
} 

predict.causalForest <- function(forest, newdata, weight.type = NULL, predict.all = FALSE, type="vector") {
  if (!inherits(forest, "causalForest")) stop("Not a legitimate \"causalForest\" object")
  individual <- sapply(forest$trees, function(tree.fit) {
    predict(tree.fit, newdata = newdata, type = "vector")
  })
  
  #replace sapply with a loop if needed
  print(dim(individual))
  if (is.null(weight.type)) {
    print("use non-weighted prediction")
    aggregate <- rowMeans(individual)
  } else {
    if (weight.type == 0){
      aggregate <- matrixStats::rowWeightedMeans(individual, w = forest$d0weights)
    } else if(weight.type == 2){
      aggregate <- matrixStats::rowWeightedMeans(individual, w = forest$d2weights)
    }
    
  }
  
  if (predict.all) {
    list(aggregate = aggregate, individual = individual)
  } else {
    aggregate
  }
}

causalForest <- function(formula, data, treatment,  
                         na.action = na.causalTree, 
                         split.Rule="CT", double.Sample =T, split.Honest=T, split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                         propensity, control, split.alpha = 0.5, cv.alpha = 0.5,
                         sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = .5,
                         mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
                         cost=F, weights=F,ncolx,ncov_sample) {
  
  # do not implement subset option of causalTree, that is inherited from rpart but have not implemented it here yet
  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2:length(vars)]
  treatmentdf <- data.frame(treatment)
  if(class(data)[1]=="data.table"){
    treatmentdt <- data.table(treatment)
    datax<-data[, ..x]
    datay<-data[, y, with=FALSE]
    data <- cbind(datax, datay, treatmentdt)
  }else if (class(data)=="data.frame"){
    data <- data[, c(x, y)]
    data <- cbind(data, treatmentdf)
  }
  
  num.obs <-nrow(data)
  
  causalForest.obj <- init.causalForest(formula=formula, data=data, treatment=treatment, weights=weights, cost=cost, num.trees=num.trees, ncov_sample=ncov_sample)
  
  sample.size <- min(sample.size.total, num.obs)
  if (double.Sample) {
    train.size <- round(sample.size.train.frac*sample.size)
    est.size <- sample.size - train.size  
  }
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)

    print(paste0("The number of total observation is ", num.obs))
    print(paste0("The length of set a is ", length(full.idx)))
    
    if(double.Sample) {
      train.idx <- full.idx[1:train.size]
      reestimation.idx <- full.idx[(train.size+1):sample.size]
    }
    
    #randomize over the covariates for splitting (both train and reestimation)
    cov_sample<-sample.int(ncolx)
    cov_sample<-cov_sample[1:ncov_sample]
    
    #modify the y=f(x) equation accordingly for this tree
    #and modify the colnames
    fsample<-""
    nextx<-""
    nameall_sample<-c()
    for (ii in 1:(ncov_sample)){
      nextxindex <- cov_sample[ii]
      nextx <- x[[nextxindex]]
      if (ii==1) {
        fsample <-nextx
        name<- nextx
      }
      if (ii>1) {
        fsample <- paste0(fsample, "+", nextx)
        name <- c(name, nextx)
      }
    }
    nameall_sample <- c(name, y, "w") #, "tau_true")
    
    #store this var subset for each tree (need it during testing/predict stage)
    causalForest.obj$cov_sample[tree.index,]<-cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.obj$nameall_sample[tree.index, ]<-nameall_sample
    causalForest.obj$fsample[[tree.index]]<-fsample
    
    if(class(data)[1]=="data.table"){
      if (double.Sample) {
        dataTree <- data.table(data[train.idx, ])
        dataEstim <- data.table(data[reestimation.idx, ])
      }else{
        dataTree <- data.table(data[full.idx, ])
      }
      #pick relevant covariates for tree
      treeRange<-c(cov_sample, (ncolx+1):ncol(dataTree))
      estimRange<-c(cov_sample, (ncolx+1):ncol(dataEstim))
      dataTree <- dataTree[, ..treeRange]
      if (double.Sample) {
        dataEstim <- dataEstim[, ..estimRange]
      }
    }else if(class(data)=="data.frame"){
      if (double.Sample) {
        dataTree <- data.frame(data[train.idx, ])
        dataEstim <- data.frame(data[reestimation.idx, ])
      }else{
        dataTree <- data.frame(data[full.idx, ])
      }
      #pick relevant covariates for tree
      dataTree <- dataTree[, c(cov_sample, (ncolx+1):ncol(dataTree))]
      if (double.Sample) {
        dataEstim <- dataEstim[, c(cov_sample, (ncolx+1):ncol(dataEstim))]
      }
    }
    
    
    #change colnames to reflect the sampled cols
    names(dataTree) <- nameall_sample
    if(double.Sample) {
      names(dataEstim) <- nameall_sample
    }
    
    #save rdata for debug here, if needed
    formula<-paste0(y, "~", fsample)
    
    if (double.Sample) {
      tree.obj <- honest.causalTree(formula, data = dataTree, 
                                    treatment = treatmentdf[train.idx, ], 
                                    est_data=dataEstim, est_treatment= treatmentdf[reestimation.idx, ],
                                    split.Rule=split.Rule, split.Honest= split.Honest, split.Bucket=split.Bucket, 
                                    bucketNum = bucketNum, 
                                    bucketMax = bucketMax, cv.option="CT", cv.Honest=T, 
                                    minsize = nodesize, 
                                    split.alpha = 0.5, cv.alpha = 0.5, xval=0, 
                                    HonestSampleSize=est.size, cp=0)

      print(colnames(dataTree))
      print(dim(dataTree))
    }else {
      tree.obj <- causalTree(formula, data = dataTree, treatment = treatmentdf[full.idx,], 
                             na.action = na.causalTree, 
                             split.Rule=split.Rule, split.Honest= split.Honest, split.Bucket=split.Bucket, 
                             bucketNum = bucketNum,
                             bucketMax = bucketMax, cv.option="CT", cv.Honest=T,
                             x = FALSE, y = TRUE,
                             split.alpha = 0.5, cv.alpha = 0.5, cv.gamma=0.5, split.gamma=0.5)
      
    }
    
    causalForest.obj$trees[[tree.index]] <- tree.obj
    causalForest.obj$inbag[full.idx, tree.index] <- 1
    if (double.Sample) {causalForest.obj$inbag.Est[reestimation.idx, tree.index] <- 1}
  }
  return (causalForest.obj)
}

consistentcausalForest <- function(formula, data, treatment,  
                                    na.action = na.causalTree, 
                                    split.Rule="CT", double.Sample =T, split.Honest=T, split.Bucket=F, bucketNum = 5,
                                    bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                                    propensity, control, split.alpha = 0.5, cv.alpha = 0.5,
                                    sample.size.total = floor(nrow(data) * 0.4), sample.size.train.frac = .5,
                                    mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
                                    cost=F, weights=F, ncolx, ncov_sample) {
  
  # do not implement subset option of causalTree, that is inherited from rpart but have not implemented it here yet
  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2:length(vars)]
  treatmentdf <- data.frame(treatment)
  if(class(data)[1] == "data.table") {
    treatmentdt <- data.table(treatment)
    datax <- data[, ..x]
    datay <- data[, y, with = FALSE]
    data <- cbind(datax, datay, treatmentdt)
  } else if(class(data) == "data.frame") {
    data <- data[, c(x, y)]
    data <- cbind(data, treatmentdf)
  }
  
  num.obs <- nrow(data)
  
  causalForest.obj <- init.causalForest(formula = formula, data = data, treatment = treatment, weights = weights, cost = cost, num.trees = num.trees*2, ncov_sample = ncov_sample)
  
  sample.size <- min(sample.size.total, num.obs)

  # set sample size for tree a and b
  sample.size.a <- floor(sample.size/2)
  sample.size.b <- sample.size - sample.size.a

  if (double.Sample) {
    train.size.a <- round(sample.size.train.frac * sample.size.a)
    est.size.a <- sample.size.a - train.size.a
    
    train.size.b <- round(sample.size.train.frac * sample.size.b)
    est.size.b <- sample.size.b - train.size.b
  
    eval.size <- num.obs - sample.size.b - sample.size.a
  }
  
  print("Building trees ...")
  
  d0.weights <- c()
  d2.weights <- c()
  tree.index <- 1
  # for (tree.index in 1:num.trees) {
  while(tree.index <= num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    # a set idx
    full.idx.a <- sample.int(num.obs, sample.size.a, replace = FALSE)

    # b set idx
    full.idx.b <- sample(seq(1, num.obs)[!(seq(1, num.obs) %in% full.idx.a)], sample.size.b, replace = FALSE)

    # eval set idx
    eval.idx <- seq(1, num.obs)[!(seq(1, num.obs) %in% c(full.idx.a, full.idx.b))]

    print("Intersection of a and b is ")
    print(intersect(full.idx.a, full.idx.b))
    print(setdiff(union(union(full.idx.a, full.idx.b) , eval.idx), seq(1, num.obs)))

    print(paste0("The number of total observation is ", num.obs))
    print(paste0("The length of set a is ", length(full.idx.a)))
    print(paste0("The length of set b is ", length(full.idx.b)))
    print(paste0("The length of set eval is ", length(eval.idx)))

    if(double.Sample) {
      train.idx.a <- full.idx.a[1:train.size.a]
      reestimation.idx.a <- full.idx.a[(train.size.a + 1):sample.size.a]

      train.idx.b <- full.idx.b[1:train.size.b]
      reestimation.idx.b <- full.idx.b[(train.size.b + 1):sample.size.b]
    }
    
    #randomize over the covariates for splitting (both train and reestimation)
    cov_sample <- sample.int(ncolx)
    cov_sample <- cov_sample[1:ncov_sample]
    
    #modify the y=f(x) equation accordingly for this tree
    #and modify the colnames
    fsample <- ""
    nextx <- ""
    nameall_sample <- c()
    for (ii in 1:(ncov_sample)) {
      nextxindex <- cov_sample[ii]
      nextx <- x[[nextxindex]]
      if (ii == 1) {
        fsample <- nextx
        name <- nextx
      }
      if (ii > 1) {
        fsample <- paste0(fsample,"+",nextx)
        name <- c(name, nextx)
      }
    }
    nameall_sample <- c(name, y, "w") #, "tau_true")
    
    #store this var subset for each tree (need it during testing/predict stage)
    causalForest.obj$cov_sample[tree.index,] <- cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.obj$nameall_sample[tree.index,] <- nameall_sample
    causalForest.obj$fsample[[tree.index]] <- fsample
    
    if(class(data)[1] == "data.table") {
      if (double.Sample) {
        dataTree.a <- data.table(data[train.idx.a, ])
        dataEstim.a <- data.table(data[reestimation.idx.a, ])

        # test set
        dataTree.b <- data.table(data[train.idx.b, ])
        dataEstim.b <- data.table(data[reestimation.idx.b, ])

        #pick relevant covariates for tree
        treeRange.a <- c(cov_sample, (ncolx + 1):ncol(dataTree.a))
        treeRange.b <- c(cov_sample, (ncolx + 1):ncol(dataTree.b))
        estimRange.a <- c(cov_sample, (ncolx + 1):ncol(dataEstim.a))
        estimRange.b <- c(cov_sample, (ncolx + 1):ncol(dataEstim.b))
        dataTree.a <- dataTree.a[, ..treeRange.a]
        dataTree.b <- dataTree.b[, ..treeRange.b]
        if (double.Sample) {
          dataEstim.a <- dataEstim.a[, ..estimRange.a]
          dataEstim.b <- dataEstim.b[, ..estimRange.b]
        }

        dataeval <- data.table(data[eval.idx, ])
        treeRange.eval <- c(cov_sample, (ncolx + 1):ncol(dataeval))
        dataeval <- dataeval[, ..treeRange.eval]

      }else{
        dataTree <- data.table(data[full.idx, ])

        #pick relevant covariates for tree
        treeRange<-c(cov_sample, (ncolx+1):ncol(dataTree))
        estimRange<-c(cov_sample, (ncolx+1):ncol(dataEstim))
        dataTree <- dataTree[, ..treeRange]
        if (double.Sample) {
          dataEstim <- dataEstim[, ..estimRange]
        }
      }
      
    }else if(class(data) == "data.frame"){
      if (double.Sample) {
        dataTree.a <- data.frame(data[train.idx.a, ])
        dataEstim.a <- data.frame(data[reestimation.idx.a, ])

        # test set
        dataTree.b <- data.frame(data[train.idx.b, ])
        dataEstim.b <- data.frame(data[reestimation.idx.b, ])

        #pick relevant covariates for tree
        dataTree.a <- dataTree.a[, c(cov_sample, (ncolx + 1):ncol(dataTree.a))]
        dataTree.b <- dataTree.b[, c(cov_sample, (ncolx + 1):ncol(dataTree.b))]
        if (double.Sample) {
          dataEstim.a <- dataEstim.a[, c(cov_sample, (ncolx + 1):ncol(dataEstim.a))]
          dataEstim.a <- dataEstim.b[, c(cov_sample, (ncolx + 1):ncol(dataEstim.b))]
        }

        dataeval <- data.frame(data[eval.idx, ])
        treeRange.eval <- c(cov_sample, (ncolx + 1):ncol(dataeval))
        dataeval <- dataeval[, c(cov_sample, (ncolx + 1):ncol(dataeval))]

      } else {
        dataTree <- data.frame(data[full.idx, ])

        #pick relevant covariates for tree
        dataTree <- dataTree[, c(cov_sample, (ncolx+1):ncol(dataTree))]
        if (double.Sample) {
          dataEstim <- dataEstim[, c(cov_sample, (ncolx+1):ncol(dataEstim))]
        }
      }

    }
    
    #change colnames to reflect the sampled cols
    if(double.Sample){
      names(dataTree.a) <- nameall_sample
      names(dataTree.b) <- nameall_sample
      names(dataEstim.a) <- nameall_sample
      names(dataEstim.b) <- nameall_sample
      names(dataeval) <- nameall_sample
    } else {
      names(dataTree) <- nameall_sample
    }
    
    
    #save rdata for debug here, if needed
    formula <- paste0(y, "~", fsample)
    
    if (double.Sample) {
      # we train two tree with given covariates using different part of data
      tree.obj.a <- honest.causalTree(formula, data = dataTree.a, 
                                    treatment = dataTree.a$w, #treatmentdf[train.idx.a, ], 
                                    est_data = dataEstim.a, est_treatment = dataEstim.a$w, #treatmentdf[reestimation.idx.a, ],
                                    split.Rule = split.Rule, split.Honest = split.Honest, split.Bucket = split.Bucket, 
                                    bucketNum = bucketNum, 
                                    bucketMax = bucketMax, cv.option = "CT", cv.Honest = T, 
                                    minsize = nodesize, 
                                    split.alpha = 0.5, cv.alpha = 0.5, xval = 0, 
                                    HonestSampleSize = est.size.a, cp = 0)

      print(colnames(dataTree.a))
      print(dim(dataTree.a))
      tree.obj.b <- honest.causalTree(formula, data = dataTree.b, 
                                    treatment = dataTree.b$w, #treatmentdf[train.idx.b, ], 
                                    est_data = dataEstim.b, est_treatment = dataTree.b$w, #treatmentdf[reestimation.idx.b, ],
                                    split.Rule = split.Rule, split.Honest = split.Honest, split.Bucket = split.Bucket, 
                                    bucketNum = bucketNum, 
                                    bucketMax = bucketMax, cv.option = "CT", cv.Honest = T, 
                                    minsize = nodesize, 
                                    split.alpha = 0.5, cv.alpha = 0.5, xval = 0, 
                                    HonestSampleSize = est.size.b, cp = 0)

      # # make prediction using two trees respectively.
      # a.prediction <- predict(tree.obj.a, newdata = dataeval, type = "vector")
      # b.prediction <- predict(tree.obj.b, newdata = dataeval, type = "vector")

      # # test the correlation between the prediction of a and b tree
      # cor.res <- cor.test(a.prediction, b.prediction)
 
      # print(colnames(dataeval)[!(colnames(dataeval) %in% c("y", "w"))])
      dataeval_x <- dataeval[, colnames(dataeval)[!(colnames(dataeval) %in% c("y", "w"))]]
      # print(dataeval_x[1:5, ])

      # tree.obj.a$nodes
      nodes.a <- tryCatch({
        get.tree.struct(tree.obj.a)
      }, error = function(e){
        NA
      })
      
      # tree.obj.b$
      nodes.b <- tryCatch({
        get.tree.struct(tree.obj.b)
      }, error = function(e){
        NA
      })

      # tree.obj.a$nodes <- get.tree.struct(tree.obj.a)
      # tree.obj.b$nodes <- get.tree.struct(tree.obj.b)

      if(is.na(nodes.a) || is.na(nodes.b)){
        # return(list(tree.obj.a, tree.obj.b, dataeval_x))
        next # a situation I can't resolve yet, same count but "improve" improve, we can't know there is a new node... 
      } else {
        tree.obj.a$nodes <- nodes.a
        tree.obj.b$nodes <- nodes.b
      }

      # return(list(tree.obj.a, tree.obj.b, dataeval_x))

      d0 <- obtain_hamming_distance(obtain_bs(tree.obj.a, dataeval_x), obtain_bs(tree.obj.b, dataeval_x))
      d2 <- obtain_d2_distance(tree.obj.a, tree.obj.b, dataeval_x)

      d0.weights <- c(d0.weights, d0)
      d2.weights <- c(d2.weights, d2)

    }else {
      tree.obj <- causalTree(formula, data = dataTree, treatment = treatmentdf[full.idx, ], 
                             na.action = na.causalTree, 
                             split.Rule = split.Rule, split.Honest = split.Honest, split.Bucket = split.Bucket, 
                             bucketNum = bucketNum,
                             bucketMax = bucketMax, cv.option = "CT", cv.Honest = T,
                             x = FALSE, y = TRUE,
                             split.alpha = 0.5, cv.alpha = 0.5, cv.gamma = 0.5, split.gamma = 0.5)
      
    }
    
    if (double.Sample){
      # if (cor.res$p.value > 0.05){
      #   next
      # } else {
      #   causalForest.obj$trees[[tree.index]] <- tree.obj.a
      #   causalForest.obj$inbag[full.idx.a, tree.index] <- 1
      #   if (double.Sample) {causalForest.obj$inbag.Est[reestimation.idx.a, tree.index] <- 1}
      # }
      
      # how to handle the tree object
      causalForest.obj$trees[[tree.index]] <- tree.obj.a
      causalForest.obj$inbag[full.idx.a, tree.index] <- 1

      causalForest.obj$trees[[tree.index + num.trees]] <- tree.obj.b
      causalForest.obj$inbag[full.idx.b, (tree.index + num.trees)] <- 1

      # print(d2.weights)
      d0.weights <- 1/d0.weights
      d2.weights <- 1/d2.weights
      print(d0.weights)
      print(d2.weights)

      # causalForest.obj$d0weights <- c(d0.weights/sum(d0.weights), d0.weights/sum(d0.weights))
      # causalForest.obj$d2weights <- c(d2.weights/sum(d2.weights), d2.weights/sum(d2.weights))

      causalForest.obj$d0weights <- c(d0.weights, d0.weights)
      causalForest.obj$d2weights <- c(d2.weights, d2.weights)

      # print(causalForest.obj$d2weights)

      if (double.Sample) {
        causalForest.obj$inbag.Est[reestimation.idx.a, tree.index] <- 1
        causalForest.obj$inbag.Est[reestimation.idx.b, (tree.index + num.trees)] <- 1
      }

    } else {
      causalForest.obj$trees[[tree.index]] <- tree.obj
      causalForest.obj$inbag[full.idx, tree.index] <- 1
      if (double.Sample) {causalForest.obj$inbag.Est[reestimation.idx, tree.index] <- 1}
    }

    tree.index <- tree.index + 1

  }
  
  print("Finished.")
  return(causalForest.obj)
}

propensityForest <- function(formula, data, treatment,  
                             na.action = na.causalTree, 
                             split.Rule="CT", split.Honest=T, split.Bucket=F, bucketNum = 5,
                             bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                             propensity = mean(treatment), control, split.alpha = 0.5, cv.alpha = 0.5,  
                             sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = 1,
                             mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
                             ncolx, ncov_sample) {
  
  # do not implement subset option of causalTree, inherited from rpart
  # do not implement weights and costs yet
  
  if (sample.size.train.frac != 1) {
    print("warning: for propensity Forest, sample.size.train.frac should be 1; resetting to 1")
    sample.size.train.frac <- 1
  }
  
  num.obs <- nrow(data)
  
  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2:length(vars)]
  treatmentdf <- data.frame(treatment)
  if (class(data)[1] == "data.table") {
    treatmentdt <- data.table(treatment)
    datax <- data[, ..x]
    datay <- data[, y, with = FALSE]
    data <- cbind(datax, datay, treatmentdt)
  }else if (class(data) == "data.frame") {
    data <- data[, c(x, y)]
    data <- cbind(data, treatmentdf)
  }
  
  causalForest.hon <- init.causalForest(formula=formula, data=data, treatment=treatment, num.trees=num.trees, weights=F, cost=F,ncov_sample=ncov_sample)
  sample.size <- min(sample.size.total, num.obs)
  train.size <- round(sample.size.train.frac*sample.size)
  
  outcomename = as.character(formula[2])
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    train.idx <- full.idx[1:train.size]
    
    cov_sample<-sample.int(ncolx)
    cov_sample<-cov_sample[1:ncov_sample]
    
    #modify the y=f(x) equation accordingly for this tree
    #and modify the colnames
    fsample <- ""
    nextx <- ""
    nameall_sample <- c()
    for (ii in 1:(ncov_sample)) {
      nextxindex <- cov_sample[ii]
      nextx <- x[[nextxindex]]
      if (ii == 1) {
        fsample <- nextx
        name <- nextx
      }
      if (ii > 1) {
        fsample <- paste0(fsample, "+", nextx)
        name <- c(name, nextx)
      }
    }
    
    #nameall_sample <- c( name,"temptemp","y", "tau_true","treattreat")
    nameall_sample <- c( name,"temptemp", y, "treattreat")
    nameall_sample_save <- c( name,  y, "w") #, "tau_true")
    
    #store this var subset for each tree (need it during testing/predict stage)
    causalForest.hon$cov_sample[tree.index,] <- cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.hon$nameall_sample[tree.index,] <- nameall_sample_save
    causalForest.hon$fsample[[tree.index]] <- fsample
    
    # rename variables as a way to trick rpart into building the tree with all the object attributes considering the outcome variable as named
    # by the input formula, even though the tree itself is trained on w.  Note that we aren't saving out this propensity tree anyway, but if
    # we decided later to try to save out the propensity trees and do something directly with the propensity scores, we would need to do something
    # more tedious like estimate the propensity tree with the original names, and then edit the attributes to replace the treatment variable name
    # with the outcome variable name for the estimate part
    dataTree <- data.frame(data[train.idx, ])
    dataTree$treattreat <- treatmentdf[train.idx, ]
    names(dataTree)[names(dataTree) == outcomename] <- "temptemp"
    names(dataTree)[names(dataTree) == "treattreat"] <- outcomename
    
    
    if (class(data)[1] == "data.table"){
      #sample covariates and pick relevant covariates for tree
      treeRange <- c(cov_sample,(ncolx + 1):ncol(dataTree))
      dataTree <- dataTree[, ..treeRange]
    }else if(class(data) == "data.frame") {
      #pick relevant covariates for tree
      dataTree <- dataTree[, c(cov_sample, (ncolx + 1):ncol(dataTree))]
    }
    
    #change colnames to reflect the sampled cols
    names(dataTree) = nameall_sample
    # names(dataEstim)=nameall_sample
    formula <- paste(y, "~", fsample, sep="")
    #one options: estimate the propensity tree with anova so that it will be type "anova" when we re-estimate
    #here: replace elements of the rpart object to make it look like anova tree, so that we'll be able to properly predict with it later, etc.
    tree.propensity <- rpart(formula = formula, data = dataTree, method = "class", 
                             control = rpart.control(cp = 0, minbucket = nodesize))
    
    # make it look like a method="anova" tree 
    tree.propensity$method <- "anova"
    tree.propensity$frame$yval2 <- NULL
    tree.propensity$functions$print <- NULL
    
    # switch the names back in the data frame so that when we estimate treatment effects, will have the right outcome variables
    names(dataTree)[names(dataTree) == y] <- "treattreat"
    names(dataTree)[names(dataTree) == "temptemp"] <- y
    tree.treatment <- estimate.causalTree(object = tree.propensity, data = dataTree, treatment = dataTree$treattreat)
    
    causalForest.hon$trees[[tree.index]] <- tree.treatment
    causalForest.hon$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.hon)
}

#' return the tree structure as grf
#' 
#' @param tree the causal tree object we get from causalForest object
get.tree.struct <- function(tree) {
  # obtain the df with split node and corresponding split value
  tmp.1 <- as.data.frame(tree$splits)
  tmp.2 <- tree$frame
  row.names(tmp.1) <- NULL
  tmp.1$var <- attributes(tree$splits)$dimnames[[1]]
  # tmp.1.rm0 <- tmp.1[tmp.1$count != 0, ]
  # tmp.1.rm0.index1 <- tmp.1.rm0 %>% group_by(count) %>% filter(row_number() == 1)
  index_tmp <- c()
  for (i in 1:length(tmp.1$count)){
    if (i == 1){
      index_tmp <- c(index_tmp, 1)
      next
    } else {
      if ((tmp.1$count[i] != tmp.1$count[i-1]) && tmp.1$count[i] != 0){
        index_tmp <- c(index_tmp, i)
      }
    }
  }
  tmp.1.rm0.index1 <- tmp.1[index_tmp, ]

  tmp.2$dep_ind <- rownames(tmp.2)
  rownames(tmp.2) <- NULL
  tmp.2.rmleaf <- tmp.2[tmp.2$var != "<leaf>", ]
  tmp.1.rm0.index1$dep_ind <- tmp.2.rmleaf$dep_ind
  tmp.1.rm0.index1$var <- NULL
  join_res <- left_join(tmp.2, tmp.1.rm0.index1, by = "dep_ind")
  
  # return a grf tree like struct
  # initial the root node
  nodes <- list()
  
  for (i in 1:length(join_res$dep_ind)) {
    nodes <- list.append(nodes, get_node(join_res$dep_ind[i], join_res))
  }
  
  return(nodes)
}

#' support function to get the node information
#' @param dep_ind the depth index of the node
#' @param df dataframe we get from first part of get.tree.struct function
get_node <- function(dep_ind, df) {
  node <- list()
  dep_ind <- as.numeric(dep_ind)
  node$dep_ind <- dep_ind
  
  dep_ind_2 <- dep_ind * 2
  if ((dep_ind_2 %in% df$dep_ind) == TRUE){
    node$is_leaf <- FALSE
  } else {
    node$is_leaf <- TRUE
  }
  
  index <- which(df$dep_ind == dep_ind)
  node$split_variable <- df$var[index]
  node$split_value <- df$index[index]
  
  if (node$is_leaf == TRUE) {
    node$left_child <- NULL
    node$right_child <- NULL
    node$leaf_stats <- df$yval[index]
  } else {
    # lc_index <- i*2
    # rc_index <- i*2+1
    # 
    # node$left_child <- lc_index
    # node$right_child <- rc_index
    
    lc_index <- which(df$dep_ind == dep_ind_2)
    rc_index <- which(df$dep_ind == (dep_ind_2 + 1))
    
    node$left_child <- lc_index
    node$right_child <- rc_index
  }

  return(node)
}

#' Obtain binary strings.
#' @param tree casual tree object.
#' @param data the data used to train the casual forest.
obtain_bs <- function(tree, data) {
  variable_label <- rep(0, length(colnames(data)))
  names(variable_label) <- colnames(data)
  
  nodes_used <- unlist(tree$nodes)[which(names(unlist(tree$nodes)) == "split_variable")]
  nodes_used <- unique(nodes_used[!(nodes_used %in% c("<leaf>"))])
  
  nodes_count <- length(nodes_used)
  for (j in 1:nodes_count) {
    variable_label[names(variable_label) == nodes_used[j]] <- 1
  }
  
  return(variable_label)
}

#' Obtain Hamming distance between two trees.
#' @param binlist1 binary strings corresponding to tree 1.
#' @param binlist2 binary strings corresponding to tree 2.
obtain_hamming_distance <- function(binlist1, binlist2) {
  mismatch_count <- 0
  for (i in 1:length(binlist1)) {
      if (binlist1[i] != binlist2[i]) {
          mismatch_count <- mismatch_count + 1
      }
  }

  d_0_t1_t2 <- mismatch_count/length(binlist1)

  return(d_0_t1_t2)
}

#' Obtain d2 distance between two trees. (Deprecated)
#' @param tree1 casual tree object 1.
#' @param tree2 casual tree object 1.
#' @param data the data used to train the casual forest.
obtain_d2_distance <- function(tree1, tree2, data) {
  d2 <- 0
  for (i in 1:length(rownames(data))) { # check all pair of data.
      y_hat1 <- predict_tree(tree1$nodes[[1]], data[i, ], tree1)[1]
      y_hat2 <- predict_tree(tree2$nodes[[1]], data[i, ], tree2)[1]
      d2 <- d2 + (y_hat1 - y_hat2)^2/length(rownames(data))
  }
  return(d2)
}

#' Obtain d2 distance between two trees (multicore edition).
#' @param tree1 casual tree object 1.
#' @param tree2 casual tree object 1.
#' @param data the data used to train the casual forest.
obtain_d2_distance_mc <- function(tree1, tree2, data) {
  d2 <- foreach(i = 1:length(
      rownames(data)), .combine = "+") %dopar% { # check all pair of data.
      y_hat1 <- predict_tree(tree1$nodes[[1]], data[i, ], tree1)[1]
      y_hat2 <- predict_tree(tree2$nodes[[1]], data[i, ], tree2)[1]
      (y_hat1 - y_hat2)^2
  }
  d2 <- d2/length(rownames(data))
  return(d2)
}

#' Make prediction using single causal tree.
#' @param node node object in causal tree object.
#' @param sample single sample from the provided dataset (single patient).
#' @param tree Causal tree object.
predict_tree <- function(node, sample, tree) {
  # Check if this is leaf node
  if (node$is_leaf == TRUE) {
    return(node$leaf_stats)
  } else {
    # Check split varibale
    split_var <- node$split_variable
    split_val <- node$split_value
    if (sample[split_var] <= split_val) {
      return(predict_tree(tree$nodes[[node$left_child]], sample, tree))
    } else {
      return(predict_tree(tree$nodes[[node$right_child]], sample, tree))
    }
  }
}