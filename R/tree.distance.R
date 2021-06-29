#' return the tree structure as grf
#' 
#' @param tree the causal tree object we get from causalForest object
get.tree.struct <- function(tree) {
  # obtain the df with split node and corresponding split value
  tmp.1 <- as.data.frame(tree$splits)
  tmp.2 <- tree$frame
  row.names(tmp.1) <- NULL
  tmp.1$var <- attributes(tree$splits)$dimnames[[1]]
  tmp.1.rm0 <- tmp.1[tmp.1$count != 0, ]
  tmp.1.rm0.index1 <- tmp.1.rm0 %>% group_by(count) %>% filter(row_number() == 1)
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

predict_sample_node <- function(node, sample, tree) {
  # Check if this is leaf node
  if (node$is_leaf == TRUE) {
    return(node$dep_ind)
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

obtain_itij <- function(tree, sample_i, sample_j) {
  
  tree_predict_node_i <- predict_sample_node(tree$nodes[[1]], sample_i, tree)
  tree_predict_node_j <- predict_sample_node(tree$nodes[[1]], sample_j, tree)

  it <- ifelse(tree_predict_node_i == tree_predict_node_j, 1, 0)

  return(it)
}

obtain_d1_distance <- function(tree1, tree2, data){
  d1 <- 0
  
  allcomb <- combinat::combn(nrow(data), 2)
  
  for (k in 1:dim(allcomb)[2]) { 
    i <- allcomb[1, k]
    j <- allcomb[2, k]
    
    # IT1
    it1_ij <- obtain_itij(tree1, data[i, ], data[j, ])
    it2_ij <- obtain_itij(tree2, data[i, ], data[j, ])

    d1 <- d1 + abs(it1_ij - it2_ij)
  }
  
  print(d1)
  d1 <- d1/choose(nrow(data),2)
  return(d1)
}

obtain_d1_distance_mc <- function(tree1, tree2, data){
  d1 <- 0
  
  allcomb <- combinat::combn(nrow(data), 2)
  
  d1 <- foreach(k = 1:dim(allcomb)[2], .combine = "+") %dopar% {
    i <- allcomb[1, k]
    j <- allcomb[2, k]
    
    # IT1
    it1_ij <- obtain_itij(tree1, data[i, ], data[j, ])
    it2_ij <- obtain_itij(tree2, data[i, ], data[j, ])
    
    abs(it1_ij - it2_ij)
  }
  
  print(d1)
  d1 <- d1/choose(nrow(data),2)
  return(d1)
}

obtain_d1_distance_star <- function(tree1, tree2, data){
  d1 <- 0
  
  allcomb <- combinat::combn(nrow(data), 2)
  
  for (k in 1:dim(allcomb)[2]) { 
    i <- allcomb[1, k]
    j <- allcomb[2, k]
    
    # IT1
    it1_ij <- obtain_itij(tree1, data[i, ], data[j, ])
    it2_ij <- obtain_itij(tree2, data[i, ], data[j, ])

    d1 <- d1 + abs(it1_ij - it2_ij)
  }
  
  print(d1)
  d1 <- d1/choose(nrow(data),2)
  d1 <- d1 * obtain_hamming_distance(obtain_bs(tree1, data), obtain_bs(tree2, data))
  return(d1)
}

obtain_d1_distance_star_mc <- function(tree1, tree2, data){
  d1 <- 0
  
  allcomb <- combinat::combn(nrow(data), 2)
  
  d1 <- foreach(k = 1:dim(allcomb)[2], .combine = "+") %dopar% {
    i <- allcomb[1, k]
    j <- allcomb[2, k]
    
    # IT1
    it1_ij <- obtain_itij(tree1, data[i, ], data[j, ])
    it2_ij <- obtain_itij(tree2, data[i, ], data[j, ])
    
    abs(it1_ij - it2_ij)
  }
  
  print(d1)
  d1 <- d1/choose(nrow(data),2)
  d1 <- d1 * obtain_hamming_distance(obtain_bs(tree1, data), obtain_bs(tree2, data))
  return(d1)
}

# obtain_d1_distance <- function(tree1, tree2, data){
#   d1 <- 0
  
#   allcomb <- combinat::combn(nrow(data), 2)
  
#   for (k in 1:dim(allcomb)[2]) { 
#     i <- allcomb[1, k]
#     j <- allcomb[2, k]
    
#     # IT1
#     tree1_predict_node_i <- predict_sample_node(tree1$nodes[[1]], data[i, ], tree1)
#     tree1_predict_node_j <- predict_sample_node(tree1$nodes[[1]], data[j, ], tree1)
    
#     IT1 <- ifelse(tree1_predict_node_i == tree1_predict_node_j, 1, 0)
#     # print(IT1)
    
#     # IT2
#     tree2_predict_node_i <- predict_sample_node(tree2$nodes[[1]], data[i, ], tree2)
#     tree2_predict_node_j <- predict_sample_node(tree2$nodes[[1]], data[j, ], tree2)
    
#     IT2 <- ifelse(tree1_predict_node_i == tree1_predict_node_j, 1, 0)

#     d1 <- d1 + abs(IT1 - IT2)
#   }
  
#   d1 <- d1/choose(nrow(data),2)
#   return(d1)
# }

# obtain_d1_distance_mc <- function(tree1, tree2, data){
#   d1 <- 0
  
#   allcomb <- combinat::combn(nrow(data), 2)
  
#   d1 <- foreach(k = 1:dim(allcomb)[2], .combine = "+") %dopar% {
#     i <- allcomb[1, k]
#     j <- allcomb[2, k]
    
#     # IT1
#     tree1_predict_node_i <- predict_sample_node(tree1$nodes[[1]], data[i, ], tree1)
#     tree1_predict_node_j <- predict_sample_node(tree1$nodes[[1]], data[j, ], tree1)
    
#     IT1 <- ifelse(tree1_predict_node_i == tree1_predict_node_j, 1, 0)
#     # print(IT1)
    
#     # IT2
#     tree2_predict_node_i <- predict_sample_node(tree2$nodes[[1]], data[i, ], tree2)
#     tree2_predict_node_j <- predict_sample_node(tree2$nodes[[1]], data[j, ], tree2)
    
#     IT2 <- ifelse(tree1_predict_node_i == tree1_predict_node_j, 1, 0)
#     abs(IT1 - IT2)
#   }
  
#   print(d1)
#   d1 <- d1/choose(nrow(data),2)
#   return(d1)
# }