# ultility functions to compute smoothed prediction


# ------ local average ---------------------
# the weight for voxel i is 1 and 1/2 for it's neighbors

compute_local_avg <- function(Y, weight) {
  # weight is a vector containing weights for each Y_i, 0 for non-exist nbrs
  # Y is a 3-by-3-by-n array
  nbrs <- which(weight != 0) # nbrs with weight 0 means they are not exist
  Y_nbrs <- Y[, , nbrs]
  #print(dim(Y_nbrs))
  #print(weight)
  if (length(dim(Y_nbrs)) != 3) { # has no nbr
    return(Y_nbrs)
  }
  weight_nbr <- weight[nbrs]
  # get chol Y
  n <- length(nbrs)
  chol_Y <- vector("list", n)
  for (i in seq_len(n)) {
    chol_Y[[i]] <- chol_C(Y_nbrs[, , i])
  }
  chol_Y <- matrix(unlist(chol_Y), 3)
  WeightedFrechetMeanCpp(chol_Y, weight_nbr, 3)
}



gen_nbr <- function(n_vox, weight = 0.5) {
  # assume the voxel is on a 2D grid
  # n_vox: number of voxels
  grid_size <- sqrt(n_vox)
  nbr_ind <- vector("list", n_vox)
  weight_l <- vector("list", n_vox)
  g_full <- graph.lattice(length = grid_size, dim = 2)
  net <- as.matrix(get.adjacency(g_full, attr=NULL))
  for (i in seq_len(n_vox)) {
    nbr_ind[[i]] <- c(i, which(net[i, ] != 0)) # including itself
    weight_l[[i]] <- c(1, rep(weight, length(which(net[i, ] != 0))))
  }
  list(nbr_ind = nbr_ind, weight_l = weight_l)
}


# computed smoothed Frechet
compute_smooth_Fr <- function(pred_Fr, weight = 0.5, 
                              seed_ = 2022) {
  # pred_Fr: a 3*3*n*nrep array, containing predicted SPD
  nrep <- dim(pred_Fr)[4]   # number of repetitions
  n <- dim(pred_Fr)[3] # number of samples
  
  nbr_weight <- gen_nbr(nrep, weight = weight)
  nbr_ind <- nbr_weight$nbr_ind
  weight_l <- nbr_weight$weight_l
  
  # call compute_local_avg to get weighted Y
  Y_weighted <- array(dim = c(3, 3, n, nrep)) 
  for (i in seq_len(n)) {
    Y_tmp <- pred_Fr[, , i, ]
    for (j in seq_len(nrep)) {
      Y_weighted[, , i, j] <- compute_local_avg(Y_tmp[,,nbr_ind[[j]]], 
                                                weight_l[[j]])
    }
  }
  Y_weighted
}


gen_SWP <- function(x, correlation, grid_size, m, seed_) {
  # generate spatially correlated SPD matrices using spatial Wishart process
  # described in the paper "Geostatistical Modeling of Positive Definite 
  # Matrices and Its Applications to Diffusion Tensor Imaging".
  # x: a n*p design matrix, n is the number of subject, p is the number of variables
  # correlation: correlation parameter that control the strength of spatial
  #              correlation
  # grid_size: a grid_size*grid_size grid is generated 
  # m: degree of freedom of the Wishart distribution of U
  # seed_: seed for reproducing results
  set.seed(seed_)
  seqx <- seq_len(grid_size)
  seqy <- seq_len(grid_size)
  s <- expand.grid(x = seqx, y = seqy) # expand.grid in package 'base'
  N <- grid_size^2 # number of SPD matrices on the grid, number of voxels
  cov <- varcov.spatial(
    coords = s,
    cov.model = "exponential",
    cov.pars = c(1, correlation)
  )$`varcov`
  
  beta11 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta22 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta33 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta21 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta31 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta32 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  n <- nrow(x) # number of row
  Y <- array(dim = c(3, 3, n, N)) # 3*3*n sample*n voxel
  for (sub in seq_len(n)) {
    # for each subject
    D1=rmvn(n = m, mu=rep(0, N), cov)
    D2=rmvn(n = m, mu=rep(0, N), cov)
    D3=rmvn(n = m, mu=rep(0, N), cov)
    
    for (v in seq_len(N)) {
      # for each voxel
      L <- matrix(0, 3, 3)
      L[1, 1] <- exp(x[sub, c(1, 5, 10)]%*%beta11[, v])
      L[2, 2] <- exp(x[sub, c(1, 5, 10)]%*%beta22[, v])
      L[3, 3] <- exp(x[sub, c(1, 5, 10)]%*%beta33[, v])
      L[2, 1] <- x[sub, c(15, 20, 25)]%*%beta33[, v]
      L[3, 1] <- x[sub, c(15, 20, 25)]%*%beta33[, v]
      L[3, 2] <- x[sub, c(15, 20, 25)]%*%beta33[, v]
      
      mm <- cbind(D1[, v], D2[, v], D3[, v])
      Y[, , sub, v] <- L%*%t(mm)%*%mm%*%t(L)/m
    }
  }
  Y
}


