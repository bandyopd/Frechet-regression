compute_dcov <- function(n, p, type_data, seed, d_y, r=0.5) {
  set.seed(seed)
  if (type_data == "SPD") {
    m = 3
    #r = 0.5
    spd_data = gendata_SPD(n, p, m, r, seed)
    X <- spd_data$x
    Y <- spd_data$Min
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[, , i]
    }
  } else if (type_data == "density") {
    m <- 20
    #r <- 0.5
    density_data <- gendata_prob(n, p, m, r, seed)
    X <- density_data$x
    Y <- density_data$Qin
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[i, ]
    }
  } else if (type_data == "spherical") {
    #r <- 0.5
    sigma_u <- 0.2
    spherical_data <- gendata_spher(n, p, r, sigma_u, seed)
    X <- spherical_data$x
    Y <- spherical_data$Y
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[i, ]
    }
  } else if (type_data == "cate") {
    m = 3
    #r = 0.5
    spd_data = gendata_SPD_cate(n, p, m, r, seed)
    X <- spd_data$x
    Y <- spd_data$Min
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[, , i]
    }
    
  } else if (type_data == "cate2") {
    m <- 3
    #r <- 0.5
    spd_data <- gendata_SPD_cate2(n, p, m, r, seed)
    X <- spd_data$x
    Y <- spd_data$Min
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[, , i]
    }
  }
  
  dcov_SPD <- rep(0, p)
  if (type_data == "spherical") {
    for (i in 1:p) {
      dcov_SPD[i] <- dcov_eu(X[, i], Y_list, d_y, geo = TRUE)
    }
  } else if (type_data == "cate") {
    for (i in 1:p) {
      x_dummy <- convert2dummy(X[, i], "snp")
      dcov_SPD[i] <- dcov_DTI_cpp(x_dummy, Y_list)
    }
  } else if (type_data == "cate2") {
    for (i in 1:p) {
      if (i == 1) {
        dcov_SPD[i] <- dcov_DTI_cpp(scale(X[, i]), Y_list)
      } else if (i == 2) {
        dcov_SPD[i] <- dcov_DTI_cpp(X[, i], Y_list)
      } else if (i == 3) {
        x_dummy <- convert2dummy(X[, i], "group")
        dcov_SPD[i] <- dcov_DTI_cpp(x_dummy, Y_list)
      } else {
        x_dummy <- convert2dummy(X[, i], "snp")
        dcov_SPD[i] <- dcov_DTI_cpp(x_dummy, Y_list)
      }
    }
    
  } else if (type_data == "SPD") {
    for (i in 1:p) {
      dcov_SPD[i] <- dcov_DTI_cpp(X[, i], Y_list)
    }
  } else if (type_data == "density") {
    for (i in 1:p) {
      dcov_SPD[i] <- dcov_eu(X[, i], Y_list, d_y)
    }
  }
  list(dcov_SPD=dcov_SPD, X=X, Y=Y_list)
  
}




compute_dcov2 <- function(X, Y, type_data, d_y) {
  n <- nrow(X)
  if (type_data == "SPD") {
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[, , i]
    }
  } else if (type_data == "density") {
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[i, ]
    }
  } else if (type_data == "spherical") {
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[i, ]
    }
  } else if (type_data == "cate") {
    
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[, , i]
    }
    
  } else if (type_data == "cate2") {
    
    Y_list <- vector("list", n)
    for (i in 1:n) {
      Y_list[[i]] <- Y[, , i]
    }
  }
  
  p <- ncol(X)
  dcov_SPD <- rep(0, p)
  if (type_data == "spherical") {
    for (i in 1:p) {
      dcov_SPD[i] <- dcov_eu_v2(X[, i], Y_list, d_y, geo = TRUE)
      #dcov_SPD[i] <- sqrt(dcov_eu_v2(X[, i], Y_list, d_y) / sqrt(dcov_one(X[, i], d_Eu) * dcov_one(Y_list, d_y)))
    }
  } else if (type_data == "cate") {
    for (i in 1:p) {
      x_dummy <- convert2dummy(X[, i], "snp")
      dcov_SPD[i] <- dcov_DTI_cpp(x_dummy, Y_list)
    }
  } else if (type_data == "cate2") {
    for (i in 1:p) {
      if (i == 1) {
        dcov_SPD[i] <- dcov_DTI_cpp(scale(X[, i]), Y_list)
      } else if (i == 2) {
        dcov_SPD[i] <- dcov_DTI_cpp(X[, i], Y_list)
      } else if (i == 3) {
        x_dummy <- convert2dummy(X[, i], "group")
        dcov_SPD[i] <- dcov_DTI_cpp(x_dummy, Y_list)
      } else {
        x_dummy <- convert2dummy(X[, i], "snp")
        dcov_SPD[i] <- dcov_DTI_cpp(x_dummy, Y_list)
      }
    }
    
  } else if (type_data == "SPD") {
    for (i in 1:p) {
      dcov_SPD[i] <- dcov_DTI_cpp(X[, i], Y_list)
    }
  } else if (type_data == "density") {
    for (i in 1:p) {
      dcov_SPD[i] <- dcov_eu_v2(X[, i], Y_list, d_y)
      #dcov_SPD[i] <- sqrt(dcov_eu_v2(X[, i], Y_list, d_y) / sqrt(dcov_one(X[, i], d_Eu) * dcov_one(Y_list, d_y)))
    }
  }
  dcov_SPD
  
}