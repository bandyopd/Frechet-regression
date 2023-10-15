# distance covariance for Y in metric space and X in eu
dcov_eu <- function(X, Y, d_y) {
  # X is an n dimensional vector
  # Y is a list with size n
  # d_x and d_y are distance functions in X, and Y
  n <- length(X)
  distance_x <- rep(0, n^2)
  distance_y <- rep(0, n^2)
  V_statistic <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      distance_x[(i-1)*n+j] <- d_Eu(X[i], X[j])
      distance_y[(i-1)*n+j] <- d_y(Y[[i]], Y[[j]])
      for (k in 1:n) {
        V_statistic = d_Eu(X[i], X[j]) * d_y(Y[[i]], Y[[k]]) / n^3 +
          V_statistic
        
      }
    }
  }
  V_statistic <- -2 * V_statistic 
  V_statistic <- sum(distance_x * distance_y) / n^2 + 
    (sum(distance_x) / n^2) * (sum(distance_y) / n^2) + V_statistic
  V_statistic
  
}

# distance covariance of one random object
dcov_one <- function(X, d_x) {
  # X is either a vector or a list
  if (typeof(X) != "list") {
    # if not list, convert to list
    n <- length(X)
    X_list <- vector("list", n)
    for (i in 1:n) {
      X_list[[i]] <- X[i]
    }
  } else {
    X_list <- X
  }
  
  A <- matrix(0, n, n)
  for (j in 1:n) {
    for (k in 1:n) {
      A[j, k] <- d_x(X_list[[j]], X_list[[k]])
    }
  }
  
  # row mean and col mean
  A_row_mean <- apply(A, 1, mean)
  A_col_mean <- apply(A, 2, mean)
  A_mean <- mean(A)
  
  for (j in 1:n) {
    for (k in 1:n) {
      A[j, k] <- A[j, k] - A_row_mean[j] - A_col_mean[k] + A_mean
    }
  }
  sum(A * A) / n^2
}


dcov_eu_v2 <- function(X, Y, d_y, geo=FALSE) {
  n <- length(X)
  A <- matrix(0, n, n)
  B <- matrix(0, n, n)
  if (geo) {
    for (j in 1:n) {
      for (k in 1:n) {
        A[j, k] <- d_Eu(X[j], X[k])
        #print(sum(Y[[j]] * Y[[k]]))
        if (abs(sum(Y[[j]] * Y[[k]])) > 1) {
          B[j, k] <- acos(round(sum(Y[[j]] * Y[[k]])))
        } else {
          B[j, k] <- d_y(Y[[j]], Y[[k]])
        }
      }
    }
  } else {
    for (j in 1:n) {
      for (k in 1:n) {
        A[j, k] <- d_Eu(X[j], X[k])
        #print(sum(Y[[j]] * Y[[k]]))
        B[j, k] <- d_y(Y[[j]], Y[[k]])
      }
    }
  }

  # row mean and col mean
  A_row_mean <- apply(A, 1, mean)
  A_col_mean <- apply(A, 2, mean)
  B_row_mean <- apply(B, 1, mean)
  B_col_mean <- apply(B, 2, mean)
  A_mean <- mean(A)
  B_mean <- mean(B)
  for (j in 1:n) {
    for (k in 1:n) {
      A[j, k] <- A[j, k] - A_row_mean[j] - A_col_mean[k] + A_mean
      B[j, k] <- B[j, k] - B_row_mean[j] - B_col_mean[k] + B_mean
    }
  }
  base::sum(A * B) / n^2
}

d_Eu <- function(x, y) {
  abs(x - y)
}

standardize <- function(x) {
  (x - mean(x)) / sd(x)
}


convert2dummy <- function(x, type) {
  # convert a vector to dummy matrix
  n <- length(x)
  if (type == "snp") {
    n_level <- 3
  } else if (type == "group") {
    n_level <- 5
  } 
  dummy_mat <- matrix(0, n, n_level-1)
  for (i in seq_len(n)) {
    dummy_mat[i, x[i]] <- 1
  }
  dummy_mat
}



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


d_vec <- function(x, y) {
  sqrt(sum((x - y)^2))
}

norm_eu <- function(x) {
  # Euclidean norm of x
  sqrt(sum(x^2))
}

d_geo <- function(x, y) {
  acos(sum(x * y))
}

compute_diff <- function(Y, Y_hat) {
  # both Y and Y_hat are 3*3*m arrays
  m <- dim(Y)[3]
  difference <- vector("double", m)
  for (i in seq_len(m)) {
    difference[i] <- cholesky_error_C(Y[, , i], Y_hat[, , i])^2
  }
  difference
}






# generate data for probability distributions ----------------------------------
# Generates data from example 5.2.2.
gendata_prob=function(n,p,m,r,seed) {
  #' @param n: sample size 
  #' @param p: number of predictors
  #' @param m: number of quantiles for each simulated density output
  #' @param r: correlation between simulated predictors
  
  # Recall the generation of Y=Qin in example 6.1.1: Qin = Vmu + Vsigma%*%t(Ninvpgrid)
  set.seed(seed)
  Qgrid=seq(1/m, 1-1/m, 1/m)
  Ninvpgrid=qnorm(Qgrid)
  x=2*pnorm(matrix(rnorm(n*p),n,p)%*%chol(r^(as.matrix(dist(1:p)))))-1
  
  Vmu0=0
  Vsigma0=3
  Vbeta=3/4
  Vgamma=1
  Vnu1=1
  Vnu2=0.5
  j0=4
  Vmu=Vmu0+Vbeta*(x[,1] + x[, 5] + x[, 9] + x[, 13] + x[,17] + x[,20])+rnorm(n)*sqrt(Vnu1)
  GamP1=(Vsigma0+Vgamma*(x[, 1]))^2/Vnu2
  GamP2=Vnu2/(Vsigma0+Vgamma*(x[, 1]))
  Vsigma=rep(0,n)
  for (j in 1:n) {Vsigma[j]=rgamma(1,shape=GamP1[j], scale=GamP2[j])}
  Qin=Vsigma%*%t(Ninvpgrid)+Vmu
  return(list(x=x, Qin=Qin))
}



# generate data for Spherical data ---------------------------------------------
gendata_spher <- function(n, p, r, tau, seed) {
  #' @param n: sample size 
  #' @param p: number of predictors
  #' @param r: correlation between simulated predictors
  #' @param tau: standard deviation of U
  #' @param seed: seed set for the random generation
  set.seed(seed)
  x=pnorm(matrix(rnorm(n*p),n,p)%*%chol(r^(as.matrix(dist(1:p)))))
  mreg <- cbind(sqrt(1 - x[, 5]^2) * cos(pi*(x[,1] + x[, 5] + x[, 9])),
                sqrt(1 - x[, 5]^2) * sin(pi*(x[,1] + x[, 5] + x[, 9])),
                x[, 5])
  u <- matrix(rnorm(n*2, 0, tau), n, 2)
  # add noise to generate Yi
  Y <- matrix(0, n, 3)
  for (i in 1:n) {
    Y[i, ] <- add_noise(mreg[i, ], u[i, ])
  }
  list(x=x, Y=Y)
}

add_noise <- function(p, u) {
  p1 <- p[1]
  p2 <- p[2]
  p3 <- p[3]
  e11 <- sqrt(p2^2 / (p1^2 + p2^2))
  e1 <- c(e11, -e11 * p1 / p2, 0)
  
  e21 <- 1 / sqrt(1 + p2^2/p1^2 + (p1^2 + p2^2)^2 / (p1^2 * p3^2))
  e2 <- c(e21, e21*p2/p1, -(p1^2 + p2^2) / (p1 * p3)*e21)
  
  v <- u[1] * e1 + u[2] * e2
  cos(norm_eu(v)) * p + sin(norm_eu(v)) * v/norm_eu(v)
  
}

# generate data for SPD data ---------------------------------------------------
# Generates data from example 5.3.2
gendata_SPD=function(n,p,m,r,seed) {
  #' @param n: sample size 
  #' @param p: number of predictors
  #' @param m: dimension of SPD matrix outputs
  #' @param r: correlation between simulated predictors
  #' @param seed: seed set for the random generation
  
  # Recall the generation of Y=Min in example 5.3.2: Min = A^tA
  #                                                  A = (mu + sigma)I + (sigma)U                                                  
  set.seed(seed)
  x=2*pnorm(matrix(rnorm(n*p),n,p)%*%chol(r^(as.matrix(dist(1:p)))))
  #sigma_x <- matrix(r, p, p)
  #diag(sigma_x) <- 1
  #x=2*pnorm(matrix(rnorm(n*p),n,p)%*%chol(sigma_x))
  
  M <- array(0,c(m,m,n))
  Vmu0=3
  Vsigma0=3 
  Vbeta=2 
  Vgamma=3
  
  Vnu1=1
  Vnu2=2  
  
  Vmu=Vbeta*((x[, 1] + x[, 5]))+rnorm(n, Vmu0,sqrt(Vnu1))
  GamP1=(Vsigma0+Vgamma*(x[, 1] + x[, 5] + x[,9] + x[ ,13] + x[, 17] + x[, 20]))^2/Vnu2
  GamP2=Vnu2/(Vsigma0+Vgamma*(x[, 1] + x[, 5] + x[,9] + x[ ,13] + x[, 17] + x[, 20]))
  I =diag(m)
  off_diag = matlab::ones(m)
  off_diag[lower.tri(off_diag, diag = TRUE)] <- 0
  
  for (j in 1:n) {
    set.seed(seed)
    Vsigma=rgamma(1,shape=GamP1[j], scale=GamP2[j]) 
    #A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    aux<-t(A)%*%A
    M[,,j]<-aux
  }
  
  return(list(x=x, Min=M))
}

# generate data for SPD data where x is categorical ----------------------------

# convert norm to categorical
norm2cate <- function(x, seed_freq) {
  # x is a 2n-by-p matrix where each row is normal
  # seed_freq: seed for generate allele_frequency
  n <- dim(x)[1] / 2 
  p <- dim(x)[2]
  x_cate <- matrix(0, n, p)
  set.seed(seed_freq)
  allele_frequency <- runif(p, 0.2, 0.8)
  quantile_C <- qnorm(allele_frequency, lower.tail = FALSE)
  for (i in seq_len(n)) {
    Z1 <- x[2*i-1, ]
    Z2 <- x[2*i, ]
    for (j in seq_len(p)) {
      # if (x[i, j] < -0.26) {
      #   x_cate[i, j] <- 1
      # } else if (x[i, j] > 0.35) {
      #   x_cate[i, j] <- 3
      # } else {
      #   x_cate[i, j] <- 2
      # }
      
      if (Z1[j] > quantile_C[j]) {
        Z1[j] <- 1
      } else {
        Z1[j] <- 0
      }
      if (Z2[j] > quantile_C[j]) {
        Z2[j] <- 1
      } else {
        Z2[j] <- 0
      }
    }
    x_cate[i, ] <- Z1 + Z2
  }
  x_cate
}

# convert a vector to a 5-level cate
vec2cate5 <- function(x) {
  # x: a standard normal vector
  x_mean <- mean(x)
  x_sd <- sd(x)
  quantiles <- qnorm(c(0.25, 0.45, 0.65, 0.85), x_mean, x_sd)
  for (i in seq_along(x)) {
    if (x[i] < quantiles[1]) {
      x[i] <- 0
    } else if (x[i] >= quantiles[1] & x[i] < quantiles[2]) {
      x[i] <- 1
    } else if (x[i] >= quantiles[2] & x[i] < quantiles[3]) {
      x[i] <- 2
    } else if (x[i] >= quantiles[3] & x[i] < quantiles[4]) {
      x[i] <- 3
    } else if (x[i] > quantiles[4]) {
      x[i] <- 4
    }
  } 
  x
}

gendata_SPD_cate <- function(n,p,m,r,seed) {
  #' @param n: sample size 
  #' @param p: number of predictors
  #' @param m: dimension of SPD matrix outputs
  #' @param r: correlation between simulated predictors
  #' @param seed: seed set for the random generation
  
  # Recall the generation of Y=Min in example 5.3.2: Min = A^tA
  #                                                  A = (mu + sigma)I + (sigma)U                                                  
  set.seed(seed)
  x=norm2cate(matrix(rnorm(2*n*p),2*n,p)%*%chol(r^(as.matrix(dist(1:p)))), seed)
  
  M <- array(0,c(m,m,n))
  Vmu0=3
  Vsigma0=3 
  Vbeta=2 
  Vgamma=3
  
  Vnu1=1
  Vnu2=2  
  
  Vmu=Vbeta*((x[, 1] + x[, 5]))+rnorm(n, Vmu0,sqrt(Vnu1))
  GamP1=(Vsigma0+Vgamma*(x[, 1] + x[, 5] + x[,9] + x[ ,13] + x[, 17] + x[, 20]))^2/Vnu2
  GamP2=Vnu2/(Vsigma0+Vgamma*(x[, 1] + x[, 5] + x[,9] + x[ ,13] + x[, 17] + x[, 20]))
  I =diag(m)
  off_diag = matlab::ones(m)
  off_diag[lower.tri(off_diag, diag = TRUE)] <- 0
  
  for (j in 1:n) {
    set.seed(seed)
    Vsigma=rgamma(1,shape=GamP1[j], scale=GamP2[j]) 
    #A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    aux<-t(A)%*%A
    M[,,j]<-aux
  }
  
  return(list(x=x, Min=M))
}

gendata_SPD_cate2 <- function(n,p,m,r,seed) {
  #' @param n: sample size 
  #' @param p: number of predictors
  #' @param m: dimension of SPD matrix outputs
  #' @param r: correlation between simulated predictors
  #' @param seed: seed set for the random generation
  
  # Recall the generation of Y=Min in example 5.3.2: Min = A^tA
  #                                                  A = (mu + sigma)I + (sigma)U                                                  
  set.seed(seed)
  x_snp <- matrix(rnorm(2*n*(p - 3)),2*n)%*%chol(r^(as.matrix(dist(1:(p-3))))) # the original data
  x_snp <- norm2cate(x_snp, seed)
  x_cont <- abs(rnorm(n)) # the continuous variable
  x_binary <- sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  x_cate5 <- vec2cate5(rnorm(n))
  x <- cbind(x_cont, x_binary, x_cate5, x_snp)
  
  M <- array(0,c(m,m,n))
  Vmu0=3
  Vsigma0=3 
  Vbeta=2 
  Vgamma=3
  
  Vnu1=1
  Vnu2=2  
  
  Vmu=Vbeta*((x[, 1] + x[, 5]))+rnorm(n, Vmu0,sqrt(Vnu1))
  GamP1=(Vsigma0+Vgamma*(x[, 1] + x[, 5] + x[,9] + x[ ,13] + x[, 17] + x[, 20]))^2/Vnu2
  GamP2=Vnu2/(Vsigma0+Vgamma*(x[, 1] + x[, 5] + x[,9] + x[ ,13] + x[, 17] + x[, 20]))
  I =diag(m)
  off_diag = matlab::ones(m)
  off_diag[lower.tri(off_diag, diag = TRUE)] <- 0
  
  for (j in 1:n) {
    set.seed(seed)
    Vsigma=rgamma(1,shape=GamP1[j], scale=GamP2[j]) 
    #A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    aux<-t(A)%*%A
    M[,,j]<-aux
  }
  
  return(list(x=x, Min=M))
}

# convert a list to a 3D array, each element in the list is a 3*3 matrix
convert_list23D <- function(L) {
  # L: a list with 3*3 matrices elements
  # output: a 3D array
  n <- length(L)
  L_3D <- array(0, c(3, 3, n))
  for (i in seq_len(n)) {
    L_3D[,, i] <- L[[i]]
  }
  L_3D
}

# compute FA of given matrix A
compute_FA <- function(A) {
  eigen_val <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lambda_1 <- eigen_val[1]
  lambda_2 <- eigen_val[2]
  lambda_3 <- eigen_val[3]
  #lambda_bar <- mean(eigen_val)
  FA <- sqrt(0.5) * sqrt((lambda_1 - lambda_2)^2 + (lambda_2 - lambda_3)^2 +
                           (lambda_3 - lambda_1)^2)
  #FA <- sqrt(1.5) * sqrt((lambda_1 - lambda_bar)^2 + (lambda_2 - lambda_bar)^2 +
  #                         (lambda_3 - lambda_bar)^2)
  FA / sqrt(lambda_1^2 + lambda_2^2 + lambda_3^2)
}


convert_3D2FA <- function(Y) {
  # Y is a 3*3*m array
  # output: a m dimensional FA vector
  m <- dim(Y)[3]
  FA_vec <- rep(0, m)
  for (i in seq_len(m)) {
    FA_vec[i] <- compute_FA(Y[,, i])
  }
  FA_vec
}


convert_list2mat <- function(L) {
  # convert a list L to a matrix
  n <- length(L)
  x <- matrix(nrow = n, ncol = length(L[[1]]))
  for (i in seq_len(n)) {
    x[i, ] <- L[[i]]
  }
  x
}
gen_var_group <- function(p, type_data) {
  # generate variable groups that pass to gglasso function
  # p: number of variables
  if (type_data == "cate") {
    group1 <- rep(1:p, each = 2)
  } else if (type_data == "cate2") {
    group1 <- c(1, 2, rep(3, times = 4), rep(4:p, each = 2))
  }
  group1
}


gen_dummy_X <- function(x, type_data) {
  # generate dummy x matrix
  # x: the original data matrix
  # type_data: "cate" or "cate2"
  n <- nrow(x)
  p <- ncol(x)
  if (type_data == "cate") {
    # "cate": all variables are SNPs
    #        so, x_dummy is n*2p
    x_dummy <- matrix(nrow = n, ncol = 2*p)
    for (i in seq_len(p)) {
      x_dummy[, (2*i-1):(2*i)] <- convert2dummy(x[, i], "snp")
    }
  } else if (type_data == "cate2") {
    # "cate2": x1 is continuous, x2 is binary
    #          x3 is 5-level, other are snp
    x_dummy <- cbind(x[, 1:2], convert2dummy(x[, 3], "group"))
    x_snp <- matrix(nrow = n, ncol = 2*(p-3))
    for (i in seq_len(p-3)) {
      x_snp[, (2*i-1):(2*i)] <- convert2dummy(x[, i+3], "snp")
    }
    x_dummy <- cbind(x_dummy, x_snp)
  }
  x_dummy
}


# compute the cholesky decomposition of a matrix and return the 6 elements of 
# the upper triangular matrix
vec_chol <- function(A) {
  # the 6 elements vector is obtained by stacking rows of the upper tri
  chol_A <- chol_C(A) # upper triangular
  c(chol_A[1, ], chol_A[2, c(2:3)], chol_A[3, 3])
}

# convert 3D array to matrix
convert_3D2mat <- function(Y) {
  # Y: a 3*3*m array
  # output: a 6 * m matrix
  m <- dim(Y)[3]
  Y_mat <- matrix(0, m, 6)
  for (i in seq_len(m)) {
    Y_mat[i, ] <- vec_chol(Y[, , i])
  }
  Y_mat
}

convert_vec2mat <- function(x) {
  # convert the 6 elements of chol decomposition to a matrix
  # x: a 6 elements vector
  # output: an upper tri matrix
  x_mat <- matrix(0, 3, 3)
  x_mat[1, ] <- x[1:3]
  x_mat[2, c(2:3)] <- x[4:5]
  x_mat[3, 3] <- x[6]
  x_mat
}


# convert chol matrix to 3D array
convert_mat23D <- function(Y) {
  # Y: a m*6 matrix, each row is the elements of chol decomposition (upper tri)
  # output: a 3*3*m array
  m <- nrow(Y)
  Y_array <- array(0, c(3, 3, m))
  for (i in seq_len(m)) {
    chol_upper <- convert_vec2mat(Y[i, ])
    Y_array[, , i] <- t(chol_upper) %*% chol_upper
  }
  Y_array
}

# get active predictor, the indices of zero rows of estimated beta
get_ind_row <- function(beta_hat) {
  which(apply(beta_hat, 1, function(x) any(x != 0)))
}

# get active predictor, the indices of nonzero element
get_nonzero_ele <- function(x) {
  # x: a p dimensional coefficient vector
  which(x != 0)
}

# select variable according to the result
# a variable is selected if it is nonzero in one column
select_var_vec <- function(x, group1) {
  # x is a (p+1) dimensional vector, the first element is intercept
  coef_mat <- x[-1] # remove intercept
  unique(group1[get_nonzero_ele(coef_mat)])
}


# select variable according to the result
# a variable is selected if it is nonzero in one column
select_var <- function(x, group1) {
  # x is a (p+1)*6 matrix, the first row is intercept
  coef_mat <- x[-1, ] # remove intercept
  unique(group1[get_ind_row(coef_mat)])
}











