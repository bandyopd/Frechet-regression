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

# map the coordinate (x, y, z) to the index (from 1 to number of voxel)
# lapply(c(1:22), function(x) apply(ROI_xyz[[x]], 2, min))
map_tri_idx_2_single_idx <- function(ROI_ind) {
  # ROI_ind: an integer, from 1 to 22
  # output: a 3D array, element is zero if no voxel at that location
  indices <- ROI_xyz[[ROI_ind]]
  n_voxel <- nrow(indices)
  # make sure x+1, y+1, z+1 still valid
  # for most regions, min(x), min(y) and min(z) are greater than 1
  # thus, x-1, y-1 and z-1 are valid
  tri_idx_2_single_idx <-
    array(0, c(max(indices[, 1])+1, max(indices[, 2])+1, max(indices[, 3])+1)) 
  for (i in seq_len(n_voxel)) {
    trans_idx <- indices[i, ]
    tri_idx_2_single_idx[trans_idx[1], trans_idx[2], trans_idx[3]] <- i
  }
  tri_idx_2_single_idx
}

# indices <- cbind(pull(DTI[[1]], "x"), pull(DTI[[1]], "y"), pull(DTI[[1]], "z"))
# x_range <- range(pull(DTI[[1]], "x")) # 24 - 67 
# y_range <- range(pull(DTI[[1]], "y")) # 20 - 87 
# z_range <- range(pull(DTI[[1]], "z")) # 30 - 57 
# tri_idx_2_single_idx <- array(0, c(68, 88, 58)) # make sure x+1, y+1, z+1 still valid
# #base_idx <- c(23, 19, 29) # subtract this such that idx goes from (1, 1, 1)
# for (i in 1:15273) {
#   trans_idx <- indices[i, ]
#   tri_idx_2_single_idx[trans_idx[1], trans_idx[2], trans_idx[3]] <- i
# }


get_neighbors <-
  function(x, y, z, indices, Y_scan,
           tri_idx_2_single_idx,
           weight_nbr = 0.5) {
    # indices: the indices of all voxel in the form of (x, y, z)
    # Y_scan: 3-by-3-by-n_vox array, where n_vox is the number of voxels
    # return a 3-by-3-by-n array of neighbors of voxel at (x, y, z)
    Y_nbrs <- array(0, c(3, 3, 27)) # 27 = 3^3 number of nbrs including itself
    weights <- rep(0, 27)
    cnt <- 1
    
    for (i in c(x - 1, x, x + 1)) {
      for (j in c(y - 1, y, y + 1)) {
        for (k in c(z - 1, z, z + 1)) {
          idx <- tri_idx_2_single_idx[i, j, k]
          if (idx != 0) {
            Y_nbrs[, , cnt] <- Y_scan[, , idx]
            if (i == x & j == y & k == z) {
              weights[cnt] <- 1
            } else {
              weights[cnt] <- weight_nbr
            }
          } else {
            Y_nbrs[, , cnt] <- matrix(0, 3, 3) # or identity matrix?
          }
          cnt = cnt + 1
          
        }
      }
    }
    #print(cnt)
    list(Y_nbrs, weights, length(which(weights != 0)))
  }



get_Y_pred <- function(idx, all_Y_pred) {
  # extracts predicted Y of all voxels for subject idx 
  # idx: a number from 1 to 169
  # all_Y_pred: a list with length n_voxel, each element is 3*3*n_sample matrix
 
  n_vox <- length(all_Y_pred)
  Y_pred <- array(0, c(3, 3, n_vox))
  for (i in seq_len(n_vox)) {
    Y_pred[, , i] <- all_Y_pred[[i]][, , idx]
  }
  Y_pred
}


compute_avg_Y <- function(idx, indices, tri_idx_2_single_idx, all_Y_pred) {
  # compute the locally weighted frechet mean for subject idx
  # idx: a number from 1 to 169
  # indices: nvox-by-3 matrix, rows are the 3-D coordinates
  # tri_idx_2_single_idx
  # all_Y_pred: a list with length n_voxel, each element is 3*3*n_sample matrix
  Y_pred <- get_Y_pred(idx, all_Y_pred)
  n <- dim(Y_pred)[3] # number of voxels
  Y_local_avg <- array(0, c(3, 3, n))
  number_nbrs <- rep(0, n)
  for (i in seq_len(n)) {
    xyz <- indices[i, ]
    nbrs <- get_neighbors(xyz[1], xyz[2], xyz[3],
                          indices, Y_pred, tri_idx_2_single_idx)
    Y_local_avg[, , i] <- compute_local_avg(nbrs[[1]], nbrs[[2]])
    number_nbrs[i] <- nbrs[[3]]
  }
  list(Y_local_avg, number_nbrs)
}



