
library(rTensor)
kruskal_product <- function(A, B, C) {
  # Get the dimensions
  R <- ncol(A)
  dim1 <- nrow(A)
  dim2 <- nrow(B)
  dim3 <- nrow(C)
  
  # Initialize the tensor
  tensor <- array(0, dim = c(dim1, dim2, dim3))
  
  # Calculate the Kruskal product
  for (r in 1:R) {
    tensor[,,] <- tensor[,,] + A[, r] %o% B[, r] %o% C[, r]
  }
  
  return(tensor)
}

kruskal_product_high_d <- function(U) {
  # Get the dimensions
  R <- ncol(U[[1]])
  d = length(U)
  dim = unlist(lapply(U, nrow))
  
  # Initialize the tensor
  tensor <- array(0, dim = dim)
  
  # Calculate the Kruskal product
  for (r in 1:R) {
    tmp_tensor = U[[1]][, r]
    for (k in 2:d) {
      tmp_tensor = tmp_tensor %o% U[[k]][, r]
    }
    tensor <- tensor + tmp_tensor
  }
  
  return(tensor)
}







# Function to generate a random CP rank R tensor
generate_random_cp_tensor <- function(dimensions, R) {
  # Create random factor matrices
  A <- matrix(runif(dimensions[1] * R, -1, 1), nrow = dimensions[1], ncol = R)
  B <- matrix(runif(dimensions[2] * R, -1, 1), nrow = dimensions[2], ncol = R)
  C <- matrix(runif(dimensions[3] * R, -1, 1), nrow = dimensions[3], ncol = R)
  
  for (r in 1:R){
    A[,r] = r*A[,r]
  }
  
  # Generate the tensor using the Kruskal product
  tensor <- kruskal_product(A, B, C)
  
  return(list(tensor = tensor, U = list(A = A, B = B, C = C)))
}


generate_random_cp_tensor_condition_number_1_1 <- function(dimensions, R) {
  # Create random factor matrices
  A <- matrix(runif(dimensions[1] * R, 0, 1), nrow = dimensions[1], ncol = R)
  B <- matrix(runif(dimensions[2] * R, 0, 1), nrow = dimensions[2], ncol = R)
  C <- matrix(runif(dimensions[3] * R, 0, 1), nrow = dimensions[3], ncol = R)
  
  A = normalize_columns(A)
  B = normalize_columns(B)
  C = normalize_columns(C)
  
  for (r in 1:R){
    A[,r] = A[,r]
  }
  
  # Generate the tensor using the Kruskal product
  tensor <- kruskal_product(A, B, C)
  
  return(list(tensor = tensor, U = list(A = A, B = B, C = C)))
}


generate_random_cp_tensor_condition_number_1 <- function(dimensions, R) {
  # Create random factor matrices
  A <- matrix(runif(dimensions[1] * R, -1, 1), nrow = dimensions[1], ncol = R)
  B <- matrix(runif(dimensions[2] * R, -1, 1), nrow = dimensions[2], ncol = R)
  C <- matrix(runif(dimensions[3] * R, -1, 1), nrow = dimensions[3], ncol = R)
  
  A = normalize_columns(A)
  B = normalize_columns(B)
  C = normalize_columns(C)
  
  for (r in 1:R){
    A[,r] = A[,r]
  }
  
  # Generate the tensor using the Kruskal product
  tensor <- kruskal_product(A, B, C)
  
  return(list(tensor = tensor, U = list(A = A, B = B, C = C)))
}


generate_random_cp_tensor_high_d <- function(dimensions, R) {
  d <- length(dimensions)
  
  # Create list of random factor matrices
  U <- vector("list", d)
  for (i in 1:d) {
    mat <- matrix(runif(dimensions[i] * R, -1, 1), nrow = dimensions[i], ncol = R)
    #mat = normalize_columns(mat)
    for (r in 1:R) {
      mat[, r] <- r^(1/d) * mat[, r]  # scale each rank-1 component
    }
    U[[i]] <- mat
  }
  
  # Create a KruskalTensor object and convert to full tensor
  ktensor <- kruskal_product_high_d(U)
  
  return(list(tensor = ktensor, U = U))
}

generate_random_cp_tensor_high_d2 <- function(dimensions, R) {
  
  d <- length(dimensions)
  
  # Create list of random factor matrices
  U <- vector("list", d)
  for (i in 1:d) {
    mat <- matrix(runif(dimensions[i] * R, -1, 1), nrow = dimensions[i], ncol = R)
    mat = normalize_columns(mat)
    for (r in 1:R) {
      mat[, r] <- r^(1/d) * mat[, r]  # scale each rank-1 component
    }
    U[[i]] <- mat
  }
  
  # Create a KruskalTensor object and convert to full tensor
  ktensor <- kruskal_product_high_d(U)
  
  return(list(tensor = ktensor, U = U))
}

matrandcong <- function(m, n, gamma) {
  # Create the target Gram matrix (inner products)
  CG <- matrix(gamma, n, n)
  diag(CG) <- 1  # Set diagonal to 1 (norms)
  
  # Cholesky decomposition
  CGR <- chol(CG)
  
  # Generate random m x n matrix and orthonormalize columns (like QR)
  X <- matrix(rnorm(m * n), m, n)

  Q <- svd(X)$u  # Q has orthonormal columns (n columns) Haar measure
  
  
  # Apply transformation to enforce desired inner products
  X_final <- Q %*% CGR
  
  return(X_final)
}


generate_random_cp_tensor_high_d_with_cor <- function(dimensions, R, cor = 0.1) {
  # Create list of random factor matrices
  d= length(dimensions)
  U <- vector("list", d)
  for (i in 1:d) {
    
    mat <- matrandcong(dimensions[i], R,cor)
    for (r in 1:R) {
      mat[, r] <- r^(1/d) * mat[, r]  # scale each rank-1 component
    }
    U[[i]] <- mat
  }
  
  # Create a KruskalTensor object and convert to full tensor
  ktensor <- kruskal_product_high_d(U)
  
  return(list(tensor = ktensor, U = U))
}

# Function to generate a random CP rank 2 tensor with given condtional number
generate_random_cp_tensor_given_condition_number <- function(dimensions, condition_number) {
  R = 2
  # Create random factor matrices
  A <- matrix(runif(dimensions[1] * R, -1,1), nrow = dimensions[1], ncol = R)
  B <- matrix(runif(dimensions[2] * R, -1,1), nrow = dimensions[2], ncol = R)
  C <- matrix(runif(dimensions[3] * R, -1,1), nrow = dimensions[3], ncol = R)
  
  A = normalize_columns(A)
  B = normalize_columns(B)
  C = normalize_columns(C)
  
  A[,1] = condition_number * A[,1]
  
  # Generate the tensor using the Kruskal product
  tensor <- kruskal_product(A, B, C)
  
  return(list(tensor = tensor, U = list(A = A, B = B, C = C)))
}



khatri_rao <- function(A, B) {
  # Initialize the result matrix
  result <- matrix(numeric(0), nrow = nrow(A) * nrow(B), ncol = ncol(A))
  
  # Compute the Khatri-Rao product column-wise
  for (i in 1:ncol(A)) {
    result[, i] <- kronecker(A[, i], B[, i])
  }
  return(result)
}

normalize_columns <- function(A) {
  # Calculate the Euclidean norm of each column
  column_norms <- apply(A, 2, function(col) sqrt(sum(col^2)))
  
  # Divide each column by its norm to normalize
  A_normalized <- sweep(A, 2, column_norms, FUN = "/")
  
  return(A_normalized)
}


my_loss = function(U1, U2){
  R = ncol(U1[[1]])
  column_permutations <- gtools::permutations(R, R)
  lambda_sign = list()
  norm_permu_list = c()
  for (i in 1:nrow(column_permutations)) {
    norm_mode_list = c()
    lambda_sign_tmp = list()
    for (m in 1:length(U1)) {
      A_normalized = normalize_columns(U1[[m]])
      B_normalized = normalize_columns(U2[[m]])
      column_norms1 <- apply(A_normalized[, column_permutations[i, ]] - B_normalized, 2, function(col) sqrt(sum(col^2)))
      column_norms2 <- apply(A_normalized[, column_permutations[i, ]] + B_normalized, 2, function(col) sqrt(sum(col^2)))
      column_norms = data.frame(column_norms1, column_norms2)
      column_norms_which_min = (-apply(column_norms, 1, which.min))*2+3
      column_norms = apply(column_norms, 1, min)
      norm_mode_list = c(norm_mode_list, max(column_norms))
      lambda_sign_tmp[[m]] = column_norms_which_min
    }
    norm_permu_list = c(norm_permu_list, max(norm_mode_list))
    lambda_sign[[i]] = lambda_sign_tmp
  }
  correct_sign = lambda_sign[[which.min(norm_permu_list)]]
  correct_permutation = column_permutations[which.min(norm_permu_list),]
  return(min(norm_permu_list))
}


permutation_and_sign = function(U1, U2){
  R = ncol(U1[[1]])
  column_permutations <- gtools::permutations(R, R)
  lambda_sign = list()
  norm_permu_list = c()
  for (i in 1:nrow(column_permutations)) {
    norm_mode_list = c()
    lambda_sign_tmp = list()
    for (m in 1:length(U1)) {
      A_normalized = normalize_columns(U1[[m]])
      B_normalized = normalize_columns(U2[[m]])
      column_norms1 <- apply(A_normalized[, column_permutations[i, ]] - B_normalized, 2, function(col) sqrt(sum(col^2)))
      column_norms2 <- apply(A_normalized[, column_permutations[i, ]] + B_normalized, 2, function(col) sqrt(sum(col^2)))
      column_norms = data.frame(column_norms1, column_norms2)
      column_norms_which_min = (-apply(column_norms, 1, which.min))*2+3
      column_norms = apply(column_norms, 1, min)
      norm_mode_list = c(norm_mode_list, max(column_norms))
      lambda_sign_tmp[[m]] = column_norms_which_min
    }
    norm_permu_list = c(norm_permu_list, max(norm_mode_list))
    lambda_sign[[i]] = lambda_sign_tmp
  }
  correct_sign = lambda_sign[[which.min(norm_permu_list)]]
  correct_permutation = column_permutations[which.min(norm_permu_list),]
  return(list(correct_sign = correct_sign, correct_permutation = correct_permutation))
}



sim_diag = function(input_tensor, R){
  p = dim(input_tensor)
  input_tensor = as.tensor(input_tensor)
  a = rnorm(p[1])
  a = matrix(a/sqrt(sum(a^2)), ncol = p[1])
  b = rnorm(p[1])
  b = matrix(b/sqrt(sum(b^2)), ncol = p[1])
  
  T1 = ttm(input_tensor, a, 1)@data[1,,]
  T2 = ttm(input_tensor, b, 1)@data[1,,]
  
  U2 = matrix(Re(eigen(T1 %*% MASS::ginv(T2))$vectors[,1:R]), ncol = R)
  U3 = matrix(Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$vectors[,1:R]), ncol = R)
  
  lambda2 = Re(eigen(T1 %*% MASS::ginv(T2))$values[1:R])
  lambda3 = Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$values[1:R])
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    d_v = c(d_v, sum(abs(lambda2[column_permutations[i,]] * lambda3 - 1)))
  }
  U2 = matrix(U2[,column_permutations[which.min(d_v),]], ncol = R)
  
  U1 = matrix(k_unfold(input_tensor, 1)@data %*% t(MASS::ginv(khatri_rao(U3, U2))), ncol = R)
  
  return_list = list(
    U = list(U1,
             U2,
             U3),
    tensor = kruskal_product(U1, U2, U3)
  )
  return(return_list)
}


sim_diag_d4 = function(input_tensor, R){
  p = dim(input_tensor)
  if (length(p) != 4) {stop("d!=4")}
  input_tensor = as.tensor(input_tensor)
  
  #the 3 and 4th mode
  
  a1 = rnorm(p[1])
  a1 = matrix(a1/sqrt(sum(a1^2)), ncol = p[1])
  a2 = rnorm(p[1])
  a2 = matrix(a2/sqrt(sum(a2^2)), ncol = p[1])
  
  b1 = rnorm(p[2])
  b1 = matrix(b1/sqrt(sum(b1^2)), ncol = p[2])
  b2 = rnorm(p[2])
  b2 = matrix(b2/sqrt(sum(b2^2)), ncol = p[2])
  
  T1 = ttl(input_tensor, list(a1,b1), 1:2)@data[1,1,,]
  T2 = ttl(input_tensor, list(a2,b2), 1:2)@data[1,1,,]
  
  U3 = matrix(Re(eigen(T1 %*% MASS::ginv(T2))$vectors[,1:R]), ncol = R)
  U4 = matrix(Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$vectors[,1:R]), ncol = R)
  
  lambda2 = Re(eigen(T1 %*% MASS::ginv(T2))$values[1:R])
  lambda3 = Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$values[1:R])
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    d_v = c(d_v, sum(abs(lambda2[column_permutations[i,]] * lambda3 - 1)))
  }
  U4 = matrix(U4[,column_permutations[which.min(d_v),]], ncol = R)
  
  
  
  #the 1 and 2th mode
  
  a1 = rnorm(p[3])
  a1 = matrix(a1/sqrt(sum(a1^2)), ncol = p[3])
  a2 = rnorm(p[3])
  a2 = matrix(a2/sqrt(sum(a2^2)), ncol = p[3])
  
  b1 = rnorm(p[4])
  b1 = matrix(b1/sqrt(sum(b1^2)), ncol = p[4])
  b2 = rnorm(p[4])
  b2 = matrix(b2/sqrt(sum(b2^2)), ncol = p[4])
  
  T1 = ttl(input_tensor, list(a1,b1), 3:4)@data[,,1,1]
  T2 = ttl(input_tensor, list(a2,b2), 3:4)@data[,,1,1]
  
  U1 = matrix(Re(eigen(T1 %*% MASS::ginv(T2))$vectors[,1:R]), ncol = R)
  U2 = matrix(Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$vectors[,1:R]), ncol = R)
  
  lambda2 = Re(eigen(T1 %*% MASS::ginv(T2))$values[1:R])
  lambda3 = Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$values[1:R])
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    d_v = c(d_v, sum(abs(lambda2[column_permutations[i,]] * lambda3 - 1)))
  }
  U2 = matrix(U2[,column_permutations[which.min(d_v),]], ncol = R)
  
  
  
  # match U1 U3
  
  unfolded_mat_4 <- rs_unfold(input_tensor, m = 4)@data
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    
    U_list = list(U1, U2, as.matrix(U3[,column_permutations[i,]]))
    V <- hadamard_list(lapply(U_list, function(x) {
      t(x) %*% x
    }))
    V_inv <- MASS::ginv(V)
    tmp <- unfolded_mat_4 %*% khatri_rao_list(U_list, 
                                                 reverse = TRUE) %*% V_inv
    U_list[[4]] <- tmp
    Z <- rTensor:::.superdiagonal_tensor(num_modes = 4, 
                                         len = R, elements = 1)
    est <- ttl(Z, U_list, ms = 1:4)
    d_v = c(d_v,fnorm(est - input_tensor))
  }
  
  U3 = matrix(U3[,column_permutations[which.min(d_v),]], ncol = R)
  
  
  U_list = list(U1, U2, U3)
  V <- hadamard_list(lapply(U_list, function(x) {
    t(x) %*% x
  }))
  V_inv <- MASS::ginv(V)
  U4 <- unfolded_mat_4 %*% khatri_rao_list(U_list, 
                                            reverse = TRUE) %*% V_inv
  
  return_list = list(
    U = list(U1,
             U2,
             U3,
             U4),
    tensor = kruskal_product_high_d(list(U1, U2, U3, U4))
  )
  return(return_list)
}








sim_diag_d5 = function(input_tensor, R){
  p = dim(input_tensor)
  if (length(p) != 5) {stop("d!=5")}
  input_tensor = as.tensor(input_tensor)
  
  #the 4 and 5th mode
  
  a1 = rnorm(p[1])
  a1 = matrix(a1/sqrt(sum(a1^2)), ncol = p[1])
  a2 = rnorm(p[1])
  a2 = matrix(a2/sqrt(sum(a2^2)), ncol = p[1])
  
  b1 = rnorm(p[2])
  b1 = matrix(b1/sqrt(sum(b1^2)), ncol = p[2])
  b2 = rnorm(p[2])
  b2 = matrix(b2/sqrt(sum(b2^2)), ncol = p[2])
  
  c1 = rnorm(p[3])
  c1 = matrix(c1/sqrt(sum(c1^2)), ncol = p[3])
  c2 = rnorm(p[3])
  c2 = matrix(c2/sqrt(sum(c2^2)), ncol = p[3])
  
  T1 = ttl(input_tensor, list(a1,b1,c1), 1:3)@data[1,1,1,,]
  T2 = ttl(input_tensor, list(a2,b2,c2), 1:3)@data[1,1,1,,]
  
  U4 = matrix(Re(eigen(T1 %*% MASS::ginv(T2))$vectors[,1:R]), ncol = R)
  U5 = matrix(Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$vectors[,1:R]), ncol = R)
  
  lambda2 = Re(eigen(T1 %*% MASS::ginv(T2))$values[1:R])
  lambda3 = Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$values[1:R])
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    d_v = c(d_v, sum(abs(lambda2[column_permutations[i,]] * lambda3 - 1)))
  }
  U5 = matrix(U5[,column_permutations[which.min(d_v),]], ncol = R)
  
  # my_loss(true_U[4:5], list(U4,U5))
  
  
  
  #the 2 and 3th mode
  
  a1 = rnorm(p[1])
  a1 = matrix(a1/sqrt(sum(a1^2)), ncol = p[1])
  a2 = rnorm(p[1])
  a2 = matrix(a2/sqrt(sum(a2^2)), ncol = p[1])
  
  b1 = rnorm(p[4])
  b1 = matrix(b1/sqrt(sum(b1^2)), ncol = p[4])
  b2 = rnorm(p[4])
  b2 = matrix(b2/sqrt(sum(b2^2)), ncol = p[4])
  
  c1 = rnorm(p[5])
  c1 = matrix(c1/sqrt(sum(c1^2)), ncol = p[5])
  c2 = rnorm(p[5])
  c2 = matrix(c2/sqrt(sum(c2^2)), ncol = p[5])
  
  T1 = ttl(input_tensor, list(a1,b1,c1), c(1,4,5))@data[1,,,1,1]
  T2 = ttl(input_tensor, list(a2,b2,c2), c(1,4,5))@data[1,,,1,1]
  
  U2 = matrix(Re(eigen(T1 %*% MASS::ginv(T2))$vectors[,1:R]), ncol = R)
  U3 = matrix(Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$vectors[,1:R]), ncol = R)
  
  lambda2 = Re(eigen(T1 %*% MASS::ginv(T2))$values[1:R])
  lambda3 = Re(eigen(t(T2) %*% MASS::ginv(t(T1)))$values[1:R])
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    d_v = c(d_v, sum(abs(lambda2[column_permutations[i,]] * lambda3 - 1)))
  }
  U3 = matrix(U3[,column_permutations[which.min(d_v),]], ncol = R)
  # my_loss(true_U[2:3], list(U2,U3))
  
  
  # match U2 U4
  
  unfolded_mat_1 <- rs_unfold(input_tensor, m = 1)@data
  
  column_permutations <- gtools::permutations(R, R)
  d_v = c()
  for (i in 1:nrow(column_permutations)) {
    
    U_list = list(as.matrix(U2[,column_permutations[i,]]), as.matrix(U3[,column_permutations[i,]]),
                  U4,U5)
    V <- hadamard_list(lapply(U_list, function(x) {
      t(x) %*% x
    }))
    V_inv <- MASS::ginv(V)
    tmp <- unfolded_mat_1 %*% khatri_rao_list(U_list, 
                                              reverse = TRUE) %*% V_inv
    U_list = c(list(tmp), U_list)
    Z <- rTensor:::.superdiagonal_tensor(num_modes = 5, 
                                         len = R, elements = 1)
    est <- ttl(Z, U_list, ms = 1:5)
    d_v = c(d_v,fnorm(est - input_tensor))
  }
  
  U2 = matrix(U2[,column_permutations[which.min(d_v),]], ncol = R)
  U3 = matrix(U3[,column_permutations[which.min(d_v),]], ncol = R)
  
  
  U_list = list(U2, U3, U4, U5)
  V <- hadamard_list(lapply(U_list, function(x) {
    t(x) %*% x
  }))
  V_inv <- MASS::ginv(V)
  U1 <- unfolded_mat_1 %*% khatri_rao_list(U_list, 
                                           reverse = TRUE) %*% V_inv
  
  return_list = list(
    U = list(U1,
             U2,
             U3,
             U4,
             U5),
    # my_loss(true_U, U)
    tensor = kruskal_product_high_d(list(U1, U2, U3, U4,U5))
  )
  return(return_list)
}


# dimensions = c(15,14,13)
# dimensions = c(R,R,R,R)
# dimensions = c(R,R,R)
# input_tensor = generate_random_cp_tensor_high_d(dimensions, R)
# true_U = input_tensor$U
# input_tensor = input_tensor$tensor
# my_loss(return_list$U, true_U)


sim_diag_new = function(input_tensor, R){
  p = dim(input_tensor)
  d = length(p)
  input_tensor = as.tensor(input_tensor)
  
  random_list_1 = list()
  for (k in 3:(d)) {
    a1 = rnorm(p[k])
    a1 = matrix(a1/sqrt(sum(a1^2)), ncol = p[k])
    random_list_1[[k-2]] = a1
  }
  
  random_list_2 = list()
  for (k in 3:(d)) {
    a1 = rnorm(p[k])
    a1 = matrix(a1/sqrt(sum(a1^2)), ncol = p[k])
    random_list_2[[k-2]] = a1
  }
  
  T1 = unfold(ttl(input_tensor, random_list_1, 3:d), 1, 2:d)@data
  T2 = unfold(ttl(input_tensor, random_list_2, 3:d), 1, 2:d)@data
  
  U1 = matrix(Re(eigen(T1 %*% MASS::ginv(T2))$vectors[,1:R]), ncol = R)
  
  other_U = MASS::ginv(U1) %*% unfold(input_tensor, 1, 2:d)@data
  
  U_list = list(U1)
  
  for (j in 2:d) {
    U_list[[j]] = matrix(NA, nrow = p[j], ncol = R)
  }
  
  for (r in 1:R) {
    tmp_tensor = drop(fold(other_U[r,, drop = F], 1, 2:d, modes=c(1,p[-1]))@data)
    
    if (d==3) {
      U_list[[2]][,r] = svd(tmp_tensor)$u[,1]
      U_list[[3]][,r] = svd(tmp_tensor)$v[,1]
    }else{
      cp_tmp = cp(as.tensor(tmp_tensor), 1)
      for (j in 2:d) {
        U_list[[j]][,r] = cp_tmp$U[[j-1]]
      }
    }
  }
  
  return_list = list(
    U = U_list
  )
  
  return(return_list)
}


h_d_new = function(input_tensor, R){
  p = dim(input_tensor)
  d = length(p)
  
  if (!is(input_tensor, "Tensor")){
    input_tensor = as.tensor(input_tensor)
  }
  
  result_hosvd = hosvd(input_tensor, rep(R, d))
  
  if (R == 1) {
    U = list()
    for (k in 1:d) {
      U[[k]] = result_hosvd$U[[k]]
    }
  }else{
    result = sim_diag_new(result_hosvd$Z@data, R)
    
    U = list()
    for (k in 1:d) {
      U[[k]] = result_hosvd$U[[k]] %*% result$U[[k]]
    }
  }
  
  
  return_list = list(
    U = U,
    tensor = kruskal_product_high_d(U)
  )
  return(return_list)
}

# dimensions = c(10,11,12,13)
# input_tensor = generate_random_cp_tensor_high_d(dimensions, R)
# true_U = input_tensor$U
# input_tensor = input_tensor$tensor
# my_loss(return_list$U, true_U)



h_d = function(input_tensor, R){
  result_hosvd = hosvd(as.tensor(input_tensor), rep(R, 3))
  result = sim_diag(result_hosvd$Z@data, R)
  
  U1 = result_hosvd$U[[1]] %*% result$U[[1]]
  U2 = result_hosvd$U[[2]] %*% result$U[[2]]
  U3 = result_hosvd$U[[3]] %*% result$U[[3]]
  
  return_list = list(
    U = list(U1,U2,U3),
    tensor = kruskal_product(U1, U2, U3)
  )
  return(return_list)
}



h_d_4 = function(input_tensor, R){
  d = 4
  
  if (!is(input_tensor, "Tensor")){
    input_tensor = as.tensor(input_tensor)
  }

  result_hosvd = hosvd(input_tensor, rep(R, d))
  result = sim_diag_d4(result_hosvd$Z@data, R)
  
  U = list()
  for (k in 1:d) {
    U[[k]] = result_hosvd$U[[k]] %*% result$U[[k]]
  }
  
  return_list = list(
    U = U,
    tensor = kruskal_product_high_d(U)
  )
  return(return_list)
}


h_d_5 = function(input_tensor, R){
  d = 5
  
  if (!is(input_tensor, "Tensor")){
    input_tensor = as.tensor(input_tensor)
  }
  
  result_hosvd = hosvd(input_tensor, rep(R, d))
  result = sim_diag_d5(result_hosvd$Z@data, R)
  
  U = list()
  for (k in 1:d) {
    U[[k]] = result_hosvd$U[[k]] %*% result$U[[k]]
  }
  
  return_list = list(
    U = U,
    # my_loss(true_U, U)
    tensor = kruskal_product_high_d(U)
  )
  return(return_list)
}




my_cp = function (tnsr,
                  initilization, 
                  num_components = NULL, max_iter = 25, tol = 1e-05,
                  test_conv = F) 
{
  #the initilization should be a list of matrices, each one represents singularvectors for one mode. 
  
  if (is.null(num_components)) 
    stop("num_components must be specified")
  stopifnot(is(tnsr, "Tensor"))
  num_modes <- tnsr@num_modes
  modes <- tnsr@modes
  U_list <- initilization
  unfolded_mat <- vector("list", num_modes)
  tnsr_norm <- fnorm(tnsr)
  for (m in 1:num_modes) {
    unfolded_mat[[m]] <- rs_unfold(tnsr, m = m)@data
    #U_list[[m]] <- matrix(rnorm(modes[m] * num_components), 
    #                      nrow = modes[m], ncol = num_components)
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est) {
    curr_resid <- fnorm(est - tnsr)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter == 1) 
      return(FALSE)
    if (abs(curr_resid - fnorm_resid[curr_iter - 1])/tnsr_norm < 
        tol) 
      return(TRUE)
    else {
      return(FALSE)
    }
  }
  pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  norm_vec <- function(vec) {
    norm(as.matrix(vec), type = "F")
  }
  
  U_all_list = list()
  
  while ((curr_iter < max_iter) && (!converged)) {
    setTxtProgressBar(pb, curr_iter)
    for (m in 1:num_modes) {
      V <- hadamard_list(lapply(U_list[-m], function(x) {
        t(x) %*% x
      }))
      V_inv <- MASS::ginv(V)
      tmp <- unfolded_mat[[m]] %*% khatri_rao_list(U_list[-m], 
                                                   reverse = TRUE) %*% V_inv
      lambdas <- apply(tmp, 2, norm_vec)
      U_list[[m]] <- sweep(tmp, 2, lambdas, "/")
      Z <- rTensor:::.superdiagonal_tensor(num_modes = num_modes, 
                                           len = num_components, elements = lambdas)
      est <- ttl(Z, U_list, ms = 1:num_modes)
    }
    
    
    if (test_conv){
      U_all_list[[curr_iter]] = U_list
    }
    
    if (CHECK_CONV(est)) {
      converged <- TRUE
      setTxtProgressBar(pb, max_iter)
    }
    else {
      curr_iter <- curr_iter + 1
    }
  }
  if (!converged) {
    setTxtProgressBar(pb, max_iter)
  }
  close(pb)
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1)/tnsr_norm)) * 
    100
  invisible(list(lambdas = lambdas, U = U_list, conv = converged, 
                 est = est, norm_percent = norm_percent, fnorm_resid = tail(fnorm_resid, 
                                                                            1), all_resids = fnorm_resid,
                 U_all_list = U_all_list))
}


our_method = function(input_tensor, R){
  initialization_result = h_d(input_tensor, R)
  tnsr= as.tensor(input_tensor)
  initilization = initialization_result$U
  cp_result = my_cp(tnsr,initilization, R)
  return(cp_result)
}



our_method_d_3 = function(input_tensor, R, tol = 0, max_iter = 10){
  initialization_result = h_d(input_tensor, R)
  tnsr= as.tensor(input_tensor)
  initilization = initialization_result$U
  cp_result = my_cp(tnsr,initilization, R, tol = tol,max_iter = max_iter)
  return(list(cp_result = cp_result,
              initialization_loss_X = 
                fnorm(as.tensor(initialization_result$tensor - input_tensor))))
}

our_method_d_4 = function(input_tensor, R, tol = 0, max_iter = 10){
  initialization_result = h_d_4(input_tensor, R)
  tnsr= as.tensor(input_tensor)
  initilization = initialization_result$U
  cp_result = my_cp(tnsr,initilization, R, tol = tol,max_iter = max_iter)
  return(list(cp_result = cp_result,
              initialization_loss_X = 
                fnorm(as.tensor(initialization_result$tensor - input_tensor))))
}



our_method_d_5 = function(input_tensor, R, tol = 0, max_iter = 10){
  initialization_result = h_d_5(input_tensor, R)
  tnsr= as.tensor(input_tensor)
  initilization = initialization_result$U
  cp_result = my_cp(tnsr,initilization, R, tol = tol, max_iter = max_iter)
  return(list(cp_result = cp_result,
              initialization_loss_X = 
                fnorm(as.tensor(initialization_result$tensor - input_tensor))))
}

our_method_new = function(input_tensor, R, tol = 0, max_iter = 10, test_conv = F){
  initialization_result = h_d_new(input_tensor, R)
  tnsr= as.tensor(input_tensor)
  initilization = initialization_result$U
  
  cp_result = my_cp(tnsr,initilization, R, tol = tol, max_iter = max_iter, test_conv = test_conv)
  
  if (test_conv){
    U_all_list = c(list(initialization_result$U), cp_result$U_all_list)
  }else{
    U_all_list = list()
  }
  
  return(list(cp_result = cp_result,
              initialization_loss_X = 
                fnorm(as.tensor(initialization_result$tensor - input_tensor)),
              U_all_list = U_all_list))
}

# dimensions = c(10,11,12,13)
# input_tensor = generate_random_cp_tensor_high_d(dimensions, R)
# true_U = input_tensor$U
# input_tensor = input_tensor$tensor
# my_loss(cp_result$U, true_U)



r1als = function(input_tensor, R, L = 50, nu = 0.1){
  # random_tensor = generate_random_cp_tensor_given_condition_number(dimensions, 2)
  # random_tensor = generate_random_cp_tensor(dimensions, R = 2)
  # random_tensor = generate_random_cp_tensor_condition_number_1(dimensions, R = 2)
  # input_tensor = random_tensor[["tensor"]]
  # input_tensor = random_tensor$tensor + noise_tensor@data * sigma
  if (!is(input_tensor, "Tensor")){
    input_tensor = as.tensor(input_tensor)
  }
  
  R1_results = list()
  all_T_abc_values = c()
  for (i in 1:L) {
    re_tmp = cp(input_tensor, 1)
    R1_results[[i]] = list(lambda = re_tmp$lambdas, U = re_tmp$U)
    re_tmp$U[[1]] = t(re_tmp$U[[1]])
    re_tmp$U[[2]] = t(re_tmp$U[[2]])
    re_tmp$U[[3]] = t(re_tmp$U[[3]])
    all_T_abc_values = c(all_T_abc_values, as.vector(ttl(input_tensor, re_tmp$U,ms= 1:3)@data))
  }
  
  all_T_abc_values = abs(all_T_abc_values)
  
  U1 = matrix(nrow = dim(input_tensor)[1], ncol = R)
  U2 = matrix(nrow = dim(input_tensor)[2], ncol = R)
  U3 = matrix(nrow = dim(input_tensor)[3], ncol = R)
  
  all_T_abc_values_remaining = data.frame(index = 1:L, values = all_T_abc_values)
  checka = 0
  for (r in 1:R) {
    r_index = all_T_abc_values_remaining$index[which.max(all_T_abc_values_remaining$values)]
    
    if (length(r_index) == 0){
      checka = 1
      break
    }
    
    U1[,r] = R1_results[[r_index]]$U[[1]]
    U2[,r] = R1_results[[r_index]]$U[[2]]
    U3[,r] = R1_results[[r_index]]$U[[3]]
    
    if (r == R){
      break
    }
    
    all_T_abc_values_remaining = all_T_abc_values_remaining[-which.max(all_T_abc_values_remaining$values),]
    
    value_table = data.frame(index = all_T_abc_values_remaining$index, 
                             inner_values = NA)
    for (i in 1:nrow(value_table)) {
      u = c(sum(U1[,r] * R1_results[[value_table$index[i]]]$U[[1]]),
            sum(U2[,r] * R1_results[[value_table$index[i]]]$U[[2]]),
            sum(U3[,r] * R1_results[[value_table$index[i]]]$U[[3]]))
      value_table$inner_values[i] = max(abs(u))
    }
    index_to_remove = value_table$index[value_table$inner_values > nu]
    all_T_abc_values_remaining = all_T_abc_values_remaining[!(all_T_abc_values_remaining$index %in% index_to_remove), ]
  }
  U = list(U1, U2, U3)
  if (checka == 1){
    U=NA
  }
  return(list(U = U))
}



# my_loss(list(random_tensor$U[[1]][,2,drop = FALSE]), list(U1[,1,drop = FALSE]))


#— Helpers for both algos —#

# mats(input_tensor, S): step 2 of Alg 3
#   unfold input_tensor into a prod(d[S])×prod(d[-S]) matrix
mats <- function(input_tensor, S) {
  dims <- dim(input_tensor) 
  dS   <- prod(dims[S])
  dO   <- prod(dims[-S])
  perm <- c(S, setdiff(seq_along(dims), S))
  Tp   <- aperm(input_tensor, perm)
  matrix(Tp, nrow = dS, ncol = dO)
}

# mat_k(x_vec, sub_dims, k_sub): step 4 of Alg 3
#   reshape x_vec of length prod(sub_dims) into array dim=sub_dims
#   then mode-k_sub-unfold it into sub_dims[k_sub]×prod(sub_dims[-k_sub])
mat_k <- function(x_vec, sub_dims, k_sub) {
  A    <- array(x_vec, dim = sub_dims)
  perm <- c(k_sub, setdiff(seq_along(sub_dims), k_sub))
  M    <- aperm(A, perm)
  matrix(M, nrow = sub_dims[k_sub], ncol = prod(sub_dims[-k_sub]))
}

# choose S to maximize min(d_S, d/d_S): step 1 of Alg 3
choose_S <- function(dims) {
  N      <- length(dims)
  allS   <- unlist(lapply(1:floor(N/2),
                          function(s) combn(N, s, simplify=FALSE)),
                   recursive = FALSE)
  best   <- NULL
  bestV  <- -Inf
  for (S in allS) {
    v <- min(prod(dims[S]), prod(dims[-S]))
    if (v > bestV) { bestV <- v; best <- S }
  }
  best
}

# right-inverse of a d×r matrix A: needed in both algos
right_inverse <- function(A) MASS::ginv(A)


#— Algorithm 3: Composite PCA (CPCA) —#
cpca <- function(input_tensor, R, S = integer(0)) {
  # random_tensor = generate_random_cp_tensor_given_condition_number(dimensions, 2)
  # input_tensor = random_tensor[["tensor"]]
  dims <- dim(input_tensor); N <- length(dims)
  # 1. pick S if not supplied
  if (length(S)==0) S <- choose_S(dims)
  
  # 2. unfold
  M    <- mats(input_tensor, S)             # d_S × (d/d_S)
  
  # 3. top-r SVD
  sv   <- svd(M, nu = R, nv = R)
  lambda_hat <- sv$d[1:R]
  Uhat <- sv$u[,1:R,drop=FALSE]
  Vhat <- sv$v[,1:R,drop=FALSE]
  
  # 4. recover each mode-k factor â_jk
  Aout <- vector("list", N)
  for (k in seq_len(N)) {
    Aout[[k]] <- matrix(0, nrow = dims[k], ncol = R)
    for (j in seq_len(R)) {
      if (k %in% S) {
        pos   <- which(S == k)
        Mkj   <- mat_k(Uhat[,j], dims[S], pos)
      } else {
        Sc    <- setdiff(seq_len(N), S)
        pos   <- which(Sc == k)
        Mkj   <- mat_k(Vhat[,j], dims[Sc], pos)
      }
      Aout[[k]][,j] <- svd(Mkj, nu=1)$u
    }
  }
  
  list(U = Aout, lambda = lambda_hat)
}

# random_tensor = generate_random_cp_tensor_given_condition_number(dimensions, 1)
# random_tensor = generate_random_cp_tensor_given_condition_number(dimensions, 2)
# my_loss(cpca(random_tensor$tensor, R = 2)$U, random_tensor$U)


#— Algorithm 4: Iterative Concurrent Orthogonalization (ICO) —#
ico <- function(input_tensor, R, a0, tol = 1e-5, max_iter = 20) {
  #input_tensor = random_tensor$tensor
  #a0 = cpca(random_tensor$tensor, R = 2)$U
  dims <- dim(input_tensor); N <- length(dims)
  A_prev <- a0
  B_prev <- lapply(A_prev, right_inverse)  # step 1 → {b̂⁽¹⁾_jk}
  m <- 0
  
  repeat {
    m <- m + 1
    A_new <- vector("list", N)
    B_new <- vector("list", N)
    
    for (k in seq_len(N)) {
      # step 4: form mode-k unfolding once
      M_k <- mats(input_tensor, k)  # dims[k] × prod(dims[-k])
      A_new[[k]] <- matrix(0, nrow = dims[k], ncol = R)
      
      for (j in seq_len(R)) {
        # step 6: contract input_tensor along all modes ≠ k with b̂⁽m⁾_jll
        others <- setdiff(seq_len(N), k)
        bvecs  <- lapply(others, function(ll) B_prev[[ll]][j, ])
        # outer product across all ll∈others
        b_arr  <- Reduce(outer, bvecs)
        b_flat <- as.vector(b_arr)
        
        v       <- M_k %*% b_flat  # ∈ ℝ^{d_k}
        A_new[[k]][, j] <- v / sqrt(sum(v^2))  # step 7
      }
      
      # step 9: update right-inverse for mode k
      B_new[[k]] <- right_inverse(A_new[[k]])
    }
    
    # step 13: check convergence of projection matrices
    max_diff <- 0
    for (k in seq_len(N)) {
      for (j in seq_len(R)) {
        P_old <- A_prev[[k]][,j] %*% t(A_prev[[k]][,j, drop = F])
        P_new <- A_new [[k]][,j] %*% t(A_new [[k]][,j, drop = F])
        d_k   <- norm(P_new - P_old, type = "2")
        max_diff <- max(max_diff, d_k)
      }
    }
    A_prev <- A_new
    B_prev <- B_new
    if (max_diff <= tol || m >= max_iter) break
  }
  
  # step 12: final λ̂⁽m⁾_j = |T ×ₖ b̂⁽m⁾_jk|
  lambda <- numeric(R)
  for (j in seq_len(R)) {
    bvecs  <- lapply(seq_len(N), function(ll) B_prev[[ll]][j, ])
    b_arr  <- Reduce(outer, bvecs)
    lambda[j] <- abs(sum(input_tensor * b_arr))
  }
  
  list(U = A_new, lambda = lambda)
  
  # my_loss(A_new, random_tensor$U)
}



r1als_2 = function(input_tensor, R, max_iter  = 20, tol  = 10^(-5)){
  # random_tensor = generate_random_cp_tensor_given_condition_number(dimensions, 2)
  # input_tensor = random_tensor[["tensor"]]
  if (!is(input_tensor, "Tensor")){
    input_tensor = as.tensor(input_tensor)
  }
  dims <- dim(input_tensor)
  components <- list(U1 = matrix(0, nrow = dims[1], ncol = R),
                     U2 = matrix(0, nrow = dims[2], ncol = R),
                     U3 = matrix(0, nrow = dims[3], ncol = R))
  
  residual <- input_tensor
  
  for (r in 1:R) {
    # HOSVD initialization
    hosvd <- hosvd(residual, ranks = dims)
    a <- hosvd$U[[1]][, 1, drop = F]
    b <- hosvd$U[[2]][, 1, drop = F]
    c <- hosvd$U[[3]][, 1, drop = F]
    
    # Power iteration
    result = my_cp(residual, initilization = list(a,b,c), num_components = 1)
    
    components$U1[, r] <- a
    components$U2[, r] <- b
    components$U3[, r] <- c
    
    # Deflate residual
    outer_prod <- result$est
    residual <- residual - outer_prod
  }
  
  return(list(U = components))
}





# random_tensor = generate_random_cp_tensor_given_condition_number(dimensions, 2)
# input_tensor = random_tensor[["tensor"]]






