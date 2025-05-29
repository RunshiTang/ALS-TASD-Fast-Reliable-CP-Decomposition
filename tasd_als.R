library(rTensor)


# input: a multi array (NOT a rTensor::tensor object)
# output: result in the format of rTensor::cp()
# part of this program is adopted from rTensor::cp()

tasd_als = function(input_tensor, R, tol = 0, max_iter = 10, test_conv = F){
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




# helper functions


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








