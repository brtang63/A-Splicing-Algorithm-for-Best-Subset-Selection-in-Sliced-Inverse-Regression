# implement sir_splicing -----------------------------------------------------
get_top_inds <- function(vec, k) order(vec, decreasing = TRUE)[1:k]

get_bottom_inds <- function(vec, k) order(vec, decreasing = FALSE)[1:k]

convert_slice <- function(vec, H) {
  breaks <- quantile(vec, probs = seq(0, 1, length.out = H + 1))
  discretized <- cut(vec, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  return(discretized)
}

sgevd <- function(M, Sigma, A, d) {
  M_A <- M[A, A]
  Sigma_A <- Sigma[A, A]
  
  s_A <- nrow(M_A)
  
  if (is.invertible(Sigma_A)) {
    eigen <- geigen_top(M_A, Sigma_A, d)
  } else {
    eigen <- geigen_top(M_A, Sigma_A + 1e-8 * diag(s_A), d)
  }
  
  if (d == 2) {
    Lambda <- diag(eigen$values)
  }
  
  if (d == 1) {
    Lambda <- eigen$values
  }
  B <- eigen$vectors
  
  return(list(Lambda = Lambda, B = B))
}

is.invertible <- function(matrix) {
  return(det(matrix) != 0)
}

geigen_top <- function(M, Sigma, d) {
  res_eig <- geigen::geigen(M, Sigma, symmetric=T)
  ind <- order(res_eig$values, decreasing = TRUE)[1:d]
  res_eig$values <- res_eig$values[ind]
  res_eig$vectors <- as.matrix(res_eig$vectors[,ind])
  
  return(res_eig)
}


sir_splicing_once <- function(M, Sigma, s, d, H = 5, A_init) {
  
  p <- ncol(M)
  S <- 1:p
  A <- A_init
  I <- setdiff(S, A) 
  
  obj_sgevd <- sgevd(M, Sigma, A, d)
  Lambda <- obj_sgevd$Lambda
  B_A <- obj_sgevd$B
  B <- matrix(0, p, d)
  B[A, ] <- B_A
  
  D_I <- (-2 * M %*% Sigma %*% B + B %*% Lambda)[I, ]
  D <- matrix(0, p, d)
  D[I, ] <- D_I
  
  L <- sum(diag(t(B) %*% M %*% B))
  
  mmax <- 50
  tmax <- min(5, s)
  ite <- 0
  for (m in 1:mmax) {
    ite <- ite + 1
    A_old <- A
    I_old <- I
    for (t in 1:tmax) {
      bi <- rowSums(B[A, , drop = F]^2)
      di <- rowSums(D[I, , drop = F]^2)
      
      S1 <- A[get_bottom_inds(bi, t)]
      S2 <- I[get_top_inds(di, t)]
      
      A_new <- union(setdiff(A, S1), S2)
      I_new <- setdiff(S, A_new)
      
      obj_sgevd <- sgevd(M, Sigma, A_new, d)
      Lambda_new <- obj_sgevd$Lambda
      B_A_new <- obj_sgevd$B
      B_new <- matrix(0, p, d)
      
      B_new[A_new, ] <- B_A_new
      
      D_I_new <- (-2 * M %*% B_new + 2 * Sigma %*% B_new %*% Lambda_new)[I_new, ]
      D_new <- matrix(0, p, d)
      D_new[I_new, ] <- D_I_new
      
      L_new <- sum(diag(t(B_new) %*% M %*% B_new))
      if (L_new > L) {
        L <- L_new
        A <- A_new
        I <- I_new
        B <- B_new
        Lambda <- Lambda_new
      }
    }
    if (setequal(A, A_old)) break
  }
  return(list(A = A, B = B, L = L, s = s))
}

sir_splicing_cv <- function(x, y, s_seq, d, K, H = 5, categorical = FALSE) {
  n <- nrow(x)
  
  if (categorical) {
    folds <- cut(sample(n), breaks=K, labels=FALSE)
  } else {
    folds <- cut(seq(1, n), breaks=K, labels=FALSE)
  }
  
  # Store results
  results <- list()
  avg_loss <- numeric(length(s_seq))
  
  # Perform K-fold cross-validation
  for(i in 1:K){
    # Split data into training and validation sets
    test_indices <- which(folds == i, arr.ind = TRUE)
    train_indices <- setdiff(seq_len(n), test_indices)
    
    x_train <- x[train_indices, ]
    y_train <- y[train_indices]
    
    x_test <- x[test_indices, ]
    y_test <- y[test_indices]
    
    # initialize
    p <- ncol(x_train)
    
    if (categorical) {
      y_trans <- y_train 
      y_unique <- as.numeric(names(table(y_trans)))
      H <- length(y_unique)
      
      phat <- sapply(y_unique, function(h) mean(y_trans == h))
      Ehat <- sapply(y_unique, function(h) colMeans(x_train[y_trans == h,]) - colMeans(x_train))
    } else {
      y_trans <- convert_slice(y_train, H)
      
      phat <- sapply(1:H, function(h) mean(y_trans == h))
      Ehat <- sapply(1:H, function(h) colMeans(x_train[y_trans == h,]) - colMeans(x_train))
    }
    
    M <- Ehat %*% diag(phat) %*% t(Ehat)
    Sigma <- cov(x_train)
    S <- 1:p
    
    if (is.invertible(Sigma)) {
      eigen <- geigen_top(M, Sigma, d)
    } else {
      eigen <- geigen_top(M, Sigma + 1e-8 * diag(p), d)
    }
    
    init_sum <- apply(eigen$vectors, 1, function(x) sum(x^2))
    A_init <- order(init_sum, decreasing = TRUE)
    
    for(s in s_seq){
      fit <- sir_splicing_once(M, Sigma, s, d, H, A_init[1:s])
      loss_s <- energy::dcor(y_test, x_test %*% fit$B)
      avg_loss[which(s_seq == s)] <- avg_loss[which(s_seq == s)] + loss_s
    }
  }
  
  avg_loss <- avg_loss / K
  best_s <- s_seq[which.max(avg_loss)[1]]
  
  # on full data
  ## initialize
  p <- ncol(x)
  
  if (categorical) {
    y_trans <- y 
    y_unique <- as.numeric(names(table(y_trans)))
    H <- length(y_unique)
    
    phat <- sapply(y_unique, function(h) mean(y_trans == h))
    Ehat <- sapply(y_unique, function(h) colMeans(x[y_trans == h,]) - colMeans(x))
  } else {
    y_trans <- convert_slice(y, H)
    
    phat <- sapply(1:H, function(h) mean(y_trans == h))
    Ehat <- sapply(1:H, function(h) colMeans(x[y_trans == h,]) - colMeans(x))
  }
  
  M <- Ehat %*% diag(phat) %*% t(Ehat)
  Sigma <- cov(x)
  S <- 1:p
  
  if (is.invertible(Sigma)) {
    eigen <- geigen_top(M, Sigma, d)
  } else {
    eigen <- geigen_top(M, Sigma + 1e-8 * diag(p), d)
  }
  
  init_sum <- apply(eigen$vectors, 1, function(x) sum(x^2))
  A_init <- order(init_sum, decreasing = TRUE)
  
  fit <- sir_splicing_once(M, Sigma, best_s, d, H, A_init[1:best_s])
  
  # Return a list containing all the results and the best sparsity level
  return(list(
    fit = fit,
    best_s = best_s,
    avg_loss = avg_loss
  ))
}

