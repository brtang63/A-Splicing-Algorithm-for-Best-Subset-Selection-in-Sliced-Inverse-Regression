library(snowfall)
source("seas/utility.R")
source("seas/seas.R")
source("sir_splicing.R")

path_results <- "results/"
if (!dir.exists(path_results)) {
  dir.create(path_results)
}

# generate data -----------------------------------------------------------
generate_data <- function(seed = 1, n = 400, p = 20, rho = 0.5, type_model = 1) {
  set.seed(seed)
  
  H <- 5
  TA <- c(1:8)
  s <- length(TA)
  
  if (type_model %in% c(1,2)) {
    TB1 <- TB2 <- rep(0, p)
    TB1[TA] <- c(rep(0.5, 4), rep(0,4))
    TB2[TA] <- c(rep(0, 4), rep(0.5,4))
    TB <- cbind(TB1, TB2)
  }
  
  if (type_model %in% c(3,4)) {
    TB <- matrix(rep(0, p), ncol = 1) 
    TB[TA, 1] <- rep(0.5, 8)
  }
  
  Sigma <- matrix(0, p, p)
  Sigma <- rho^(abs(row(Sigma) - col(Sigma)))
  diag(Sigma) <- 1
  
  x <- MASS::mvrnorm(n, rep(0, p), Sigma)
  
  epsilon1 <- rnorm(n) 
  epsilon2 <- rcauchy(n) 
  
  y <- switch (type_model,
               (x %*% TB1) / (0.5 + (1.5 + x %*% TB2)^2) + 0.2 * epsilon1,
               (x %*% TB1) * exp(x %*% TB2 + 0.5 * epsilon1),
               x %*% TB + epsilon1,
               exp(x %*% TB) + epsilon1
  )
  
  return(list(x = x, y = y, TA = TA, TB = TB))
}

# simulation functions ---------------------------------------------------------------
sim_splicing <- function(x, y, s, d, TA, TB, H = 5) {
  time <- system.time(
    fit_cv <- sir_splicing_cv(x, y, s, d, K = 5)
  )[3]
  fit <- fit_cv$fit
  B <- fit$B
  A <- fit$A
  p <- nrow(B)
  corr <- cancor(B, TB)$cor[1]
  
  size <- length(A)
  TP <- length(intersect(A, TA))
  FP <- length(setdiff(A, TA))
  FN <- length(setdiff(TA, A))
  TPR <- TP / (TP + FN)
  FPR <- FP / (p - length(TA))
  
  return(c(TPR = TPR, FPR = FPR, corr = corr, size = size, time = time))
} 

sim_lasso_sir <- function(x, y, d, TA, TB){
  time <- system.time(fit <- LassoSIR::LassoSIR(x, y, H = 5, choosing.d = "manual", no.dim = d))[3]
  B <- fit$beta
  A <- which(abs(rowSums(B)) > 0)
  p <- nrow(B)
  corr <- cancor(B, TB)$cor[1]
  
  size <- length(A)
  TP <- length(intersect(A, TA))
  FP <- length(setdiff(A, TA))
  FN <- length(setdiff(TA, A))
  TPR <- TP / (TP + FN)
  FPR <- FP / (p - length(TA))
  
  return(c(TPR = TPR, FPR = FPR, corr = corr, size = size, time = time))
}

sim_seas <- function(x, y, d, TA, TB) {
  time <- system.time(fit <- cv.seas(x, y, d = d))[3]
  B <- fit$beta
  A <- which(abs(rowSums(B)) > 0)
  p <- nrow(B)
  corr <- cancor(B, TB)$cor[1]
  
  size <- length(A)
  TP <- length(intersect(A, TA))
  FP <- length(setdiff(A, TA))
  FN <- length(setdiff(TA, A))
  TPR <- TP / (TP + FN)
  FPR <- FP / (p - length(TA))
  
  return(c(TPR = TPR, FPR = FPR, corr = corr, size = size, time = time))
}

sim_once_safe <- function(seed = 1, n, p, type_model = 1, rho = 0) {
  H <- 5
  smax <- 15
  
  if (type_model %in% c(1,2))
    d <- 2
  if (type_model %in% c(3,4))
    d <- 1
  
  data <- generate_data(seed = seed, n = n, p = p, rho = rho, type_model = type_model)
  x <- data$x
  y <- data$y
  TA <- data$TA
  TB <- data$TB
  
  res_splicing <- tryCatch({
    sim_splicing(x, y, (d+1):smax, d, TA, TB)
  }, error = function(e) {
    rep(NA, 5)
  })
  
  res_lasso_sir <- tryCatch({
    sim_lasso_sir(x, y, d, TA, TB)
  }, error = function(e) {
    rep(NA, 5)
  })
  
  res_seas <- tryCatch({
    sim_seas(x, y, d, TA, TB)
  }, error = function(e) {
    rep(NA, 5)
  })
  
  return(rbind(res_splicing, res_seas, res_lasso_sir))
}


# simulation --------------------------------------------------------------

M <- 200
sfInit(parallel = TRUE, cpus = 20) 
sfLibrary(geigen)
sfLibrary(LassoSIR)
sfLibrary(msda)
sfLibrary(RSpectra)
sfLibrary(energy)
sfExportAll()

p_seq <- 600
n_seq <- 500 * 2:6
rho_seq <- c(0, 0.5)
type_model_seq <- c(1, 2, 3, 4)

for (p in p_seq) {
  for (type_model in type_model_seq) {
    for (rho in rho_seq) {
      for (n in n_seq) {
        res <- sfLapply(1:M, sim_once_safe, n = n, p = p, type_model = type_model, rho = rho)
        res <- simplify2array(res)
        save(res, file = paste0(path_results, paste("res", type_model, n, p, rho, sep = "_"), ".rda"))
        cat(n, "done!", "\n")
      }
    }
  }
}



