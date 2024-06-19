library(glmnet)
library(ROCR)
library(msda)
library(energy)
library(pROC)
library(spls)
library(dplyr)
library(tidyverse)
library(GGally)
library(scales)
library(ggpubr)

source("source_realdata_multi.R")
source("seas/utility.R")
source("seas/seas.R")
source("sir_splicing.R")

# evaluation functions ---------------------------------------------------------------

evaluate_logistic_regression <- function(lp_train, y_train, lp_test, y_test) {
  
  lp_train_matrix <- as.matrix(lp_train)
  lp_test_matrix <- as.matrix(lp_test)
  
  glmnet_model <- glmnet(lp_train_matrix, y_train, family = "multinomial", lambda = 0)
  
  probabilities <- predict(glmnet_model, newx = lp_test_matrix, type = "response", s = 0)[,,1]
  class_labels <- sort(unique(y_train))
  predictions <- sapply(apply(probabilities, 1, which.max), function(x) class_labels[x])
  misclass_rate <- mean(predictions != y_test)
  
  auc_values <- sapply(class_labels, function(class) {
    class_index <- which(class_labels == class)
    roc_obj <- roc(as.numeric(y_test == class), probabilities[, class_index])
    auc(roc_obj)
  })
  
  mean_auc <- mean(auc_values)
  
  return(c(mce = misclass_rate, auc = mean_auc))
}

sim_splicing <- function(x_train, y_train, x_test, y_test, s, d, H = 5) {
  time <- system.time(
    fit_cv <- sir_splicing_cv(x_train, y_train, s, d, K = 5, categorical = T)
  )[3]
  plot(fit_cv$avg_loss)
  fit <- fit_cv$fit
  B <- fit$B
  A <- fit$A
  size <- length(A)
  
  lp_train <- x_train %*% B
  lp_test <- x_test %*% B
  
  eval <- evaluate_logistic_regression(lp_train, y_train, lp_test, y_test)
  mce <- eval[1]
  auc <- eval[2]
  
  return(list(mce = mce, auc = auc, size = size, time = time, lp_train = lp_train))
} 


sim_lasso_sir <- function(x_train, y_train, x_test, y_test, d){
  time <- system.time(fit <- LassoSIR::LassoSIR(x_train, y_train, categorical = TRUE, choosing.d = "manual", no.dim = d))[3]
  B <- fit$beta
  A <- which(abs(rowSums(B)) > 0)
  size <- length(A)
  
  lp_train <- x_train %*% B
  lp_test <- x_test %*% B
  
  eval <- evaluate_logistic_regression(lp_train, y_train, lp_test, y_test)
  mce <- eval[1]
  auc <- eval[2]
  
  return(list(mce = mce, auc = auc, size = size, time = time, lp_train = lp_train))
}

sim_seas <- function(x_train, y_train, x_test, y_test, d, H = 5) {
  time <- system.time(fit <- cv.seas(x_train, y_train, d = d, categorical = T))[3]
  B <- fit$beta
  A <- which(abs(rowSums(B)) > 0)
  size <- length(A)
  
  lp_train <- x_train %*% B
  lp_test <- x_test %*% B
  
  eval <- evaluate_logistic_regression(lp_train, y_train, lp_test, y_test)
  mce <- eval[1]
  auc <- eval[2]
  
  return(list(mce = mce, auc = auc, size = size, time = time, lp_train = lp_train))
}

sim_once <- function(seed, xall, y) {
  set.seed(seed)
  taskid <- seed
  
  x_train <- xall
  y_train <- y
  x_test <- xall
  y_test <- y
  
  smax <- 40
  d <- 2
  
  set.seed(1)
  res_seas <- sim_seas(x_train, y_train, x_test, y_test, d)
  set.seed(1)
  res_lasso_sir <- sim_lasso_sir(x_train, y_train, x_test, y_test, d)
  set.seed(1)
  res_splicing <- sim_splicing(x_train, y_train, x_test, y_test, (d+1):smax, d)
  
  return(list(res_splicing = res_splicing, res_seas = res_seas, res_lasso_sir = res_lasso_sir))
}



# pre-processing   --------------------------------------------------------------

set.seed(1)
data("lymphoma")
x <- lymphoma$x
y <- lymphoma$y
y <- as.numeric(y) 

y_dummy <- sapply(unique(y), function(i){
  as.numeric(y == i)
})

# Variable screen: keep 200 variables.
dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
  dcor(y_dummy, x[,i])
})
ord <- order(dist_cor, decreasing = TRUE)[1:100]
x <- x[,ord, drop=FALSE]

X <- x <- x %>% as.matrix() 



# comparison --------------------------------------------------------------

res <- sim_once(1, xall = X, y = y)
save(res, file = "results/res_realdata.rda")

# visualization ----------------------------------------------------------------
draw_plot <- function(lp) {
  
  glmnet_model <- glmnet(lp, y, family = "multinomial", lambda = 0)
  
  
  df_map <- tibble(
    y = c(0,1,2),
    Type = factor(c("DLBCL", "FL", "CLL"))
  )
  
  tbl <- tibble(lp1 = lp[,1], lp2 = lp[,2], y = y) %>%
    inner_join(df_map) %>%
    select(-y)
  
  lp1lim <- scales::expand_range(c(min(tbl$lp1), max(tbl$lp1)), mul = 0.1)
  lp2lim <- scales::expand_range(c(min(tbl$lp2), max(tbl$lp2)), mul = 0.1)
  
  # Create sequences within the expanded ranges
  lp1_seq <- seq(lp1lim[1], lp1lim[2], length.out = 300)
  lp2_seq <- seq(lp2lim[1], lp2lim[2], length.out = 300)
  
  grid_range <- expand.grid(lp1 = lp1_seq, lp2 = lp2_seq)
  probabilities <- predict(glmnet_model, newx = as.matrix(grid_range), type = "response")
  
  predicted_class <- apply(probabilities[,,1], 1, which.max)
  grid_range$predicted_class <- predicted_class
  
  p_scatter <- ggplot() +
    geom_raster(data = grid_range, aes(x = lp1, y = lp2, fill = factor(predicted_class)), alpha = 0.2) +
    geom_point(data = tbl, aes(x = lp1, y = lp2, color = Type, shape = Type), size = 2.5) +
    scale_fill_manual(values = hue_pal()(3)[c(2,3,1)], guide=F) +
    theme_minimal() +
    labs(x = "Linear Predictor 1", y = "Linear Predictor 2")+
    scale_x_continuous(limits = lp1lim, expand=c(0,0)) +
    scale_y_continuous(limits = lp2lim, expand=c(0,0)) 
  
  return(p_scatter)
}

p_splicing <- draw_plot(res$res_splicing$lp_train) +
  labs(title = "Splicing-SIR") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14)) 

p_seas <- draw_plot(res$res_seas$lp_train) +
  labs(title = "SEAS-SIR")+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))

p_lasso <- draw_plot(res$res_lasso_sir$lp_train) +
  labs(title = "LASSO-SIR")+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))

p_scatter <- ggarrange(p_splicing, p_seas, p_lasso, common.legend = T, legend = "right")
ggsave(p_scatter, filename = "results/p_scatter_realdata.pdf", width = 12, height = 8)