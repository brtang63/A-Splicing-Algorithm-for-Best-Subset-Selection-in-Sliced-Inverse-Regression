library(tidyverse)
library(cubelyr)
library(ggpubr)
library(grid)
library(ggthemes)

type_model <- 1 
# type_model <- 2
# type_model <- 3 
# type_model <- 4

path_results <- "results/"
p_seq <- c(600)
n_seq <- 500 * 2:6

rho_seq <- c(0, 0.5)
M <- 200

convert_to_arrays <- function(input_list) {
  for (i in seq_along(input_list)) {
    if (!is.array(input_list[[i]])) {
      input_list[[i]] <- array(NA, dim = c(4, 5))
      colnames(input_list[[i]]) <- c("TPR","FPR","corr", "size", "time.elapsed")
      rownames(input_list[[i]]) <- c("res_splicing","res_lasso_sir", "res_seas")
    }
  }
  return(input_list)
}


processed_res <- NULL
for (rho in rho_seq) {
  for (p in p_seq) {
    for (n in n_seq) {
      load(paste0(path_results, paste("res", type_model, n, p, rho, sep = "_"), ".rda"))
      if (is.array(res)) tmp <- res else {
        tmp <- simplify2array(convert_to_arrays(res))  
        print(1)
      }  
      dimnames(tmp)[[3]] <- 1:M
      names(dimnames(tmp)) <- c("Method", "Measure", "Replicates")
      tmp <- tmp %>%
        as.tbl_cube(met_name = "value") %>%
        as_tibble %>%
        add_column(n = n, p = p, rho = rho) %>%
        relocate(value)
      processed_res <- rbind(processed_res, tmp)
    }
  }
}

size_text <- 24
wd <- 12
ht <- 8
palette_hue <- (scales::hue_pal())(10)
size_list <- c(16, 17, 15, 3, 7, 8, 9, 21, 25, 24)

processed_res$rho <- factor(processed_res$rho, levels = c("0", "0.5"))
processed_res$Measure <- factor(processed_res$Measure, levels = c("TPR", "FPR", "corr", "size", "time.elapsed"))
processed_res <- processed_res %>%
  drop_na() %>%
  mutate(Method = recode(Method,
                         "res_splicing" = "Splicing-SIR",
                         "res_lasso_sir" = "LASSO-SIR",
                         "res_seas" = "SEAS-SIR"))  %>%
  mutate(rho = recode(rho,
                      "0" = "Independent",
                      "0.5" = "Correlated"))
processed_res$Method <- factor(processed_res$Method, levels = c("Splicing-SIR", "SEAS-SIR", "LASSO-SIR"))

tbl_active <- processed_res %>%
  filter(Measure %in% c("TPR")) %>%
  group_by(Method, p, n, rho, value) %>%
  dplyr::summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  complete(Method, p, n, rho, value, fill = list(N = 0, freq = 0)) %>%
  filter(value == 1) %>%
  add_column(ratetype = "active set")

tbl_inactive <- processed_res %>%
  filter(Measure %in% c("FPR")) %>%
  group_by(Method, p, n, rho, value) %>%
  dplyr::summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  complete(Method, p, n, rho, value, fill = list(N = 0, freq = 0)) %>%
  filter(value == 0) %>%
  add_column(ratetype = "inactive set")

tbl_exact <- processed_res %>%
  filter(Measure %in% c("TPR", "FPR")) %>%
  pivot_wider(names_from = Measure, values_from = value) %>%
  mutate(TR = TPR - FPR) %>%
  mutate(TPR = NULL, TNR = NULL) %>%
  group_by(Method, p, n, rho, TR) %>%
  dplyr::summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  complete(Method, p, n, rho, TR, fill = list(N = 0, freq = 0)) %>%
  filter(TR == 1) %>%
  dplyr::rename(value = TR) %>%
  add_column(ratetype = "exact")

## tibble of support recovery rate
tbl <- bind_rows(tbl_active, tbl_inactive, tbl_exact) %>%
  mutate(ratetype = factor(ratetype, levels = c("active set", "inactive set", "exact")))

p_rate <- tbl %>%
  ggplot(., aes(x = n, y = freq, color = Method)) + 
  facet_grid(cols = vars(rho), rows = vars(ratetype), switch = "y") +
  geom_point(aes(shape = Method)) + geom_line() + 
  theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = size_text),
        legend.key.size = unit(1, "cm"),
        panel.spacing = unit(2, "lines")) +
  ylab("Subset recovery probability") + 
  scale_color_manual(values = palette_hue) + scale_shape_manual(values = size_list) +
  scale_x_continuous(
    # breaks = c(1000, 2000, 3000),
    expression(n)) +
  guides(color=guide_legend(nrow =1))

# ggsave(filename = "model1.pdf", p1, width = 10, height = 6)
p_size <- processed_res %>%
  filter(Measure %in% c("size")) %>%
  mutate(n = as.factor(n)) %>%
  ggplot(., aes(x = n, y = value, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~rho) +
  scale_fill_manual(values = palette_hue) + scale_shape_manual(values = size_list) +
  theme_bw() +
  theme(text = element_text(size = size_text),
        panel.spacing = unit(2, "lines")) +
  geom_hline(yintercept = 8, color = "red")+
  scale_color_manual(values = palette_hue) +
  ylab("Sparsity level") +
  guides(color=guide_legend(nrow =1))

p_cor <- processed_res %>%
  filter(Measure %in% c("corr")) %>%
  group_by(n, p, rho, Method) %>%
  summarise(med = median(value), .groups = "drop") %>%
  ggplot(., aes(x = n, y = med)) +
  facet_wrap(~rho)+
  geom_point(aes(color = Method)) +
  geom_line(aes(color = Method)) +
  theme_bw()+
  theme(text = element_text(size = size_text),
        panel.spacing = unit(2, "lines")) +
  scale_color_manual(values = palette_hue) +
  ylab("Canonical correlation") +
  guides(fill=guide_legend(nrow =1))

p_0 <- ggarrange(p_rate, ggarrange(p_cor, p_size, ncol = 1, legend = "none", labels = c("b", "c"), font.label = list(size = 24)), ncol = 2, common.legend = T, labels = c("a", ""), font.label = list(size = 24), legend = "bottom") 
ggsave(p_0, filename = paste0("results/performance_model", type_model, ".pdf"), width = 16, height = 9)
