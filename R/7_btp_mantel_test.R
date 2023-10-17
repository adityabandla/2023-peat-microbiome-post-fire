## name: 7_btp_mantel_test.R
## purpose: to examine associations between pore water quality and 
## homogeneous selection across peat depth
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## start of code run on NSCC (National Super Computing Center, Singapore)
## load packages
suppressPackageStartupMessages({
  library(ape)
  library(iCAMP)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
})

## source functions
source("tools.r")

## run Mantel tests
btp_turnover_bacteria <- read.csv("2020_btp_bacteria_process_importance_each_turnover.csv")
btp_turnover_archaea <- read.csv("2020_btp_archaea_process_importance_each_turnover.csv")
treatment <- "^B_"
turnover <- btp_turnover_archaea

## subset turnovers corresponding to group and process
y.3col <- turnover %>%
  filter(str_detect(samp1, treatment), str_detect(samp2, treatment)) %>%
  select(samp1, samp2, HoS)

## import environmental data
env <- read.csv("2020_btp_depth_env_data.csv") %>%
  select(-one_of("psf_type", "depth", "TDS", "salinity")) %>%
  pivot_longer(cols = !contains("sample"), names_to = "variable", values_to = "value")

## log transform all variable except pH
env_log <- read.csv("2020_btp_depth_env_data.csv", header = TRUE) %>%
  select(-one_of("psf_type", "depth", "salinity", "TDS")) %>%
  mutate_at(vars(-pH, -sample), log) %>%
  pivot_longer(cols = !contains("sample"), names_to = "variable", values_to = "value")

## compute variable means for each pair of samples
## untransformed data
env_mean <- y.3col %>%
  select(samp1, samp2) %>%
  inner_join(., env, by = c("samp1" = "sample"), multiple = "all", relationship = "many-to-many") %>%
  inner_join(., env, by = c("samp2" = "sample"), suffix = c("_x", "_y"), multiple = "all", relationship = "many-to-many") %>%
  filter(variable_x == variable_y) %>%
  group_by(samp1, samp2, variable_x, variable_y) %>%
  mutate(mean = mean(c(value_x, value_y))) %>%
  ungroup %>%
  select(samp1, samp2, variable_x, mean) %>%
  mutate(variable_x = paste0("m_", variable_x)) %>%
  pivot_wider(names_from = "variable_x", values_from = mean)

## log transformed data
env_log_mean <- y.3col %>%
  select(samp1, samp2) %>%
  inner_join(., env, by = c("samp1" = "sample"), multiple = "all", relationship = "many-to-many") %>%
  inner_join(., env, by = c("samp2" = "sample"), suffix = c("_x", "_y"), multiple = "all", relationship = "many-to-many") %>%
  filter(variable_x == variable_y) %>%
  group_by(samp1, samp2, variable_x, variable_y) %>%
  mutate(mean = mean(c(value_x, value_y))) %>%
  ungroup %>%
  select(samp1, samp2, variable_x, mean) %>%
  mutate(variable_x = paste0("m_", variable_x)) %>%
  pivot_wider(names_from = "variable_x", values_from = mean)

## compute variable differences for each pair of samples
## untransformed data
env_diff <-  y.3col %>%
  select(samp1, samp2) %>%
  inner_join(., env, by = c("samp1" = "sample"), multiple = "all", relationship = "many-to-many") %>%
  inner_join(., env, by = c("samp2" = "sample"), suffix = c("_x", "_y"), multiple = "all", relationship = "many-to-many") %>%
  filter(variable_x == variable_y) %>%
  group_by(samp1, samp2, variable_x, variable_y) %>%
  mutate(difference = abs(value_x - value_y)) %>%
  ungroup %>%
  select(samp1, samp2, variable_x, difference) %>%
  mutate(variable_x = paste0("d_", variable_x)) %>%
  pivot_wider(names_from = "variable_x", values_from = difference)

## log transformed data
env_log_diff <-  y.3col %>%
  select(samp1, samp2) %>%
  inner_join(., env, by = c("samp1" = "sample"), multiple = "all", relationship = "many-to-many") %>%
  inner_join(., env, by = c("samp2" = "sample"), suffix = c("_x", "_y"), multiple = "all", relationship = "many-to-many") %>%
  filter(variable_x == variable_y) %>%
  group_by(samp1, samp2, variable_x, variable_y) %>%
  mutate(difference = abs(value_x - value_y)) %>%
  ungroup %>%
  select(samp1, samp2, variable_x, difference) %>%
  mutate(variable_x = paste0("d_", variable_x)) %>%
  pivot_wider(names_from = "variable_x", values_from = difference)

## combine tables 
mdenv <- inner_join(env_mean, env_diff)
mdenv_log <- inner_join(env_log_mean, env_log_diff)

## cycle through depth-wise
y.3col_surface <- y.3col %>%
  filter(str_detect(samp1, paste0(treatment, "A_")), str_detect(samp2, paste0(treatment, "A_")))

## cycle through depth-wise
y.3col_mid <- y.3col %>%
  filter(str_detect(samp1, paste0(treatment, "P_")), str_detect(samp2, paste0(treatment, "P_")))

## cycle through depth-wise
y.3col_deep <- y.3col %>%
  filter(str_detect(samp1, paste0(treatment, "B_")), str_detect(samp2, paste0(treatment, "B_")))

## cycle through depth-wise
mdenv_surface <- mdenv %>% filter(str_detect(samp1, paste0(treatment, "A_")), str_detect(samp2, paste0(treatment, "A_")))
mdenv_log_surface <- mdenv_log %>% filter(str_detect(samp1, paste0(treatment, "A_")), str_detect(samp2, paste0(treatment, "A_")))

## cycle through depth-wise
mdenv_mid <- mdenv %>% filter(str_detect(samp1, paste0(treatment, "P_")), str_detect(samp2, paste0(treatment, "P_")))
mdenv_log_mid <- mdenv_log %>% filter(str_detect(samp1, paste0(treatment, "P_")), str_detect(samp2, paste0(treatment, "P_")))

## cycle through depth-wise
mdenv_deep <- mdenv %>% filter(str_detect(samp1, paste0(treatment, "B_")), str_detect(samp2, paste0(treatment, "B_")))
mdenv_log_deep <- mdenv_log %>% filter(str_detect(samp1, paste0(treatment, "B_")), str_detect(samp2, paste0(treatment, "B_")))

## run Mantel test
mdenv_X <- mdenv_deep
y.3col_X <- y.3col_deep

mcm <- t(sapply(
  3:ncol(mdenv_X),
  function(i) {
    x.3col <- mdenv_X[, c(1, 2, i)]
    mci <- mcMantel(
      y.3col = y.3col_X, 
      y.ids = NULL,
      x.3col = x.3col, 
      grp.rand = NULL,
      grp.const = NULL, 
      try.time = 5,
      method = "pearson", 
      permutations = 999
    )
    outi <- c(r = mci$statistic, p = mci$signif)
  }
))

rownames(mcm) = colnames(mdenv)[3:ncol(mdenv)]
colnames(mcm) = c("r","P")