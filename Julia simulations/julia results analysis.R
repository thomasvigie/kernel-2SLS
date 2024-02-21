# This script analyzes the results from the Julia simulations
library("tidyverse")

source("auxiliary functions.R")
source("DGP parameters.R")

# Load the results
results <- read_csv("first results (most recent).csv")

MC_b <- results[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "bias" ) 

MC_v <- results[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "variance" )

MC_MSE <- results[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "MSE" )

MC_medAD <- results[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "medAD" )

MC_medbias <- results[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "medbias" )
