# pairwise PERMANOVA
# Lukas Damm
# RPTU Kaiserslautern-Landau
# 24.07.2025

library(dplyr)
library(vegan)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(devtools)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", lib = my_lib)
library(pairwiseAdonis, lib.loc = my_lib)

#load mapping files containing sample IDs and coordinates

factor_total_layers <- read.csv("Stations_coordinates_total.csv")
factor_total_layers$core_id <- as.factor(factor_total_layers$core_id)

env <- read.csv("Stations_coordinates_bulk_total_comb.csv")

#load OTU tables

OTU_table_total_sep_paper <- read.csv("OTU_table_total_final.csv") %>% column_to_rownames("X") %>% t() %>% as.data.frame()
OTU_table_comb_paper <- read.csv("OTU_table_all_comb_final.csv", sep = "\t") %>% t() %>% as.data.frame()


#calculate Bray-Curtis distance

gen_dist_total_sep_paper <- vegdist(OTU_table_total_sep_paper, method = "bray") %>% as.matrix()
gen_dist_comb_paper <- vegdist(OTU_table_comb_paper, method = "bray") %>% as.matrix()


#perform pairwise adonis

set.seed(1)
pairwise_total_bulk <- pairwise.adonis2(gen_dist_comb_paper ~ Extraction_method + Longitude + Latitude + Longitude * Latitude, data = env, nperm = 9999, by = "terms")
print(pairwise_total_bulk)
write.table(pairwise_total_bulk, "permanova_pairwise_total_bulk.txt", sep = "\t")


set.seed(2)
pairwise_total_layers <- pairwise.adonis2(gen_dist_total_sep_paper ~ Layer + Longitude + Latitude + Longitude * Latitude, data = factor_total_layers, strata = factor_total_layers$core_id, nperm = 9999, by = "terms")
print(pairwise_total_layers)
write.csv(pairwise_total_layers, "permanova_pairwise_total_layers.csv")


