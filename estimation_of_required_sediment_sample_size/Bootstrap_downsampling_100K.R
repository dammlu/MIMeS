# original code from Verena rubel https://github.com/verubel/MIMeS/blob/main/Bulk_vs_total_eDNA/estimate_req_sample_n_bootstrapping.R
# bootstrap combination of samples to reach disired level of diversity (ASV richness)
# Verena Rubel 
# RPTU Kaiserslautern Landau
# 27.02.2024
# new: account for increase of sequencing depth by combining ASV lists of all samples
# idea: rarefy all tables to x/1-x/19 so 1x-19x sample combination leads to read counts at same level
# code modified by Lukas Damm
# RPTU Kaiserslautern-Landau
# 10.07.2025

library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(tidyverse)
library(data.table)



#load OTU table
OTU_table <- read.csv("aggML_no_rfy.csv") #%>% column_to_rownames("X")
OTU_table <- OTU_table[, c(13:50)] %>% column_to_rownames("OTU_ID")
head(OTU_table)
colnames(OTU_table)
OTU_table_t <- t(OTU_table)

#rarefy dataset to 100.000 reads

rarefaction_level <- 100000
OTU_table_rarefied <- rrarefy(OTU_table_t, rarefaction_level)
OTU_table_final <- t(OTU_table_rarefied)

#split total and bulk samples
OTU_table_bulk <- OTU_table_final[,c(1:18)]
OTU_table_total <- OTU_table_final[,c(19:37)]



#check min read count
OTU_table_bulk_sum <- as.data.frame(colSums(OTU_table_bulk))
colnames(OTU_table_bulk_sum)[1]<-"reads"
min(OTU_table_bulk_sum)
max(OTU_table_bulk_sum)
OTU_table_total_sum <- as.data.frame(colSums(OTU_table_total))
colnames(OTU_table_total_sum)[1]<-"reads"
min(OTU_table_total_sum)
max(OTU_table_total_sum)
OTU_table_bulk_t <- t(OTU_table_bulk) %>% as.data.frame()
OTU_table_total_t <- t(OTU_table_total) %>% as.data.frame()

#remove zero values from total set
nonzero_columns_total <- colSums(OTU_table_total_t) != 0
rarefied_data_total <- OTU_table_total_t[, nonzero_columns_total]

#specify the first rarefaction level for the bootstrap draw
first_rarefaction_level <- 100000


#specify division factors
division_factors <- 1:19


#calculate rarefaction levels based on division factors

rarefaction_levels <- first_rarefaction_level / division_factors

#print rarefaction levels

print("Rarefaction levels:")
print(rarefaction_levels)

# create loop for rarefaction on different levels
for (i in seq_along(division_factors)) {
  
  rarefaction_level <- first_rarefaction_level / division_factors[i]
  
  #perform rarefaction
  set.seed(1)
  rarefied_data <- rrarefy(rarefied_data_total, rarefaction_level)
  
  #remove zero values
  nonzero_columns <- colSums(rarefied_data) != 0
  rarefied_data <- rarefied_data[, nonzero_columns]
  
  #save the rarefied data into a single file
  filename <- paste0("rarefaction_100K_result_", i, ".csv")
  write.csv(rarefied_data, file = filename, row.names = TRUE)
  
}

#remove zero values
nonzero_columns_bulk <- colSums(OTU_table_bulk_t) != 0
rarefied_data_bulk <- OTU_table_bulk_t[, nonzero_columns_bulk]
rarefied_bulk_t <- t(rarefied_data_bulk) %>% as.data.frame()
write.csv(rarefied_bulk_t, "bulk_bootstrap_100K.csv", row.names = TRUE)
