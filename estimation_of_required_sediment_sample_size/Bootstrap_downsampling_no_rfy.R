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


#split total and bulk samples
OTU_table_bulk <- OTU_table[,c(1:18)]
OTU_table_total <- OTU_table[,c(19:37)]



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

#specify division factors
division_factors <- 1:19


#no subsampling, no rarefaction!!!!!

# Loop over division factors to save multiple identical files without rarefaction
for (i in seq_along(division_factors)) {
  
  # Use the full combined OTU table as is, no rarefaction
  data_to_save <- rarefied_data_total
  
  # Optionally remove OTUs (columns) with zero counts across all samples
  nonzero_columns <- colSums(data_to_save) != 0
  data_to_save <- data_to_save[, nonzero_columns]
  
  # Save the data to a CSV file with unique file name
  filename <- paste0("no_rarefaction_result_", i, ".csv")
  write.csv(data_to_save, file = filename, row.names = TRUE)
}




bootstrap_bulk_t <- t(OTU_table_bulk_t) %>% as.data.frame()
write.csv(bootstrap_bulk_t, "bulk_bootstrap.csv", row.names = TRUE)

#remove zero values
nonzero_columns_bulk <- colSums(OTU_table_bulk_t) != 0
final_data_bulk <- OTU_table_bulk_t[, nonzero_columns_bulk]
final_bulk_t <- t(final_data_bulk) %>% as.data.frame()
write.csv(final_bulk_t, "bulk_bootstrap_no_rfy.csv", row.names = TRUE)
