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
OTU_table <- read.csv("OTU_table_all_final.csv") %>% column_to_rownames("X")
head(OTU_table)
colnames(OTU_table)


#split total and bulk samples
OTU_table_bulk <- OTU_table[,c(1:17)]
OTU_table_total <- OTU_table[,c(18:66)]



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

# rarefy total samples prior to bootstrap draw. Do not rarefy because of low read counts

#rarefaction_level_total <- 15005

#OTU_table_total_rare <- rrarefy(OTU_table_total_t, rarefaction_level_total)

#combine rarefied layers
OTU_table_total_combined <- OTU_table_total %>% 
  mutate(S10 = S10_1 + S10_2 + S10_3) %>%  
  mutate(S11 = S11_1 + S11_2 + S11_3) %>%
  mutate(S13 = S13_1 + S13_2) %>%
  mutate(S16 = S16_1 + S16_2 + S16_3) %>%
  mutate(S17 = S17_1 + S17_3) %>%
  mutate(S21 = S21_1 + S21_2 + S21_3) %>%
  mutate(S23 = S23_1 + S23_3) %>%
  mutate(S24 = S24_1 + S24_2 + S24_3) %>%
  mutate(S25 = S25_1 + S25_2 + S25_3) %>%
  mutate(S26 = S26_1 + S26_2 + S26_3) %>%
  mutate(S27 = S27_1 + S27_3) %>%
  mutate(S03 = S03_2 + S03_3) %>%
  mutate(S33 = S33_1 + S33_2 + S33_3) %>%
  mutate(S34 = S34_1 + S34_2 + S34_3) %>%
  mutate(S38 = S38_3) %>%
  mutate(S04 = S04_1 + S04_2 + S04_3) %>%
  mutate(S42 = S42_1 + S42_2 + S42_3) %>%
  mutate(S43 = S43_1 + S43_2 + S43_3) %>%
  mutate(S09 = S09_1+ S09_3) 

OTU_table_total_combined <- OTU_table_total_combined[,c(50:68)]

colSums(OTU_table_total_combined)
mean(OTU_table_total_comb_sum)
OTU_table_total_combined_t <- t(OTU_table_total_combined)

#specify the first rarefaction level for the bootstrap draw
first_rarefaction_level <- 6680 #lowest read counts of bulk samples


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
  rarefied_data <- rrarefy(OTU_table_total_combined_t, rarefaction_level)
  
  #remove zero values
  nonzero_columns <- colSums(rarefied_data) != 0
  rarefied_data <- rarefied_data[, nonzero_columns]
    
  #save the rarefied data into a single file
  filename <- paste0("rarefaction_result_", i, ".csv")
  write.csv(rarefied_data, file = filename, row.names = TRUE)
  
}



#rarefaction of bulk samples
rarefaction_level_bulk <- 6680

set.seed(1)
rarefied_bulk <- rrarefy(OTU_table_bulk_t, rarefaction_level_bulk)
#remove zero values
nonzero_columns_bulk <- colSums(rarefied_bulk) != 0
rarefied_data_bulk <- rarefied_bulk[, nonzero_columns_bulk]
rarefied_bulk_t <- t(rarefied_data_bulk)
write.csv(rarefied_bulk_t, "rarefied_bulk_bootstrap.csv", row.names = TRUE)

