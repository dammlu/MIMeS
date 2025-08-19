# Lukas Damm
# RPTU Kaiserslautern-Landau



library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)

#  prepare OTU table from rarefied ASV table

#load in OTU table from Sahar, add OTU names column beforehand____________________________________________________________________________________

#load OTU (97% sim) table
OTU_table_rfy <- read.csv("OTU_97%_rarefied_OTU_final.csv", sep=",")
head(OTU_table_rfy)
colnames(OTU_table_rfy)
OTU_table_final_rfy <- OTU_table_rfy[, c(16:86)]

taxo_table_rfy <- OTU_table_rfy[, c(6:12,16)]

#load OTU table 97% with all taxo groups

OTU_table_taxo_all <- read.csv("OTU_97_taxo.csv", sep=",")
OTU_table_taxo_all <- OTU_table_taxo_all[, c(16:86)]

OTU_table_taxo_all_bulk <- OTU_table_taxo_all[, c(1,55:71)] %>% column_to_rownames("OTU")
OTU_table_taxo_all_total <- OTU_table_taxo_all[, c(1:54)] %>% column_to_rownames("OTU") %>%  mutate_at(vars(1:53), as.numeric)

OTU_table_taxo_all_bulk_0 <- OTU_table_taxo_all_bulk[which(rowSums(OTU_table_taxo_all_bulk)>0),]
OTU_table_taxo_all_total_0 <- OTU_table_taxo_all_total[which(rowSums(OTU_table_taxo_all_total)>0),]

#OTU_table_curation__________________________________________________________________________________________________________________________________


OTU_table_form_rfy <- OTU_table_final_rfy %>% column_to_rownames("OTU") %>%  mutate_at(vars(1:70), as.numeric)
OTU_table_form_t_rfy <- as.data.frame(t(OTU_table_form_rfy)) %>% mutate_all(~ifelse(is.na(.),0,.))

OTU_table_taxo_all_form <- OTU_table_taxo_all %>% column_to_rownames("OTU") %>%  mutate_at(vars(1:70), as.numeric)
OTU_table_taxo_all_form <- as.data.frame(t(OTU_table_taxo_all_form)) %>% mutate_all(~ifelse(is.na(.),0,.))

#filter out OTUs with 0 reads
OTU_table_form_t_0_rfy <- OTU_table_form_t_rfy[, colSums(OTU_table_form_t_rfy) !=0] %>% t()

OTU_table_taxo_all_form_0 <- OTU_table_taxo_all_form[, colSums(OTU_table_taxo_all_form) !=0] %>% t()


#get read and OTU counts
reads_after_OTU_clustering <- colSums(OTU_table_form_t_0_rfy)
write.csv(reads_after_OTU_clustering, "reads_after_clustering.csv")

reads_after_OTU_all_taxo <- colSums(OTU_table_taxo_all_form_0)
write.csv(reads_after_OTU_all_taxo, "reads_after_OTU_taxo_all.csv")

OTU_table_form_t_rfy_count <- as.data.frame(t(OTU_table_form_t_0_rfy))
count_positive <- function(row) {return(sum(row > 0))}
count_pos <- apply(OTU_table_form_t_rfy_count, 1, count_positive)
write.csv(count_pos, "OTUs_after_clustering.csv")

OTU_table_taxo_all_count <- as.data.frame(t(OTU_table_taxo_all_form_0))
count_positive <- function(row) {return(sum(row > 0))}
count_pos <- apply(OTU_table_taxo_all_count, 1, count_positive)
write.csv(count_pos, "OTUs_after_clust_taxo_all.csv")

#split total and bulk samples
OTU_table_total_rfy <- as.data.frame(OTU_table_form_t_0_rfy[, c(1:53)])
OTU_table_bulk_rfy <- as.data.frame(OTU_table_form_t_0_rfy[, c(54:70)])


#bulk
OTU_table_sum_bulk_rfy <- as.data.frame(colSums(OTU_table_bulk_rfy))
colnames(OTU_table_sum_bulk_rfy)[1]<-"reads"
OTU_table_sum_bulk_rfy <- rownames_to_column(OTU_table_sum_bulk_rfy, var = "sample")

#total
OTU_table_sum_total_rfy <- as.data.frame(colSums(OTU_table_total_rfy))
colnames(OTU_table_sum_total_rfy)[1]<-"reads"
OTU_table_sum_total_rfy <- rownames_to_column(OTU_table_sum_total_rfy, var = "sample")


#load in mapping file
env_rfy <- read.csv("mapping_file.csv")

#join mapping file and OTU file
OTU_mapped_df_bulk_rfy <- as.data.frame(OTU_table_sum_bulk_rfy) %>%
  inner_join(env_rfy, by = "sample") %>% 
  select(final_sample, reads, sampling) %>% 
  filter(sampling !="NA")

OTU_mapped_df_total_rfy <- as.data.frame(OTU_table_sum_total_rfy) %>%
  inner_join(env_rfy, by = "sample") %>% 
  select(final_sample, reads, sampling) %>% 
  filter(sampling !="NA")


#get min and max number of reads
min(OTU_mapped_df_bulk_rfy["reads"])
max(OTU_mapped_df_bulk_rfy["reads"])
OTU_mapped_df_row_bulk_rfy <- column_to_rownames(OTU_mapped_df_bulk_rfy, var = "final_sample")


min(OTU_mapped_df_total_rfy["reads"])
max(OTU_mapped_df_total_rfy["reads"])
OTU_mapped_df_row_total_rfy <- column_to_rownames(OTU_mapped_df_total_rfy, var = "final_sample")


OTU_mapped_df_bulk_total_rfy <- as.data.frame(t(OTU_mapped_df_row_total_rfy)) %>% cbind(t(OTU_mapped_df_row_bulk_rfy)) 
OTU_mapped_df_bulk_total_final_rfy <- as.data.frame(t(OTU_mapped_df_bulk_total_rfy)) %>% rownames_to_column("final_sample") %>% mutate_at(vars(reads), as.numeric)


#plot read counts
plot_read_counts_rfy<- ggplot(OTU_mapped_df_bulk_total_final_rfy, aes(x = sampling, y = reads, group = sampling))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("No. of reads")+ 
  geom_point(aes(colour = sampling, shape = sampling))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")
plot_read_counts_rfy
ggsave("plot_read_counts_OTUs_rfy.pdf",plot_read_counts_rfy, width = 8, height = 10)  



#join mapping file
OTU_table_bulk_map_rfy <- as.data.frame(t(OTU_table_bulk_rfy)) %>% rownames_to_column("sample") %>%
  inner_join(env_rfy, by = "sample") %>% column_to_rownames("final_sample") %>% select(-sample, -sampling, -new, -layer)

OTU_table_total_map_rfy  <- as.data.frame(t(OTU_table_total_rfy)) %>% rownames_to_column("sample") %>%
  inner_join(env_rfy, by = "sample") %>% column_to_rownames("final_sample") %>% select(-sample, -sampling, -new, -layer)

OTU_table_bulk_t_rfy <- t(OTU_table_bulk_rfy)
OTU_table_total_t_rfy <- t(OTU_table_total_rfy)



count_positive_rfy <- function(row) {return(sum(row > 0))}
count_pos_bulk_rare_rfy <- apply(OTU_table_bulk_rfy, 1, count_positive_rfy)
count_pos_total_rare_rfy <- apply(OTU_table_total_rfy, 1, count_positive_rfy)
write.csv(count_pos_bulk_rare_rfy, "OTUs_bulk_rfy.csv")
write.csv(count_pos_total_rare_rfy, "OTUs_total_rfy.csv")


#combine total and bulk samples
OTU_table_total_row_rfy <- as.data.frame(t(OTU_table_total_map_rfy)) %>% rownames_to_column("OTU")
OTU_table_bulk_row_rfy <- as.data.frame(t(OTU_table_bulk_map_rfy)) %>% rownames_to_column("OTU")
OTU_table_rare_all_rfy <- full_join(OTU_table_bulk_row_rfy,OTU_table_total_row_rfy, by = "OTU")
OTU_table_rare_all_final_rfy <- OTU_table_rare_all_rfy %>% mutate_at(vars(-OTU), ~ifelse(is.na(.), 0, .))





#combine columns of layers

OTU_table_total_rfy_t <- as.data.frame(t(OTU_table_total_map_rfy))
OTU_table_total_combined_rfy <- OTU_table_total_rfy_t %>% 
  mutate(S10 = S10_1 + S10_2 + S10_3) %>%  
  mutate(S11 = S11_1 + S11_2 + S11_3) %>%
  mutate(S13 = S13_1 + S13_2 + S13_2) %>%
  mutate(S16 = S16_1 + S16_2 + S16_3) %>%
  mutate(S17 = S17_1 + S17_2 + S17_3) %>%
  mutate(S21 = S21_1 + S21_2 + S21_3) %>%
  mutate(S23 = S23_1 + S23_2 + S23_3) %>%
  mutate(S24 = S24_1 + S24_2 + S24_3) %>%
  mutate(S25 = S25_1 + S25_2 + S25_3) %>%
  mutate(S26 = S26_1 + S26_2 + S26_3) %>%
  mutate(S27 = S27_1 + S27_3) %>%
  mutate(S03 = S03_2 + S03_3) %>%
  mutate(S33 = S33_1 + S33_2 + S33_3) %>%
  mutate(S34 = S34_1 + S34_2 + S34_3) %>%
  mutate(S38 = S38_1 + S38_2 + S38_3) %>%
  mutate(S04 = S04_1 + S04_2 + S04_3) %>%
  mutate(S42 = S42_1 + S42_2 + S42_3) %>%
  mutate(S43 = S43_1 + S43_2 + S43_3) %>%
  mutate(S09 = S09_1  + S09_3) 

OTU_table_total_combined_final_rfy <- OTU_table_total_combined_rfy[,c(54:72)]

#join_taxonomy______________________________________________________________________________________________________________________________________

OTU_table_total_tax_rfy <- OTU_table_total_row_rfy %>% left_join(taxo_table_rfy)
OTU_table_bulk_tax_rfy <- OTU_table_bulk_row_rfy %>% left_join(taxo_table_rfy)
OTU_table_rare_all_tax_rfy <- OTU_table_rare_all_final_rfy %>% left_join(taxo_table_rfy)
OTU_table_total_combined_tax_rfy <- OTU_table_total_combined_final_rfy %>% rownames_to_column("OTU") %>% left_join(taxo_table_rfy)

#safe final OTU table with joined taxonomy
write.csv(OTU_table_total_tax_rfy, "OTU_table_tax_total_rfy.csv")
write.csv(OTU_table_bulk_tax_rfy, "OTU_table_tax_bulk_rfy.csv")
write.csv(OTU_table_rare_all_tax_rfy, "OTU_table_tax_all_rfy.csv")
write.csv(OTU_table_total_combined_tax_rfy, "OTU_table_total_comb_tax_rfy.csv")
