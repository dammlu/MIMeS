# Lukas Damm
# RPTU Kaiserslautern-Landau

library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(tidyverse)
library(data.table)
library(worrms)
library(devtools)
library(testit)

devtools::install_github("janhoo/worms")
devtools::install_github("pmartinezarbizu/dada2pp")

install.packages("remotes")
remotes::install_github("ropensci/worrms")

#load ASV table
ASV_table <- read.csv("ASVTable_nochim.tsv", sep="\t")
head(ASV_table)
colnames(ASV_table)


#load taxonomy table
taxonomy_table_DZMB <- read.csv("Taxonomy_table_DZMB_all.txt", header = TRUE, sep=",")
head(taxonomy_table_DZMB)

#remove duplicate ASVs from taxonomy_table
taxonomy_duplicated_DZMB <- taxonomy_table_DZMB[!duplicated(taxonomy_table_DZMB$asv),]  
write.csv(taxonomy_duplicated_DZMB,"Taxonomy_Table_DZMB_all_final.csv", row.names = FALSE)

taxonomy_final_DZMB <- read.csv("Taxonomy_Table_DZMB_all_final.csv", check.names = TRUE, sep=",")
taxonomy_final_DZMB <- as.data.frame(taxonomy_final_DZMB)

mapping_file <- read.csv("mapping_file_thesis.csv")


# get taxonomy from worms

unique(taxonomy_final_DZMB$assigning)

tax_worms <- wormsbynames(taxonomy_final_DZMB$Species)
tax_worms_curated <- tax_worms[,c(1,3,9,13:18)]
#tax_worms_curated1 <- tax_worms_curated[rowSums(!is.na(tax_worms_curated)) > 0, ]
tax_worms_curated1 <- rename(tax_worms_curated, c("scientificname" = "Species")) 
tax_worms_curated1 <- as.data.frame(tax_worms_curated1)
tax_worms_curated1$Species <- as.character(tax_worms_curated1$Species)
taxonomy_final_DZMB$Species <- as.character(taxonomy_final_DZMB$Species)
# save taxo_worms as csv
write.csv(tax_worms_curated1, "tax_worms_curated.csv")

tax_join <- left_join(taxonomy_final_DZMB, tax_worms_curated1, by = "Species")


tax_join_final <- tax_join[!duplicated(tax_join$asv),]  
write.csv(tax_join_final,"taxonomy_all_worms_joined_final.csv", row.names = FALSE)

taxonomy_join_final <- read.csv("taxonomy_all_worms_joined_final.csv", check.names = TRUE, sep=",")
taxonomy_join_final <- as.data.frame(taxonomy_join_final)

taxonomy_join_final_right <- taxonomy_join_final[,c(1:6)]
taxonomy_join_final_col <- taxonomy_join_final[,c(6,8:13)]
taxonomy_join_final_col1 <- taxonomy_join_final_col[rowSums(!is.na(taxonomy_join_final_col)) > 0, ]
taxonomy_join_final_merge <- left_join(taxonomy_join_final_right, taxonomy_join_final_col1, by = "AphiaID")
taxonomy_join_final_merge_dupli <- taxonomy_join_final_merge[!duplicated(taxonomy_join_final_merge$asv),]  
taxonomy_join_final_merge_NA <-  taxonomy_join_final_merge_dupli %>% filter(rowSums(is.na(.[6:12])) <= 2)

write.csv(taxonomy_join_final_merge_NA,"taxonomy_final_worms_taxo_plots.csv", row.names = FALSE)
write.csv(taxonomy_join_final_merge_dupli, "taxonomy_final_worms_statistic.csv", row.names = FALSE)

taxo_final <- read.csv("taxonomy_all_worms_curated_final.csv",  header = TRUE)

#_remove singletons_____________________________________________________________________________
#transform table column to rows
ASV_table_form <- ASV_table %>% column_to_rownames("asv")


#remove singletons
ASV_table_no_singletons <- ASV_table_form[which(rowSums(ASV_table_form)>1),]
ASV_table_no_singletons_t <- t(ASV_table_no_singletons)

#get read and ASV counts
reads.no.singletons <- colSums(ASV_table_no_singletons)
write.csv(reads.no.singletons, "reads_singletons_removed.csv")

count_positive <- function(row) {return(sum(row > 0))}
count_pos <- apply(ASV_table_no_singletons_t, 1, count_positive)
write.csv(count_pos, "ASVs_singletons_removed.csv")


#filter out samples with less than 10K reads
ASV_table_singletons_filter10 <-ASV_table_no_singletons[, which(colSums(ASV_table_no_singletons)>=10000)]
ASV_table_singletons_filter10_t <- t(ASV_table_singletons_filter10)
samples.out.filter10k <- print(ASV_table_no_singletons[, which(colSums(ASV_table_no_singletons)<=10000)])

#split total and bulk samples
ASV_table_singletons_filter_total <- ASV_table_singletons_filter10[,c(1:53)]
ASV_table_singletons_filter_bulk <- ASV_table_singletons_filter10[,c(54:70)]



#bulk
ASV_table_singletons_sum_bulk <- as.data.frame(colSums(ASV_table_singletons_filter_bulk))
colnames(ASV_table_singletons_sum_bulk)[1]<-"reads"
ASV_table_singletons_sum_bulk <- rownames_to_column(ASV_table_singletons_sum_bulk, var = "sample")

#total
ASV_table_singletons_sum_total <- as.data.frame(colSums(ASV_table_singletons_filter_total))
colnames(ASV_table_singletons_sum_total)[1]<-"reads"
ASV_table_singletons_sum_total <- rownames_to_column(ASV_table_singletons_sum_total, var = "sample")



#load in mapping file
mapping_file_old <- read.csv("mapping_file_old.csv")

#join mapping file and ASV file
ASV_mapped_df_single_bulk <- as.data.frame(ASV_table_singletons_sum_bulk) %>%
  inner_join(mapping_file_old, by = "sample") %>% 
  select(sample,final_sample, new, reads, sampling) %>% 
  filter(sampling !="NA")

ASV_mapped_df_single_total <- as.data.frame(ASV_table_singletons_sum_total) %>%
  inner_join(mapping_file_old, by = "sample") %>% 
  select(sample, final_sample, new, reads, sampling) %>% 
  filter(sampling !="NA")
 


#get min and max number of reads
min(ASV_mapped_df_single_bulk["reads"])
max(ASV_mapped_df_single_bulk["reads"])
ASV_mapped_df_single_row_bulk <- column_to_rownames(ASV_mapped_df_single_bulk, var = "new")
which(ASV_mapped_df_single_bulk["reads"]<=10000)
samples.out.single.bulk <- print(rownames(ASV_mapped_df_single_row_bulk)[which(ASV_mapped_df_single_bulk["reads"]<=10000)])

min(ASV_mapped_df_single_total["reads"])
max(ASV_mapped_df_single_total["reads"])
ASV_mapped_df_single_row_total <- column_to_rownames(ASV_mapped_df_single_total, var = "new")
which(ASV_mapped_df_single_total["reads"]<=10000)
samples.out.single.total <- print(rownames(ASV_mapped_df_single_row_total)[which(ASV_mapped_df_single_total["reads"]<=10000)])

ASV_mapped_df_single_bulk_total <- as.data.frame(t(ASV_mapped_df_single_row_total)) %>% cbind(t(ASV_mapped_df_single_row_bulk)) 
ASV_mapped_df_single_bulk_total_final <- as.data.frame(t(ASV_mapped_df_single_bulk_total)) %>% rownames_to_column("new") %>% mutate_at(vars(reads), as.numeric)


#plot read counts
plot_read_counts_sampling_single <- ggplot(ASV_mapped_df_single_bulk_total_final, aes(x = sampling, y = reads, group = sampling))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("No. reads")+ 
  geom_point(aes(colour = sampling, shape = sampling))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")
plot_read_counts_sampling_single
ggsave("plot_read_counts_single.pdf",plot_read_counts_sampling_single, width = 8, height = 10)  



# plot rarecurves (a la verena) before rarefaction

#join mapping file
ASV_table_singletons_filter_bulk_map <- as.data.frame(t(ASV_table_singletons_filter_bulk)) %>% rownames_to_column("sample") %>%
  inner_join(mapping_file_old, by = "sample") %>% column_to_rownames("final_sample") %>% select(-sample, -sampling, -new, -layer)

ASV_table_singletons_filter_total_map <- as.data.frame(t(ASV_table_singletons_filter_total)) %>% rownames_to_column("sample") %>%
  inner_join(mapping_file_old, by = "sample") %>% column_to_rownames("final_sample") %>% select(-sample, -sampling, -new, -layer)

ASV_table_singletons_filter_bulk_t <- t(ASV_table_singletons_filter_bulk)
ASV_table_singletons_filter_total_t <- t(ASV_table_singletons_filter_total)


pdf(file = "rarecurve_Mangan18_single_before_rfy_bulk.pdf",width=9, height=5)
rarecurve(ASV_table_singletons_filter_bulk_map, step = 50, col = "black", xlab = "No. reads", ylab = "No. of ASVs",  
          label = TRUE, cex = 0.3)
dev.off()


pdf(file = "rarecurve_Mangan18_single_before_rfy_total.pdf",width=9, height=5)
rarecurve(ASV_table_singletons_filter_total_map, step = 50, col = "black", xlab = "No. reads", ylab = "No. of ASVs",  
          label = TRUE, cex = 0.3)
dev.off()



#rarefaction
rowSums(ASV_table_singletons_filter_bulk_map)
min(rowSums(ASV_table_singletons_filter_bulk_map))
rarefaction_level_bulk <- 15005
ASV_table_singletons_rare_bulk <- rrarefy(ASV_table_singletons_filter_bulk_map, rarefaction_level_bulk)


rowSums(ASV_table_singletons_filter_total_map)
min(rowSums(ASV_table_singletons_filter_total_map))
rarefaction_level_total <- 15005
ASV_table_singletons_rare_total <- rrarefy(ASV_table_singletons_filter_total_map, rarefaction_level_total)

rowSums(ASV_table_singletons_filter10_t)
min(rowSums(ASV_table_singletons_filter10_t))
rarefaction_level_all <- 15005
ASV_table_rarefied_total_bulk <- rrarefy(ASV_table_singletons_filter10_t, rarefaction_level_all)
ASV_table_rarefied_total_bulk_t <- as.data.frame(t(ASV_table_rarefied_total_bulk)) %>% rownames_to_column("asv")
write.csv(ASV_table_rarefied_total_bulk_t, "ASV_table_rarefied_total_bulk.csv")


# check if rarefaction was successfull
rowSums(ASV_table_rarefied_total_bulk)
rowSums(ASV_table_singletons_rare_total)



count_positive <- function(row) {return(sum(row > 0))}
count_pos_bulk_rare <- apply(ASV_table_singletons_rare_bulk, 1, count_positive)
count_pos_total_rare <- apply(ASV_table_singletons_rare_total, 1, count_positive)
write.csv(count_pos_bulk_rare, "ASVs_after_rfy_bulk.csv")
write.csv(count_pos_total_rare, "ASVs_after_rfy_total.csv")


count_positive <- function(row) {return(sum(row > 0))}
count_pos_all_rare <- apply(ASV_table_rarefied_total_bulk, 1, count_positive)

write.csv(count_pos_all_rare, "ASVs_after_rarefaction_all.csv")



# plot rarecurves (a la verena) after rarefaction
ASV_table_singletons_rare_bulk_0 <- ASV_table_singletons_rare_bulk[, colSums(ASV_table_singletons_rare_bulk) != 0]
pdf(file = "rarecurve_Mangan18_single_after_rfy_bulk.pdf",width=9, height=5)
rarecurve(ASV_table_singletons_rare_bulk_0, step = 50, col = "black", xlab = "No. reads", ylab = "No. of ASVs",  
          label = TRUE, cex = 0.3)
dev.off()


pdf(file = "rarecurve_Mangan18_single_after_rfy_total.pdf",width=9, height=5)
rarecurve(ASV_table_singletons_rare_total, col = "black", step = 50, xlab = "No. reads", ylab = "No. of ASVs",  
          label = TRUE, cex = 0.3)
dev.off()




# combine columns of layers
ASV_table_total_rare_t <- as.data.frame(t(ASV_table_singletons_rare_total))
ASV_table_total_combined <- ASV_table_total_rare_t %>% 
  mutate(S10 = S10_1 + S10_2 + S10_3) %>%  
  mutate(S11 = S11_1 + S11_2 + S11_3) %>%
  mutate(S13 = S13_1 + S13_2) %>%
  mutate(S16 = S16_1  + S16_2 + S16_3) %>%
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
  mutate(S09 = S09_1 + S09_3) 

ASV_table_total_grouped <- ASV_table_total_combined[,c(54:72)]
colSums(ASV_table_total_grouped)



#plot rarecurves of combined samples after rarefaction
ASV_table_total_grouped_t <- as.data.frame(t(ASV_table_total_grouped))
#filter out ASVs with 0 reads
ASV_table_total_grouped_t_0 <- ASV_table_total_grouped_t[, colSums(ASV_table_total_grouped_t) != 0]
pdf(file = "rarecurve_Mangan18_after_rfy_total_grouped_layers.pdf",width=9, height=5)
rarecurve(ASV_table_total_grouped_t_0, col = "black", step = 50, xlab = "No. reads", ylab = "No. of ASVs",  
          label = TRUE, cex = 0.3)
dev.off()

#combine total and bulk samples
ASV_table_total_grouped_row <- as.data.frame(t(ASV_table_total_grouped_t_0)) %>% rownames_to_column("asv")
ASV_table_bulk_row <- as.data.frame(t(ASV_table_singletons_rare_bulk_0)) %>% rownames_to_column("asv")
ASV_table_rare_all <- full_join(ASV_table_bulk_row,ASV_table_total_grouped_row, by = "asv")
ASV_table_rare_all_final <- ASV_table_rare_all %>%
  mutate_at(vars(-asv), ~ifelse(is.na(.), 0, .))

write.csv(ASV_table_rare_all_final, "ASV_table_rare_all_final.csv")

