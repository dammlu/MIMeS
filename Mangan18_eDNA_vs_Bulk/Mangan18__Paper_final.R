# Lukas Damm
# RPTU Kaiserslautern-Landau


library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(car)
library(ggsignif)
my_lib <- "/path/to/my/lib"


#calculate diversity with OTU table (ASV table rarefied before clustering into OTUs)

#load OTU tables
#save OTU table for diversity measures
OTU_table_bulk_rfy <- read.csv("OTU_table_tax_bulk_rfy.csv",  header = T) %>% column_to_rownames("OTU")
OTU_table_total_rfy<- read.csv("OTU_table_tax_total_rfy.csv", header = T) %>% column_to_rownames("OTU")
OTU_table_all_rfy <- read.csv("OTU_table_tax_all_rfy.csv", header = T) %>% column_to_rownames("OTU")
OTU_table_total_comb_rfy <- read.csv("OTU_table_total_comb_tax_rfy.csv", header = T) %>% column_to_rownames("OTU")


OTU_table_bulk_final_t_rfy <- OTU_table_bulk_rfy[, c(2:18)] %>% t()
OTU_table_total_final_t_rfy <- OTU_table_total_rfy[, c(2:13,15:19,21:40,42:54)] %>% t()
OTU_table_all_final_t_rfy <- OTU_table_all_rfy[, c(2:30,31:36,38:58,60:71)] %>% t()
OTU_table_total_comb_final_t_rfy <- OTU_table_total_comb_rfy[, c(2:20)] %>% t()

colSums(OTU_table_total_rfy$S17_2)

OTU_table_bulk_final_1 <- as.data.frame(t(OTU_table_bulk_final_t_rfy)) %>% rownames_to_column("OTU")
OTU_table_total_comb_1 <- as.data.frame(t(OTU_table_total_comb_final_t_rfy)) %>% rownames_to_column("OTU")
OTU_table_all_comb <- OTU_table_bulk_final_1 %>% left_join(OTU_table_total_comb_1) %>% column_to_rownames("OTU") %>% t() %>% as.data.frame()

OTU_table_bulk_final <- OTU_table_bulk_final_t_rfy[which(rowSums(OTU_table_bulk_final_t_rfy)>0),]
OTU_table_total_final <- OTU_table_total_final_t_rfy[which(rowSums(OTU_table_total_final_t_rfy)>0),]
OTU_table_all_final <- OTU_table_all_final_t_rfy[which(rowSums(OTU_table_all_final_t_rfy)>0),]
OTU_table_total_comb_final <- OTU_table_total_comb_final_t_rfy[which(rowSums(OTU_table_total_comb_final_t_rfy)>0),]
OTU_table_all_comb_final <- OTU_table_all_comb[which(rowSums(OTU_table_all_comb)>0),]

OTU_table_bulk_final <- OTU_table_bulk_final[, colSums(OTU_table_bulk_final) != 0] %>% t() %>% as.data.frame()
OTU_table_total_final <- OTU_table_total_final[, colSums(OTU_table_total_final) != 0] %>% t() %>% as.data.frame()
OTU_table_all_final <- OTU_table_all_final[, colSums(OTU_table_all_final) != 0] %>% t() %>% as.data.frame()
OTU_table_total_comb_final <- OTU_table_total_comb_final[, colSums(OTU_table_total_comb_final) != 0] %>% t() %>% as.data.frame()
OTU_table_all_comb_final <- OTU_table_all_comb_final[, colSums(OTU_table_all_comb_final) !=0] %>% t()


write.csv(OTU_table_bulk_final, "OTU_table_bulk_final.csv")
write.csv(OTU_table_total_final, "OTU_table_total_final.csv")
write.csv(OTU_table_all_final, "OTU_table_all_final.csv")
write.csv(OTU_table_total_comb_final, "OTU_table_total_comb_final.csv")
write.csv(OTU_table_all_comb_final, "OTU_table_all_comb_final.csv")


#get ASV counts after rarefaction and before OTU clustering

ASV_table_post_rarefy <- read.csv("ASV_table_rarefied_total_bulk.csv")

ASV_post_rarefy_total <- ASV_table_post_rarefy[, c(2:55)] %>% column_to_rownames("asv")
ASV_post_rarefy_bulk <- ASV_table_post_rarefy[, c(2,56:72)] %>% column_to_rownames("asv")

ASV_post_rfy_total_final <- ASV_post_rarefy_total[, colSums(ASV_post_rarefy_total) !=0]
ASV_post_rfy_total_final <- ASV_post_rfy_total_final[which(rowSums(ASV_post_rfy_total_final)>0),]

ASV_post_rfy_bulk_final <- ASV_post_rarefy_bulk[, colSums(ASV_post_rarefy_bulk) !=0]
ASV_post_rfy_bulk_final <- ASV_post_rfy_bulk_final[which(rowSums(ASV_post_rfy_bulk_final)>0),]

# load mapping file 
env_rfy <- read.csv("env1.csv")
env_comb <- read.csv("env_comb.csv")

__________________________________________________________________________________________________________________________________________________________________________________

#### alpha diversity calculations ####


#calculate shannon index
shannon_bulk_rfy <- diversity(OTU_table_bulk_final, index="shannon")
shannon_df_bulk_rfy <- as.data.frame(shannon_bulk_rfy)

shannon_total_rfy <- diversity(OTU_table_total_final, index="shannon")
shannon_df_total_rfy <- as.data.frame(shannon_total_rfy)

shannon_total_comb_rfy <- diversity(OTU_table_total_comb_final, index="shannon")
shannon_df_total_comb_rfy <- as.data.frame(shannon_total_comb_rfy)

#calculate simpson index
simpson_bulk_rfy <- diversity(OTU_table_bulk_final, index="simpson")
simpson_df_bulk_rfy <- as.data.frame(simpson_bulk_rfy)

simpson_total_rfy <- diversity(OTU_table_total_final, index="simpson")
simpson_df_total_rfy <- as.data.frame(simpson_total_rfy)

simpson_total_comb_rfy <- diversity(OTU_table_total_comb_final, index="simpson")
simpson_df_total_comb_rfy <- as.data.frame(simpson_total_comb_rfy)

#calculate otu richness
otu_richness_bulk_rfy <- specnumber(OTU_table_bulk_final)
otu_richness_df_bulk_rfy <- as.data.frame(otu_richness_bulk_rfy)

otu_richness_total_rfy <- specnumber(OTU_table_total_final)
otu_richness_df_total_rfy <- as.data.frame(otu_richness_total_rfy)

otu_richness_total_comb_rfy <- specnumber(OTU_table_total_comb_final)
otu_richness_df_total_comb_rfy <- as.data.frame(otu_richness_total_comb_rfy)

OTU_div_stats_bulk_rfy <- as.data.frame(shannon_df_bulk_rfy) %>% cbind(simpson_df_bulk_rfy) %>% cbind(otu_richness_df_bulk_rfy) %>% rownames_to_column("sample")
OTU_div_stats_total_rfy <- as.data.frame(shannon_df_total_rfy) %>% cbind(simpson_df_total_rfy) %>% cbind(otu_richness_df_total_rfy) %>% rownames_to_column("sample")
OTU_div_stats_final_rfy <- as.data.frame(t(OTU_div_stats_bulk_rfy)) %>% cbind(t(OTU_div_stats_total_rfy))
OTU_div_stats_final_t_rfy <- t(OTU_div_stats_final_rfy)
rownames(OTU_div_stats_final_t_rfy) <- NULL
colnames(OTU_div_stats_final_t_rfy) <- c("sample", "shannon", "simpson", "otu_richness")
OTU_div_stats_final_data_rfy <- as.data.frame(OTU_div_stats_final_t_rfy)
OTU_div_stats_final_data_rfy <- OTU_div_stats_final_data_rfy %>% column_to_rownames("sample")

OTU_div_stats_total_comb_rfy <- as.data.frame(shannon_df_total_comb_rfy) %>% cbind(simpson_df_total_comb_rfy) %>% cbind(otu_richness_df_total_comb_rfy) %>% rownames_to_column("sample")
OTU_div_stats_final_comb_rfy <- as.data.frame(t(OTU_div_stats_bulk_rfy)) %>% cbind(t(OTU_div_stats_total_comb_rfy))
OTU_div_stats_final_comb_t_rfy <- t(OTU_div_stats_final_comb_rfy)
rownames(OTU_div_stats_final_comb_t_rfy) <- NULL
colnames(OTU_div_stats_final_comb_t_rfy) <- c("sample", "shannon", "simpson", "otu_richness")
OTU_div_stats_final_data_comb_rfy <- as.data.frame(OTU_div_stats_final_comb_t_rfy)
OTU_div_stats_final_data_comb_rfy <- OTU_div_stats_final_data_comb_rfy %>% column_to_rownames("sample")


OTU_div_map_final_data_rfy <- as.data.frame(OTU_div_stats_final_data_rfy) %>% rownames_to_column("sample") %>%
  inner_join(env_rfy, by = "sample") %>% 
  select(sample, shannon, simpson, otu_richness, sampling, layer) 
OTU_div_map_final_data_rfy <- OTU_div_map_final_data_rfy %>% mutate_at(vars(shannon, simpson, otu_richness), as.numeric)

OTU_div_map_final_data_comb_rfy <- as.data.frame(OTU_div_stats_final_data_comb_rfy) %>% rownames_to_column("sample") %>%
  inner_join(env_comb, by = "sample") %>% 
  select(sample, shannon, simpson, otu_richness, sampling) 
OTU_div_map_final_data_comb_rfy <- OTU_div_map_final_data_comb_rfy %>% mutate_at(vars(shannon, simpson, otu_richness), as.numeric)

#plot alpha diversity

library(RColorBrewer)
palette2 <- brewer.pal(3, "Set2")[1:2]

plot_div_1_rfy <- ggplot(OTU_div_map_final_data_rfy, aes(x=final_sample, y=otu_richness, fill=sampling))+
  geom_bar(stat="identity")+
  facet_grid(. ~ sampling, scales="free", space = "free_x")+
  theme_light()+
  scale_fill_manual(values = c(palette2))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 8),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_line(size = 0.5, color="black"), 
        strip.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), legend.position = "none")+ labs(x = "Sample", y = "rarefied OTU richness")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
plot_div_1_rfy
ggsave("barplot_otu_richness_final.pdf", plot= plot_div_1_rfy,width = 8, height = 7)


plot_div_2_rfy <- ggplot(OTU_div_map_final_data_rfy, aes(x=final_sample, y=simpson, fill=sampling))+
  geom_bar(stat="identity")+
  facet_grid(. ~ sampling, scales="free", space = "free_x")+
  theme_light()+
  scale_fill_manual(values = c(palette2))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 8),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_line(size = 0.5, color="black"), 
        strip.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), legend.position = "none")+ labs(x = "Sample", y = "Simpson index value", size = 12)+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
plot_div_2_rfy
ggsave("barplot_div_simpson_rfy.pdf", plot= plot_div_2_rfy,width = 8, height = 7)


plot_div_3_rfy <- ggplot(OTU_div_map_final_data_rfy, aes(x=final_sample, y=shannon, fill=sampling))+
  geom_bar(stat="identity")+
  facet_grid(. ~ sampling, scales="free", space = "free_x")+
  theme_light()+
  scale_fill_manual(values = c(palette2))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 8),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_line(size = 0.5, color="black"), 
        strip.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), legend.position = "none",)+ labs(x = "Sample", y = "Shannon index value", size = 12)+ 
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
plot_div_3_rfy
ggsave("barplot_div_shannon_rfy.pdf", plot= plot_div_3_rfy,width = 8, height = 7)

#plot alpha diversity of total combined layers and bulk of rarefied OTU table

plot_div_1_comb_rfy <- ggplot(OTU_div_map_final_data_comb_rfy, aes(x=final_sample, y=otu_richness, fill=sampling))+
  geom_bar(stat="identity")+
  facet_grid(. ~ sampling, scales="free", space = "free_x")+
  theme_light()+
  scale_fill_manual(values = c(palette2))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 10),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_line(size = 0.5, color="black"), 
        strip.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), legend.position = "none")+ labs(x = "Sample", y = "rarefied OTU richness")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
plot_div_1_comb_rfy
ggsave("barplot_otu_richness_comb_rfy.pdf", plot= plot_div_1_comb_rfy,width = 5, height = 6)


plot_div_2_comb_rfy <- ggplot(OTU_div_map_final_data_comb_rfy, aes(x=final_sample, y=simpson, fill=sampling))+
  geom_bar(stat="identity")+
  facet_grid(. ~ sampling, scales="free", space = "free_x")+
  theme_light()+
  scale_fill_manual(values = c(palette2))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 10),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_line(size = 0.5, color="black"), 
        strip.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), legend.position = "none")+ labs(x = "Sample", y = "Simpson index value", size = 12)+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
plot_div_2_comb_rfy
ggsave("barplot_div_simpson_comb_rfy.pdf", plot= plot_div_2_comb_rfy,width = 5, height = 6)


plot_div_3_comb_rfy <- ggplot(OTU_div_map_final_data_comb_rfy, aes(x=final_sample, y=shannon, fill=sampling))+
  geom_bar(stat="identity")+
  facet_grid(. ~ sampling, scales="free", space = "free_x")+
  theme_light()+
  scale_fill_manual(values = c(palette2))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 10),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_line(size = 0.5, color="black"), 
        strip.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12))+
  theme(panel.grid = element_blank(), legend.position = "none",)+ labs(x = "Sample", y = "Shannon index value", size = 12)+ 
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
plot_div_3_comb_rfy
ggsave("barplot_div_shannon_comb_rfy.pdf", plot= plot_div_3_comb_rfy,width = 5, height = 6)

# boxplot alpha div of separate SeDNA layers

OTU_div_map_final_data_rfy_sep <- OTU_div_map_final_data_rfy[18:67, ]

palette3 <- brewer.pal(3, "Set3")

box_plot_shannon_rfy <- ggplot(OTU_div_map_final_data_rfy_sep, aes(x = layer, y = shannon, group = layer))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("Shannon index value")+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")+
  geom_signif(comparisons = list(c("0-1 cm","1-2 cm")),
              map_signif_level = TRUE, y_position = 3.0, textsize = 3)+
  geom_signif(comparisons = list(c("0-1 cm","2-3 cm")),
              map_signif_level = TRUE, y_position = 3.2, textsize = 3)+
  geom_signif(comparisons = list(c("1-2 cm","2-3 cm")),
              map_signif_level = TRUE, y_position = 3.5, textsize = 3)
box_plot_shannon_rfy
ggsave("boxplot_layers_shannon_rfy.pdf",box_plot_shannon_rfy, width = 5, height = 7)


box_plot_simpson_rfy <- ggplot(OTU_div_map_final_data_rfy_sep, aes(x = layer, y = simpson, group = layer))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("Simpson index value")+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")+
  geom_signif(comparisons = list(c("0-1 cm","1-2 cm")),
            map_signif_level = TRUE, y_position = 1.1, textsize = 3)+
  geom_signif(comparisons = list(c("0-1 cm","2-3 cm")),
              map_signif_level = TRUE, y_position = 1.2, textsize = 3)+
  geom_signif(comparisons = list(c("1-2 cm","2-3 cm")),
              map_signif_level = TRUE, y_position = 1.3, textsize = 3)
box_plot_simpson_rfy
ggsave("boxplot_layers_simpson_rfy.pdf",box_plot_simpson_rfy, width = 5, height = 7)


box_plot_richness_rfy <- ggplot(OTU_div_map_final_data_rfy_sep, aes(x = layer, y = otu_richness, group = layer))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("rarefied OTU richness")+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")+
  geom_signif(comparisons = list(c("0-1 cm","1-2 cm")),
              map_signif_level = TRUE, y_position = 190, textsize = 3)+
  geom_signif(comparisons = list(c("0-1 cm","2-3 cm")),
              map_signif_level = TRUE, y_position = 220, textsize = 3)+
  geom_signif(comparisons = list(c("1-2 cm","2-3 cm")),
              map_signif_level = TRUE, y_position = 240, textsize = 3)
box_plot_richness_rfy
ggsave("boxplot_layers_richness_rfy.pdf",box_plot_richness_rfy, width = 5, height = 7)

alpha_div_layers <- (box_plot_richness_rfy + box_plot_shannon_rfy + box_plot_simpson_rfy) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold"))
alpha_div_layers
ggsave("alpha_div_layers.pdf",alpha_div_layers, width = 10, height = 8)


# statistic for alpha div measures ####

#test for sign. diff. between layers

OTU_div_map_final_data_rfy$layer <- as.factor(OTU_div_map_final_data_rfy$layer)

sediment_layers_data_rfy <- OTU_div_map_final_data_rfy %>% filter(layer %in% c("0-1 cm","1-2 cm","2-3 cm"))

normality_tests_layers_rfy <- sediment_layers_data_rfy %>% 
  group_by(layer) %>%
  summarise(shannon_p_value = shapiro.test(shannon)$p.value,
            simspon_p_value = shapiro.test(simpson)$p.value,
            otu_richness_p_value = shapiro.test(otu_richness)$p.value)


levene_shannon_layers_rfy <- leveneTest(shannon ~ layer, data = sediment_layers_data_rfy)
levene_simpson_layers_rfy <- leveneTest(simpson ~ layer, data = sediment_layers_data_rfy)
levene_otu_richness_layers_rfy <- leveneTest(otu_richness ~ layer, data = sediment_layers_data_rfy)

kruskal_shannon_layers_rfy <- kruskal.test(shannon ~ layer, data = sediment_layers_data_rfy)
kruskal_simpson_layers_rfy <- kruskal.test(simpson ~ layer, data = sediment_layers_data_rfy)
kruskal_otu_richness_layers_rfy <- kruskal.test(otu_richness ~ layer, data = sediment_layers_data_rfy)

#test for sign. diff. between layers including bulk

normality_tests_sampling_rfy <- OTU_div_map_final_data_rfy %>% 
  group_by(layer) %>%
  summarise(shannon_p_value = shapiro.test(shannon)$p.value,
            simspon_p_value = shapiro.test(simpson)$p.value,
            otu_richness_p_value = shapiro.test(otu_richness)$p.value)


levene_shannon_sampling_rfy <- leveneTest(shannon ~ layer, data = OTU_div_map_final_data_rfy)
levene_simpson_sampling_rfy <- leveneTest(simpson ~ layer, data = OTU_div_map_final_data_rfy)
levene_otu_richness_sampling_rfy <- leveneTest(otu_richness ~ layer, data = OTU_div_map_final_data_rfy)

kruskal_shannon_sampling_rfy <- kruskal.test(shannon ~ layer, data = OTU_div_map_final_data_rfy)
kruskal_simpson_sampling_rfy <- kruskal.test(simpson ~ layer, data = OTU_div_map_final_data_rfy)
kruskal_otu_richness_sampling_rfy <- kruskal.test(otu_richness ~ layer, data = OTU_div_map_final_data_rfy)

#post-hoc test for sign. Kruskal-Wallis

pairwise_wilcox_shannon_rfy <- pairwise.wilcox.test(OTU_div_map_final_data_rfy$shannon, OTU_div_map_final_data_rfy$layer, p.adjust.method = "BH")
print(pairwise_wilcox_shannon_rfy)
pairwise_wilcox_otu_richness_rfy <- pairwise.wilcox.test(OTU_div_map_final_data_rfy$otu_richness, OTU_div_map_final_data_rfy$layer, p.adjust.method = "BH")
print(pairwise_wilcox_otu_richness_rfy)
pairwise_wilcox_simpson_rfy <- pairwise.wilcox.test(OTU_div_map_final_data_rfy$simpson, OTU_div_map_final_data_rfy$layer, p.adjust.method = "BH")
print(pairwise_wilcox_simpson_rfy)

# boxplot alpha div combined layers rarefied on ASV level


box_plot_shannon_comb_rfy <- ggplot(OTU_div_map_final_data_comb_rfy, aes(x = sampling, y = shannon, group = sampling))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("Shannon index value")+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")+
  geom_signif(comparisons = list(c("SeDNA","Bulk")),
              map_signif_level = TRUE, y_position = 4.5, textsize = 3)
box_plot_shannon_comb_rfy
ggsave("boxplot_layers_shannon_comb_rfy.pdf",box_plot_shannon_comb_rfy, width = 5, height = 7)


box_plot_simpson_comb_rfy <- ggplot(OTU_div_map_final_data_comb_rfy, aes(x = sampling, y = simpson, group = sampling))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("Simpson index value")+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")+
  geom_signif(comparisons = list(c("SeDNA","Bulk")),
              map_signif_level = TRUE, y_position = 1.0, textsize = 3)
box_plot_simpson_comb_rfy
ggsave("boxplot_layers_simpson_comb_rfy.pdf",box_plot_simpson_comb_rfy, width = 5, height = 7)


box_plot_richness_comb_rfy <- ggplot(OTU_div_map_final_data_comb_rfy, aes(x = sampling, y = otu_richness, group = sampling))+ 
  geom_boxplot(fill = "grey")+
  theme_bw()+
  ylab("rarefied OTU richness")+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18), axis.title.x = element_blank())+
  theme(legend.position = "none")+
  geom_signif(comparisons = list(c("SeDNA","Bulk")),
              map_signif_level = TRUE, y_position = 550, textsize = 3)
box_plot_richness_comb_rfy
ggsave("boxplot_layers_richness_comb_rfy.pdf",box_plot_richness_comb_rfy, width = 5, height = 7)


##### add species accumulation curve

OTU_table_bulk_final_t <- OTU_table_bulk_final %>% t() %>% as.data.frame
OTU_table_total_comb_final_t <- OTU_table_total_comb_final

specaccum_bulk <- specaccum(OTU_table_bulk_final_t, method = "random")
specaccum_total <- specaccum(OTU_table_total_comb_final_t, method = "random")

df_bulk_spec <- data.frame(Sites = specaccum_bulk$sites, Richness = specaccum_bulk$richness, SD = specaccum_bulk$sd, Method = "Bulk")
df_total_spec <- data.frame(Sites = specaccum_total$sites, Richness = specaccum_total$richness, SD = specaccum_total$sd, Method = "SeDNA")

df_combined <- rbind(df_bulk_spec, df_total_spec)
label_positions <- df_combined %>%
  group_by(Method) %>%
  slice_tail(n = 1)

species_accum_plot <- ggplot(df_combined, aes(x = Sites, y = Richness, color = Method, fill = Method)) +
  geom_line(size = 1.2, color = "black") +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD, fill = Method, color = Method), alpha = 0.2) +
  scale_fill_manual(values = c("Method1" = "lightgrey", "Method2" = "gainsboro")) +
  scale_color_manual(values = c("Method1" = "darkgrey", "Method2" = "grey40")) +
  labs(x = "Number of Samples", y = "OTU Richness", size = 18) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))+
  geom_text_repel(data = label_positions,
                  aes(label = Method),
                  color = "black",
                  nudge_x = 5,
                  direction = "y",
                  hjust = 0,
                  segment.color = NA,
                  size = 6) +
  theme(legend.position = "none")
species_accum_plot
ggsave("species_accum_plot.pdf",species_accum_plot, width = 5, height = 7)

#combine alpha div plots 

combined_alpha_div <- (box_plot_richness_comb_rfy + box_plot_shannon_comb_rfy + box_plot_simpson_comb_rfy + species_accum_plot) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_alpha_div
ggsave("combined_alpha_div.pdf",combined_alpha_div, width = 10, height = 14)

___________________________________________________________________________________________________________________________________________________________

#### beta diversity ####
library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(pracma)
library(ggrepel)
library(stringr)
library(cluster)
library(reshape2)

#load mapping_file
env_rfy <- read.csv("env1.csv") 


# NMDS

# get table without taxonomic annotation

tab_notax_bulk_rfy <- OTU_table_bulk_final 
tab_notax_total_rfy <- OTU_table_total_final
tab_notax_total_out_rfy <- tab_notax_total_rfy[, c(1:4,6:50)] #S11_2 out

tab_notax_all_rfy <- read.csv("OTU_table_all_final.csv") %>% column_to_rownames("X")
tab_notax_all_out_rfy <- tab_notax_all_rfy[, c(1:21,23:66)]


tab_notax_all_comb <- read.csv("OTU_table_all_comb_final.csv") %>% column_to_rownames("X")

#total
set.seed(1)
tab_root_total_rfy <- nthroot((t(tab_notax_total_out_rfy)), 4)
otu_rare_total <-tab_root_total_rfy[,which(colSums(tab_root_total_rfy)>0)]
otu_rare_total <-tab_root_total_rfy[which(rowSums(tab_root_total_rfy)>0),]

set.seed(666)
bray_total_sep_rfy <- vegdist(otu_rare_total, method = 'bray')
nmds_FM_total_sep_rfy <-metaMDS(bray_total_sep_rfy, k=2)
saveRDS(nmds_FM_total_sep_rfy, "nmds_FM_total_sep_rfy.rds")
nmds_FM_total_sep_rfy <- readRDS("nmds_FM_total_sep_rfy.rds") 
nmds_df_total_sep_rfy<-as.data.frame(nmds_FM_total_sep_rfy$points) %>%
  rownames_to_column("sample_new") %>%
  left_join(env_rfy) 


#hull data for total separate layers

hull.data_rfy <- nmds_df_total_sep_rfy %>% 
  group_by(layer) %>%
  summarize(NMDS1 = MDS1[chull(MDS1,MDS2)],
            NMDS2 = MDS2[chull(MDS1,MDS2)])

nmds_total_sep_rfy <-ggplot()+
  geom_point(data = nmds_df_total_sep_rfy,aes(MDS1,MDS2, color = layer))+
  geom_polygon(data = hull.data_rfy, aes(x = NMDS1, y = NMDS2, fill = layer, group = layer), alpha = 0.30) + # add the convex hulls
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x = "NMDS1", y = "NMDS2")+
  geom_text_repel(data = nmds_df_total_sep_rfy,aes(MDS1,MDS2, label = sample))+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_total_sep_rfy$stress, digits=4))) +
  scale_fill_discrete(name = "Layer") + scale_color_discrete(name = "Layer") + scale_shape_discrete(name = "Layer")
nmds_total_sep_rfy
ggsave("nmds_total_sep_layers_rfy_hulls.pdf",nmds_total_sep_rfy,height=8,width=9)

#try with conf intervals
nmds_total_sep_rfy_conf <-ggplot()+
  geom_point(data = nmds_df_total_sep_rfy,aes(MDS1,MDS2, color = layer))+
  stat_ellipse(data = nmds_df_total_sep_rfy, aes(MDS1,MDS2, color = layer, fill = layer), level = 0.95, geom = "polygon", alpha = 0.2)+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x = "NMDS1", y = "NMDS2")+
  #geom_text_repel(data = nmds_df_total_sep_rfy,aes(MDS1,MDS2, label = sample))+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_total_sep_rfy$stress, digits=4))) +
  scale_fill_discrete(name = "Layer") + scale_color_discrete(name = "Layer") + scale_shape_discrete(name = "Layer")
nmds_total_sep_rfy_conf
ggsave("nmds_total_sep_layers_rfy_conf.pdf",nmds_total_sep_rfy_conf,height=6,width=7)


#bulk and total NMDS
set.seed(1)
tab_root_all_rfy <- nthroot((t(tab_notax_all_out_rfy)), 4)
otu_rare_all <-tab_root_all_rfy[,which(colSums(tab_root_all_rfy)>0)]
otu_rare_all <-tab_root_all_rfy[which(rowSums(tab_root_all_rfy)>0),]

set.seed(666)
bray_total_bulk_rfy <- vegdist(otu_rare_all, method = 'bray')
nmds_FM_all_rfy  <- metaMDS(bray_total_bulk_rfy, k=2)
saveRDS(nmds_FM_all_rfy , "nmds_FM_all_rfy.rds")
nmds_FM_all_rfy  <- readRDS("nmds_FM_all_rfy.rds") 
nmds_df_FM_all_rfy  <-as.data.frame(nmds_FM_all_rfy$points) %>%
  rownames_to_column("final_sample") %>%
  left_join(env_rfy) 

hull.data_rfy  <- nmds_df_FM_all_rfy  %>%
  group_by(extraction) %>%
  summarize(NMDS1 = MDS1[chull(MDS1, MDS2)],
            NMDS2 = MDS2[chull(MDS1, MDS2)])


nmds_all_rfy   <-ggplot()+
  geom_point(data = nmds_df_FM_all_rfy, aes(MDS1,MDS2, color = extraction))+
  
  geom_polygon(data = hull.data_rfy  , aes(x = NMDS1, y = NMDS2, fill = extraction, group = extraction), alpha = 0.30) +
  #geom_text_repel(data = nmds_df_FM_all_rfy  ,aes(MDS1,MDS2, label = final_sample), max.overlaps = 20)+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(fill = "Extraction method", shape = "Extraction method")+
  labs(x = "NMDS1", y = "NMDS2")+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_all_rfy$stress, digits=4)))+
  scale_fill_manual( values = palette2)
nmds_all_rfy
ggsave("nmds_total_bulk_layers_rfy_bray_4throot.pdf",nmds_all_rfy ,height=8,width=9)


#now with 95% conf intervals

nmds_all_rfy_conf   <-ggplot()+
  geom_point(data = nmds_df_FM_all_rfy  ,aes(MDS1,MDS2, color = extraction))+
  stat_ellipse(data = nmds_df_FM_all_rfy, aes(MDS1,MDS2, color = extraction, fill = extraction), level = 0.95, geom = "polygon", alpha = 0.2)+
  #geom_text_repel(data = nmds_df_FM_all_rfy  ,aes(MDS1,MDS2, label = final_sample), max.overlaps = 20)+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(fill = "Extraction method", shape = "Extraction method", color = "Extraction method")+
  labs(x = "NMDS1", y = "NMDS2")+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_all_rfy$stress, digits=4)))+
  scale_fill_manual( values = palette2) + scale_colour_manual( values = palette2)
nmds_all_rfy_conf
ggsave("nmds_total_bulk_layers_rfy_bray_4throot_conf.pdf",nmds_all_rfy_conf ,height=6,width=7)



#bulk and total NMDS combined layers
set.seed(1)
tab_root_all_comb <- nthroot((t(tab_notax_all_comb)), 4)
otu_rare_all_comb <-tab_root_all_comb[,which(colSums(tab_root_all_comb)>0)]
otu_rare_all_comb <-tab_root_all_comb[which(rowSums(tab_root_all_comb)>0),]

set.seed(666)
bray_total_bulk_comb <- vegdist(otu_rare_all_comb, method = 'bray')
nmds_FM_all_comb  <- metaMDS(bray_total_bulk_comb, k=2)
saveRDS(nmds_FM_all_comb , "nmds_FM_all_comb.rds")
nmds_FM_all_comb  <- readRDS("nmds_FM_all_comb.rds") 
nmds_df_FM_all_comb  <-as.data.frame(nmds_FM_all_comb$points) %>%
  rownames_to_column("final_sample") %>%
  left_join(env_comb) 

hull.data_comb  <- nmds_df_FM_all_comb %>%
  group_by(extraction) %>%
  summarize(NMDS1 = MDS1[chull(MDS1, MDS2)],
            NMDS2 = MDS2[chull(MDS1, MDS2)])


nmds_all_comb   <-ggplot()+
  geom_point(data = nmds_df_FM_all_comb  ,aes(MDS1,MDS2,))+
  
  geom_polygon(data = hull.data_comb  , aes(x = NMDS1, y = NMDS2, fill = extraction, group = extraction), alpha = 0.30) +
  geom_text_repel(data = nmds_df_FM_all_comb  ,aes(MDS1,MDS2, label = final_sample), max.overlaps = 20)+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(fill = "Sampling", shape = "Sampling")+
  labs(x = "NMDS1", y = "NMDS2")+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_all_comb$stress, digits=4)))+
  scale_fill_manual( values = palette2)
nmds_all_comb
ggsave("nmds_total_bulk_comb_bray_4throot.pdf",nmds_all_comb ,height=5,width=6)

#also with 95% conf intervals

nmds_all_comb_conf   <-ggplot()+
  geom_point(data = nmds_df_FM_all_comb  ,aes(MDS1,MDS2, color = extraction))+
  stat_ellipse(data = nmds_df_FM_all_comb, aes(MDS1,MDS2, color = extraction, fill = extraction), level = 0.95, geom = "polygon", alpha = 0.2)+
  #geom_text_repel(data = nmds_df_FM_all_comb  ,aes(MDS1,MDS2, label = final_sample), max.overlaps = 20)+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(fill = "Extraction method", shape = "Extraction method", color = "Extraction method")+
  labs(x = "NMDS1", y = "NMDS2")+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_all_comb$stress, digits=4)))+
  scale_fill_manual( values = palette2) + scale_colour_manual( values = palette2)+
  theme(legend.position = "none")
nmds_all_comb_conf
ggsave("nmds_total_bulk_comb_bray_4throot_conf.pdf",nmds_all_comb_conf,height=5,width=6)

#nmds based on presence/absence

#tab notax all combine

otu_rare_all_comb_binary <- ifelse(tab_notax_all_comb > 0,1,0) %>% t()

set.seed(888)
jaccard_total_bulk_comb_binary <- vegdist(otu_rare_all_comb_binary, method = 'jaccard', binary = T)
nmds_FM_all_comb_binary  <- metaMDS(jaccard_total_bulk_comb_binary, k=2)
saveRDS(nmds_FM_all_comb_binary , "nmds_FM_all_comb_binary.rds")
nmds_FM_all_comb_binary  <- readRDS("nmds_FM_all_comb_binary.rds") 
nmds_df_FM_all_comb_binary  <-as.data.frame(nmds_FM_all_comb_binary$points) %>%
  rownames_to_column("final_sample") %>%
  left_join(env_comb) 

hull.data_comb_binary  <- nmds_df_FM_all_comb_binary %>%
  group_by(extraction) %>%
  summarize(NMDS1 = MDS1[chull(MDS1, MDS2)],
            NMDS2 = MDS2[chull(MDS1, MDS2)])


nmds_all_comb_binary  <- ggplot()+
  geom_point(data = nmds_df_FM_all_comb_binary,aes(MDS1,MDS2, color = extraction, fill = extraction))+
  stat_ellipse(data = nmds_df_FM_all_comb_binary, aes(MDS1,MDS2, group = extraction, color = extraction, fill = extraction), level = 0.95, geom = "polygon", alpha = 0.2)+
  theme_light()+
  theme(panel.grid = element_blank(), legend.title = element_text(size = 16),  # Legend title size
        legend.text = element_text(size = 14))+
  labs(fill = "Extraction method", color = "Extraction method")+
  labs(x = "NMDS1", y = "NMDS2")+
  annotate("text", Inf, Inf,hjust=+1,vjust=+1.2,
           label=paste0("stress = ",round(nmds_FM_all_comb_binary$stress, digits=4)))+
  scale_fill_manual(values = palette2)+
    scale_color_manual(values = palette2)
nmds_all_comb_binary
ggsave("nmds_total_bulk_comb_jaccard.pdf", nmds_all_comb_binary,height=5,width=6)



#combine beta div plots
combined_beta_div <- (nmds_all_comb_conf + nmds_all_comb_binary) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_beta_div
ggsave("combined_beta_div.pdf",combined_beta_div, width = 12, height = 8)


______________________________________________________________________________________________________________________________________________________________________________

#get number of shared OTUs bulk and total

count_positive_rfy <- function(row) {return(sum(row > 0))}
count_pos_bulk_rare_rfy <- apply(tab_notax_bulk_rfy, 1, count_positive_rfy)
count_pos_total_rare_rfy <- apply(tab_notax_total_out_rfy, 1, count_positive_rfy)
write.csv(count_pos_bulk_rare_rfy, "OTUs_bulk_rfy.csv")
write.csv(count_pos_total_rare_rfy, "OTUs_total_rfy.csv")


#draw pairwise venn
install.packages("VennDiagram", lib = my_lib)
library(VennDiagram, lib.loc = my_lib)



#venn for rarefied OTUs
otus_bulk_rfy <- as.data.frame(count_pos_bulk_rare_rfy) %>% rownames_to_column("OTU")
otus_total_rfy <- as.data.frame(count_pos_total_rare_rfy) %>% rownames_to_column("OTU")

rows_bulk_rfy <- otus_bulk_rfy[otus_bulk_rfy$count_pos_bulk_rare_rfy >= 1, ] %>% as.data.frame() 
rows_total_rfy <- otus_total_rfy[otus_total_rfy$count_pos_total_rare_rfy >= 1, ] %>% as.data.frame()
samples_venn_bulk_rfy <- rows_bulk_rfy[, c(1)] 
samples_venn_total_rfy <- rows_total_rfy[, c(1)] 

pdf("venn_plot_rfy1.pdf", width = 11, height = 9)
grid.newpage()
venn.plot.rfy1 <- draw.pairwise.venn(length(samples_venn_bulk_rfy),
                                    length(samples_venn_total_rfy),
                                    length(intersect(samples_venn_bulk_rfy, samples_venn_total_rfy)),
                                    fill = (palette2),
                                    alpha = 0.5,
                                    category = c("Bulk", "SeDNA"), scaled = F,
                                    cat.dist = c(0.04, 0.04),
                                    cat.cex = 1.8,
                                    cex = 2.5)

grid.text(
  label = paste0(round(percent_bulk_only_bulk, 1), "% reads"),
  x = 0.2, y = 0.40,  # adjust these x,y coords to fit nicely inside left circle
  gp = gpar(fontsize = 16, col = "black")
)
grid.text(
  label = paste0(round(percent_total_only_total, 1), "% reads"),
  x = 0.8, y = 0.40,  # adjust as needed for right circle
  gp = gpar(fontsize = 16, col = "black")
)
grid.text(
  label = paste0("Bulk: ", round(percent_shared_bulk, 1), "% reads"),
  x = 0.50, y = 0.40,
  gp = gpar(fontsize = 16, col = "black")
)
grid.text(
  label = paste0("SeDNA: ", round(percent_shared_total, 1), "% reads"),
  x = 0.50, y = 0.35,
  gp = gpar(fontsize = 16, col = "black")
)
dev.off()
grid.draw(venn.plot.rfy1)   

#calculate relative abundance of core OTUs and shared
otus_bulk_only <- setdiff(samples_venn_bulk_rfy, samples_venn_total_rfy)
otus_total_only <- setdiff(samples_venn_total_rfy, samples_venn_bulk_rfy)
otus_shared <- intersect(samples_venn_bulk_rfy, samples_venn_total_rfy)
bulk_only_tab <- tab_notax_bulk_rfy %>% rownames_to_column("OTU") %>% filter(OTU %in% otus_bulk_only)
total_only_tab <- tab_notax_total_rfy %>% rownames_to_column("OTU") %>% filter(OTU %in% otus_total_only)

shared_in_bulk <- tab_notax_bulk_rfy %>% rownames_to_column("OTU") %>% filter(OTU %in% otus_shared)
shared_in_total <- tab_notax_total_rfy %>% rownames_to_column("OTU") %>% filter(OTU %in% otus_shared)

total_reads_bulk <- sum(as.matrix(tab_notax_bulk_rfy))
total_reads_total <- sum(as.matrix(tab_notax_total_rfy))

bulk_only_reads_bulk <- sum(as.matrix(tab_notax_bulk_rfy[rownames(tab_notax_bulk_rfy) %in% otus_bulk_only, ]))
total_only_reads_total <- sum(as.matrix(tab_notax_total_rfy[rownames(tab_notax_total_rfy) %in% otus_total_only, ]))
shared_reads_bulk <- sum(as.matrix(tab_notax_bulk_rfy[rownames(tab_notax_bulk_rfy) %in% otus_shared, ]))
shared_reads_total <- sum(as.matrix(tab_notax_total_rfy[rownames(tab_notax_total_rfy) %in% otus_shared, ]))
percent_bulk_only_bulk <- (bulk_only_reads_bulk / total_reads_bulk) * 100
percent_total_only_total <- (total_only_reads_total / total_reads_total) * 100
percent_shared_bulk <- (shared_reads_bulk / total_reads_bulk) * 100
percent_shared_total <- (shared_reads_total / total_reads_total) * 100




#calculate overlap of OTUs across layers for the rarefied dataset

OTU_layer_one_rfy <- OTU_table_total_rfy[, c(2,5,8,10,13,16,19,22,25,28,31,35,38,41,44,47,50)]
OTU_layer_two_rfy <- OTU_table_total_rfy[, c(3,6,9,11,14,17,20,23,26,29,33,36,39,42,45,48)]
OTU_layer_three_rfy <- OTU_table_total_rfy[, c(4,7,12,15,18,21,24,27,30,32,34,37,40,43,46,49)]
OTU_bulk_venn_rfy <- OTU_table_bulk_rfy[, c(2:18)]


OTU_layer_one_final_rfy <- OTU_layer_one_rfy[which(rowSums(OTU_layer_one_rfy)>= 1) ,] %>% rownames_to_column("OTU")
OTU_layer_two_final_rfy <- OTU_layer_two_rfy[which(rowSums(OTU_layer_two_rfy)>= 1) ,] %>% rownames_to_column("OTU")
OTU_layer_three_final_rfy <- OTU_layer_three_rfy[which(rowSums(OTU_layer_three_rfy)>= 1) ,] %>% rownames_to_column("OTU")
OTU_bulk_venn_final_rfy <- OTU_bulk_venn_rfy[which(rowSums(OTU_bulk_venn_rfy )>= 1) ,] %>% rownames_to_column("OTU")
OTU_table_layer_comb_venn <- OTU_table_total_comb_final[which(rowSums(OTU_table_total_comb_final)>= 1) ,] %>% rownames_to_column("OTU")
write.csv(OTU_layer_one_final_rfy, "OTU_layer_one_final_rfy.csv")
write.csv(OTU_layer_two_final_rfy, "OTU_layer_two_final_rfy.csv")
write.csv(OTU_layer_three_final_rfy, "OTU_layer_three_final_rfy.csv")
write.csv(OTU_bulk_venn_final_rfy, "OTU_bulk_venn_rfy.csv")
write.csv(OTU_table_layer_comb_venn, "OTU_layer_comb_venn_rfy.csv")

____________________________________________________________________________________________________________________________________________________________________

#### taxonomy plots ####

library(RColorBrewer)


taxo_class_bulk <- OTU_table_bulk_rfy %>% rownames_to_column("OTU") %>% select(1,22) %>% sub("\\|.*", "", taxo_class_bulk$Class) %>% as.data.frame()
taxo_class_total <- OTU_table_total_rfy %>% rownames_to_column("OTU") %>% select(1,58) %>% as.data.frame()
write.csv(taxo_class_bulk, "taxo_class_bulk.csv")
write.csv(taxo_class_total, "taxo_class_total.csv")

taxo_class_bulk <- read.csv("taxo_class_bulk.csv") %>% select(2,3)
taxo_class_total <- read.csv("taxo_class_total.csv") %>% select(2,3)

bulk_taxonomy_class <- OTU_table_bulk_final %>% rownames_to_column("OTU") %>% left_join(taxo_class_bulk)
total_taxonomy_class <- OTU_table_total_final %>% rownames_to_column("OTU") %>% left_join(taxo_class_total)

total_taxonomy_class_comb <- OTU_table_total_comb_final %>% rownames_to_column("OTU") %>% left_join(taxo_class_total)

# include colors for taxo plot class

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
env_colors <- c(cbPalette[3], cbPalette[4],cbPalette[5],cbPalette[2],cbPalette[7])

##colors for taxo plot
colPalette2 <- rev(brewer.pal(12, "Paired"))
colPalette3 <- colPalette2[-6]
colPalette4 <- brewer.pal(10, "Spectral")
palette5 <- c(colPalette2, colPalette4)
palette5[16] <- "#808080"
pal_a <- rev(brewer.pal(8, "Accent"))
palette5 <- c(colPalette2, pal_a)
palette5[19] <- "#808080"
palette5[21] <- "#000000"
palette5[17] <- "#40e0d0"
palette5 <- c(palette5, palette5)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector2 <- col_vector
col_vector2[31] <- "#333333"



#get taxo plot class
all_classes <- unique(c(class_df_bulk_in_single$Class30, class_df_total_in_single$Class30))
color_class <- setNames(col_vector2, all_classes)



OTU_table_class_bulk <- bulk_taxonomy_class[,c(2:19)] %>% 
  group_by(Class) %>% 
  summarize_each(funs(sum)) %>% 
  as.data.frame() %>%
  filter(Class!="") %>%
  column_to_rownames("Class") %>% 
  mutate(total=rowSums(.)) %>%
  filter(total>0)
OTU_table_class_bulk <- OTU_table_class_bulk[order(-OTU_table_class_bulk$total), ]
top30_bulk <- rownames(OTU_table_class_bulk)[1:30]
classes_per_sample <- apply(OTU_table_class_bulk > 0, 2, sum)



# no single samples
# additionally, all not top 30 have still to be included!
OTU_table_class_bulk <- OTU_table_class_bulk  %>% decostand(MARGIN=2, method = "total")#%>% select(-total)
write.csv(OTU_table_class_bulk, "relab_OTU_table_class_bulk.csv")
class_bulk_df_single <- OTU_table_class_bulk %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  gather("Class", "count", -sample) %>%
  mutate(Class30=ifelse(Class %in% top30_bulk, Class, "zOther")) %>%
  left_join(env_comb)
class_df_bulk_in_NOTHER_single <- class_bulk_df_single %>% select(Class30, count, sampling, final_sample) %>% group_by(Class30, count, sampling, final_sample) %>% summarise_each(funs(mean)) %>%
  filter(Class30!="zOther")
class_df_bulk_in_OTHER_single <- class_bulk_df_single %>% select(Class30, count, sampling, final_sample) %>% group_by(Class30, count, sampling, final_sample) %>% summarise_each(funs(sum)) %>%
  filter(Class30=="zOther")
class_df_bulk_in_single <- rbind(class_df_bulk_in_NOTHER_single,class_df_bulk_in_OTHER_single)
class_df_bulk_in_single$count <- as.numeric(class_df_bulk_in_single$count)

plot_taxo_class_bulk_single<-ggplot(class_df_bulk_in_single, aes(x=final_sample, y=count, fill=Class30)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  labs(fill="Class", x="Sample", y="Relative abundance of assigned reads [%]")+
  scale_fill_manual(values=color_class)
plot_taxo_class_bulk_single
ggsave("taxo_class_bulk.pdf",plot_taxo_class_bulk_single, width = 10, height = 8)


OTU_table_class_total <- total_taxonomy_class[2:52] %>% 
  group_by(Class) %>% 
  summarize_each(funs(sum)) %>% 
  as.data.frame() %>%
  filter(Class!="") %>%
  column_to_rownames("Class") %>% 
  mutate(total=rowSums(.))%>%
  filter(total>0)
OTU_table_class_total <- OTU_table_class_total[order(-OTU_table_class_total$total), ]
top30_total <- rownames(OTU_table_class_total)[1:30]
classes_per_sample <- apply(OTU_table_class_total > 0, 2, sum)


OTU_table_class_total <- OTU_table_class_total %>% select(-total) %>% decostand(MARGIN=2, method = "total")
write.csv(OTU_table_class_total, "relab_OTU_table_class_total.csv")
class_total_df_single <- OTU_table_class_total %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  gather("Class", "count", -sample) %>%
  mutate(Class30=ifelse(Class %in% top30_total, Class, "z.Other")) %>%
  left_join(env_rfy)
class_df_total_in_NOTHER_single <- class_total_df_single %>% select(Class30, count, sampling, sample, layer) %>% group_by(Class30, count, sampling, sample, layer) %>% summarise_each(funs(mean)) %>%
  filter(Class30!="z.Other")
class_df_total_in_OTHER_single <- class_total_df_single %>% select(Class30, count, sampling, sample, layer) %>% group_by(Class30, count, sampling, sample, layer) %>% summarise_each(funs(sum)) %>%
  filter(Class30=="z.Other")
class_df_total_in_single <- rbind(class_df_total_in_NOTHER_single,class_df_total_in_OTHER_single)
class_df_total_in_single$count <- as.numeric(class_df_total_in_single$count)

plot_taxo_class_total_single<-ggplot(class_df_total_in_single, aes(x=sample, y=count, fill=Class30)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
  labs(fill="Class", x="Sample", y="Relative abundance of assigned reads [%]")+
  facet_grid(. ~ layer, scales = "free_x")+
  scale_fill_manual(values=color_class)
plot_taxo_class_total_single
ggsave("taxo_class_total_layers.pdf",plot_taxo_class_total_single, width = 10, height = 8)

#taxo plot combined layers

OTU_table_class_total_comb <- total_taxonomy_class_comb[2:21] %>% 
  group_by(Class) %>% 
  summarize_each(funs(sum)) %>% 
  as.data.frame() %>%
  filter(Class!="") %>%
  column_to_rownames("Class") %>% 
  mutate(total=rowSums(.))%>%
  filter(total>0)
OTU_table_class_total_comb <- OTU_table_class_total_comb[order(-OTU_table_class_total_comb$total), ]
top30_total_comb <- rownames(OTU_table_class_total_comb)[1:30]
classes_per_sample_comb <- apply(OTU_table_class_total > 0, 2, sum)


OTU_table_class_total_comb <- OTU_table_class_total_comb  %>% decostand(MARGIN=2, method = "total")#%>% select(-total)
write.csv(OTU_table_class_total_comb, "relab_OTU_table_class_total.csv")
class_total_df_single_comb <- OTU_table_class_total_comb %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  gather("Class", "count", -sample) %>%
  mutate(Class30=ifelse(Class %in% top30_total_comb, Class, "z.Other")) %>%
  left_join(env_comb)
class_df_total_in_NOTHER_single_comb <- class_total_df_single_comb %>% select(Class30, count, sampling, sample) %>% group_by(Class30, count, sampling, sample) %>% summarise_each(funs(mean)) %>%
  filter(Class30!="z.Other")
class_df_total_in_OTHER_single_comb <- class_total_df_single_comb %>% select(Class30, count, sampling, sample) %>% group_by(Class30, count, sampling, sample) %>% summarise_each(funs(sum)) %>%
  filter(Class30=="z.Other")
class_df_total_in_single_comb <- rbind(class_df_total_in_NOTHER_single_comb,class_df_total_in_OTHER_single_comb)
class_df_total_in_single_comb$count <- as.numeric(class_df_total_in_single_comb$count)

plot_taxo_class_total_single_comb<-ggplot(class_df_total_in_single_comb, aes(x=sample, y=count, fill=Class30)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  labs(fill="Class", x="Sample", y="Relative abundance of assigned reads [%]")+
  scale_fill_manual(values=color_class)
plot_taxo_class_total_single_comb
ggsave("taxo_class_total_comb.pdf",plot_taxo_class_total_single_comb, width = 10, height = 8)


#get taxo plot phyla

taxo_phyla_bulk <- OTU_table_bulk_rfy %>% rownames_to_column("OTU") %>% select(1,21) %>% as.data.frame() #%>% sub("\\|.*", "", taxo_phyla_bulk$Phylum) 
taxo_phyla_total <- OTU_table_total_rfy %>% rownames_to_column("OTU") %>% select(1,57) %>% as.data.frame()
write.csv(taxo_phyla_bulk, "taxo_phyla_bulk.csv")
write.csv(taxo_phyla_total, "taxo_phyla_total.csv")

taxo_phyla_bulk <- read.csv("taxo_phyla_bulk.csv") %>% select(2,3)
taxo_phyla_total <- read.csv("taxo_phyla_total.csv") %>% select(2,3)

bulk_taxonomy_phyla <- OTU_table_bulk_final %>% rownames_to_column("OTU") %>% left_join(taxo_phyla_bulk)
total_taxonomy_phyla <- OTU_table_total_final %>% rownames_to_column("OTU") %>% left_join(taxo_phyla_total)


total_taxonomy_phyla_comb <- OTU_table_total_comb_final %>% rownames_to_column("OTU") %>% left_join(taxo_phyla_total)

#get taxo plot phyla
all_phyla <- unique(c(phyla_df_bulk_in_single$Phyla30, phyla_df_total_in_single$Phyla30))
color_phyla <- setNames(col_vector2, all_phyla)



OTU_table_phyla_bulk <- bulk_taxonomy_phyla[,c(2:19)] %>% 
  group_by(Phylum) %>% 
  summarize_each(funs(sum)) %>% 
  as.data.frame() %>%
  filter(Phylum!="") %>%
  column_to_rownames("Phylum") %>% 
  mutate(total=rowSums(.)) %>%
  filter(total>0)
OTU_table_phyla_bulk <- OTU_table_phyla_bulk[order(-OTU_table_phyla_bulk$total), ]
top30_bulk <- rownames(OTU_table_phyla_bulk)[1:30]
phyla_per_sample <- apply(OTU_table_phyla_bulk > 0, 2, sum)



# no single samples
# additionally, all not top 30 have still to be included!
OTU_table_phyla_bulk <- OTU_table_phyla_bulk %>% select(-total) %>% decostand(MARGIN=2, method = "total")
write.csv(OTU_table_phyla_bulk, "relab_OTU_table_phyla_bulk.csv")
phyla_bulk_df_single <- OTU_table_phyla_bulk %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  gather("Phylum", "count", -sample) %>%
  mutate(Phyla30=ifelse(Phylum %in% top30_bulk, Phylum, "zOther")) %>%
  left_join(env_comb)
phyla_df_bulk_in_NOTHER_single <- phyla_bulk_df_single %>% select(Phyla30, count, extraction, final_sample) %>% group_by(Phyla30, count, extraction, final_sample) %>% summarise_each(funs(mean)) %>%
  filter(Phyla30!="zOther")
phyla_df_bulk_in_OTHER_single <- phyla_bulk_df_single %>% select(Phyla30, count, extraction, final_sample) %>% group_by(Phyla30, count, extraction, final_sample) %>% summarise_each(funs(sum)) %>%
  filter(Phyla30=="zOther")
phyla_df_bulk_in_single <- rbind(phyla_df_bulk_in_NOTHER_single,phyla_df_bulk_in_OTHER_single)
phyla_df_bulk_in_single$count <- as.numeric(phyla_df_bulk_in_single$count)

plot_taxo_phyla_bulk_single<-ggplot(phyla_df_bulk_in_single, aes(x=final_sample, y=count, fill=Phyla30)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  labs(fill="Phylum", x="Sample", y="Relative abundance of assigned reads [%]")+
  scale_fill_manual(values=color_phyla)
plot_taxo_phyla_bulk_single
ggsave("taxo_phyla_bulk.pdf",plot_taxo_phyla_bulk_single, width = 10, height = 8)


OTU_table_phyla_total <- total_taxonomy_phyla[2:52] %>% 
  group_by(Phylum) %>% 
  summarize_each(funs(sum)) %>% 
  as.data.frame() %>%
  filter(Phylum!="") %>%
  column_to_rownames("Phylum") %>% 
  mutate(total=rowSums(.))%>%
  filter(total>0)
OTU_table_phyla_total <- OTU_table_phyla_total[order(-OTU_table_phyla_total$total), ]
top30_total <- rownames(OTU_table_phyla_total)[1:30]
phyla_per_sample <- apply(OTU_table_phyla_total > 0, 2, sum)


OTU_table_phyla_total <- OTU_table_phyla_total %>% select(-total) %>% decostand(MARGIN=2, method = "total")
write.csv(OTU_table_phyla_total, "relab_OTU_table_phyla_total.csv")
phyla_total_df_single <- OTU_table_phyla_total %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  gather("Phylum", "count", -sample) %>%
  mutate(Phyla30=ifelse(Phylum %in% top30_total, Phylum, "z.Other")) %>%
  left_join(env_rfy)
phyla_df_total_in_NOTHER_single <- phyla_total_df_single %>% select(Phyla30, count, sampling, sample, layer) %>% group_by(Phyla30, count, sampling, sample, layer) %>% summarise_each(funs(mean)) %>%
  filter(Phyla30!="z.Other")
phyla_df_total_in_OTHER_single <- phyla_total_df_single %>% select(Phyla30, count, sampling, sample, layer) %>% group_by(Phyla30, count, sampling, sample, layer) %>% summarise_each(funs(sum)) %>%
  filter(Phyla30=="z.Other")
phyla_df_total_in_single <- rbind(phyla_df_total_in_NOTHER_single,phyla_df_total_in_OTHER_single)
phyla_df_total_in_single$count <- as.numeric(phyla_df_total_in_single$count)

plot_taxo_phyla_total_single<-ggplot(phyla_df_total_in_single, aes(x=sample, y=count, fill=Phyla30)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
  labs(fill="Phylum", x="Sample", y="Relative abundance of assigned reads [%]")+
  facet_grid(. ~ layer, scales = "free_x")+
  scale_fill_manual(values=color_phyla)
plot_taxo_phyla_total_single
ggsave("taxo_phyla_total_layers.pdf",plot_taxo_phyla_total_single, width = 10, height = 8)

#taxo plot combined layers

OTU_table_phyla_total_comb <- total_taxonomy_phyla_comb[2:21] %>% 
  group_by(Phylum) %>% 
  summarize_each(funs(sum)) %>% 
  as.data.frame() %>%
  filter(Phylum!="") %>%
  column_to_rownames("Phylum") %>% 
  mutate(total=rowSums(.))%>%
  filter(total>0)
OTU_table_phyla_total_comb <- OTU_table_phyla_total_comb[order(-OTU_table_phyla_total_comb$total), ]
top30_total_comb <- rownames(OTU_table_phyla_total_comb)[1:30]
phyla_per_sample_comb <- apply(OTU_table_phyla_total > 0, 2, sum)


OTU_table_phyla_total_comb <- OTU_table_phyla_total_comb  %>% select(-total) %>% decostand(MARGIN=2, method = "total")
write.csv(OTU_table_phyla_total_comb, "relab_OTU_table_phyla_total.csv")
phyla_total_df_single_comb <- OTU_table_phyla_total_comb %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  gather("Phylum", "count", -sample) %>%
  mutate(Phyla30=ifelse(Phylum %in% top30_total_comb, Phylum, "z.Other")) %>%
  left_join(env_comb)
phyla_df_total_in_NOTHER_single_comb <- phyla_total_df_single_comb %>% select(Phyla30, count, extraction, sample) %>% group_by(Phyla30, count, extraction, sample) %>% summarise_each(funs(mean)) %>%
  filter(Phyla30!="z.Other")
phyla_df_total_in_OTHER_single_comb <- phyla_total_df_single_comb %>% select(Phyla30, count, extraction, sample) %>% group_by(Phyla30, count, extraction, sample) %>% summarise_each(funs(sum)) %>%
  filter(Phyla30=="z.Other")
phyla_df_total_in_single_comb <- rbind(phyla_df_total_in_NOTHER_single_comb,phyla_df_total_in_OTHER_single_comb)
phyla_df_total_in_single_comb$count <- as.numeric(phyla_df_total_in_single_comb$count)

plot_taxo_phyla_total_single_comb<-ggplot(phyla_df_total_in_single_comb, aes(x=sample, y=count, fill=Phyla30)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  labs(fill="Phylum", x="Sample", y="Relative abundance of assigned reads [%]")+
  scale_fill_manual(values=color_phyla)+
  theme(legend.position = "none")
plot_taxo_phyla_total_single_comb
ggsave("taxo_phyla_total_comb.pdf",plot_taxo_phyla_total_single_comb, width = 10, height = 8)

#combine plots

combined_plot_phyla_samples <- (plot_taxo_phyla_total_single_comb + plot_taxo_phyla_bulk_single) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_plot_phyla_samples
ggsave("combined_phyla_plot_samples.pdf", combined_plot_phyla_samples, width = 10, height = 8)


#### get phyla plot but summarized per extraction method #####

phyla_total_long <- total_taxonomy_phyla_comb[, c(2:21)] %>%
  pivot_longer(cols = starts_with("S"),  
               names_to = "Sample", 
               values_to = "Abundance_tot")

phyla_summary_total <- phyla_total_long %>%
  group_by(Phylum) %>%
  summarise(total_abundance_tot = sum(Abundance_tot), .groups = "drop") %>%
  mutate(rel_abundance_tot = total_abundance_tot / sum(total_abundance_tot))

phyla_bulk_long <- bulk_taxonomy_phyla[, c(2:19)] %>%
  pivot_longer(cols = starts_with("SMO"),  
               names_to = "Sample", 
               values_to = "Abundance_bulk")

phyla_summary_bulk <- phyla_bulk_long %>%
  group_by(Phylum) %>%
  summarise(total_abundance_bulk = sum(Abundance_bulk ), .groups = "drop") %>%
  mutate(rel_abundance_bulk  = total_abundance_bulk  / sum(total_abundance_bulk ))

phyla_summary_combined <- phyla_summary_bulk %>% left_join(phyla_summary_total) %>%
  mutate_all(~replace(., is.na(.), 0)) 

phyla_summary_combined_long <- phyla_summary_combined %>%
  pivot_longer(cols = starts_with("rel_abundance_bulk") | starts_with("rel_abundance_tot"),
               names_to = "Type",
               values_to = "Abundance") %>%
  mutate(Type = ifelse(grepl("bulk", Type), "Bulk", "SeDNA"))

plot_taxo_phyla_extraction <-ggplot(phyla_summary_combined_long, aes(x= Type, y= Abundance, fill= Phylum)) +
  geom_bar(stat="identity", position="fill") +
  theme_light()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 14), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  labs(fill="Phylum", x="Extraction method", y="Relative abundance of assigned reads [%]")+
  scale_fill_manual(values=color_phyla)
plot_taxo_phyla_extraction
ggsave("taxo_phyla_extraction_method.pdf",plot_taxo_phyla_extraction, width = 10, height = 8)

#### get stacked bar plot but use otu richness instead of relative read abundance

# Convert to presence/absence 
total_taxonomy_pa <- total_taxonomy_phyla_comb
total_taxonomy_pa[, 2:(ncol(total_taxonomy_pa)-1)] <- ifelse(total_taxonomy_pa[, 2:(ncol(total_taxonomy_pa)-1)] > 0, 1, 0)


# Bulk
bulk_taxonomy_pa <- bulk_taxonomy_phyla
bulk_taxonomy_pa[, 2:(ncol(bulk_taxonomy_pa)-1)] <- ifelse(bulk_taxonomy_pa[, 2:(ncol(bulk_taxonomy_pa)-1)] > 0, 1, 0)

# SeDNA OTU richness per phylum
phyla_total_otu_richness <- total_taxonomy_pa %>%
  group_by(Phylum) %>%
  summarise(otu_richness_tot = n(), .groups = "drop")

# Bulk OTU richness per phylum
phyla_bulk_otu_richness <- bulk_taxonomy_pa %>%
  group_by(Phylum) %>%
  summarise(otu_richness_bulk = n(), .groups = "drop")

#combine and calculate relative OTU richness
phyla_richness_combined <- phyla_bulk_otu_richness %>% 
  left_join(phyla_total_otu_richness, by = "Phylum") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  mutate(rel_otu_richness_bulk = otu_richness_bulk / sum(otu_richness_bulk),
         rel_otu_richness_tot = otu_richness_tot / sum(otu_richness_tot))


phyla_richness_combined_long <- phyla_richness_combined %>%
  pivot_longer(cols = starts_with("rel_otu_richness"), 
               names_to = "Type", 
               values_to = "Richness") %>%
  mutate(Type = ifelse(grepl("bulk", Type), "Bulk", "SeDNA"))

plot_otu_richness_phyla_extraction <- ggplot(phyla_richness_combined_long, 
                                             aes(x=Type, y=Richness, fill=Phylum)) +
  geom_bar(stat="identity", position="fill") +
  theme_light() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=14), 
        axis.text.y = element_text(size=10), 
        axis.title = element_text(size=18),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=12)) +
  labs(fill="Phylum", x="Extraction method", y="Relative OTU richness [%]") +
  scale_fill_manual(values=color_phyla)

plot_otu_richness_phyla_extraction
ggsave("otu_richness_phyla_extraction_method.pdf", plot_otu_richness_phyla_extraction, width = 10, height = 8)

library(patchwork)

combined_plot_phyla <- (plot_taxo_phyla_extraction + plot_otu_richness_phyla_extraction) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a") & theme(legend.position = 'right', plot.tag = element_text(size = 18, face = "bold"))
combined_plot_phyla
ggsave("combined_phyla_extraction_method.pdf", combined_plot_phyla, width = 10, height = 8)

#### get otu and read counts ####

OTU_table_all_final_t <- t(OTU_table_all_final)
reads_all <- rowSums(OTU_table_all_final_t) %>% as.data.frame()
colnames(reads_all) <- "reads"

count_positive_otu <- function(row) {return(sum(row > 0))}
OTU_counts <- apply(OTU_table_all_final_t, 1, count_positive_otu) %>% as.data.frame()
colnames(OTU_counts) <- "OTUs"


#get taxo information

phyla_in_bulk <- OTU_table_phyla_bulk %>% rownames_to_column("Phylum") 
phyla_in_bulk <- phyla_in_bulk[, c(1,19)]
phyla_in_bulk <- column_to_rownames(phyla_in_bulk, "Phylum") 
phyla_in_bulk_perc <- phyla_in_bulk %>% mutate(colsum_total = sum(total), percentage = sprintf("%.2f", (total / colsum_total) * 100))
phyla_in_bulk_final <- phyla_in_bulk_perc %>% rownames_to_column("Phylum")
write.csv(phyla_in_bulk_final, "phyla_in_bulk_final.csv")

phyla_in_total <- OTU_table_phyla_total_comb %>% rownames_to_column("Phylum") 
phyla_in_total <- phyla_in_total[, c(1,21)]
phyla_in_total <- column_to_rownames(phyla_in_total, "Phylum") 
phyla_in_total_perc <- phyla_in_total %>% mutate(colsum_total = sum(total), percentage = sprintf("%.2f", (total / colsum_total) * 100))
phyla_in_total_final <- phyla_in_total_perc %>% rownames_to_column("Phylum")
write.csv(phyla_in_total_final, "phyla_in_total_final.csv")


classes_in_bulk <- OTU_table_class_bulk %>% rownames_to_column("Class") 
classes_in_bulk <- classes_in_bulk[, c(1,19)]
classes_in_bulk <- column_to_rownames(classes_in_bulk, "Class") 
classes_in_bulk_perc <- classes_in_bulk %>% mutate(colsum_total = sum(total), percentage = sprintf("%.2f", (total / colsum_total) * 100))
classes_in_bulk_final <- classes_in_bulk_perc %>% rownames_to_column("Class")
write.csv(classes_in_bulk_final, "classes_in_bulk_final.csv")

classes_in_total <- OTU_table_class_total_comb %>% rownames_to_column("Class") 
classes_in_total <- classes_in_total[, c(1,21)]
classes_in_total <- column_to_rownames(classes_in_total, "Class") 
classes_in_total_perc <- classes_in_total %>% mutate(colsum_total = sum(total), percentage = sprintf("%.2f", (total / colsum_total) * 100))
classes_in_total_final <- classes_in_total_perc %>% rownames_to_column("Class")
write.csv(classes_in_total_final, "classes_in_total_final.csv")

#order level
taxo_order_total <- OTU_table_total_comb_rfy %>% 
  select(2:20,24) %>%
  filter(rowSums(.[, 1:19]) > 0) %>%
  as.data.frame() %>%
  rownames_to_column("OTU") 
write.csv(taxo_order_total, "taxo_order_total.csv")

taxo_order_bulk <- OTU_table_bulk_rfy %>% rownames_to_column("OTU") %>% select(1,23) %>% as.data.frame()
write.csv(taxo_order_bulk, "taxo_order_bulk.csv")

taxo_order_bulk <- OTU_table_bulk_rfy %>% 
  select(2:18,22) %>%
  filter(rowSums(.[, 1:17]) > 0) %>%
  as.data.frame() %>%
  rownames_to_column("OTU") 
write.csv(taxo_order_bulk, "taxo_order_bulk.csv")
