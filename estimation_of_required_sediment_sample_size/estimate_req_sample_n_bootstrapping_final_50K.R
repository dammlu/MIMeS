library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)


OTU_table_bulk <- read.csv("bulk_bootstrap_50K.csv")  %>% column_to_rownames("X") %>% t()%>% as.data.frame()

mapping <- read.csv("env_comb.csv")

# Create a vector with file names
file_names <- sprintf("rarefaction_50K_result_%d.csv", 1:19)



# Loop through each file, read it into a data frame, and assign a custom name
for (k in 1:19) {
  file_path <- file.path("/work/RPTU-MIMeS/Paper_eDNA_bulk/OTU_clust_no_rfy/bootstrap_rfy_50K"
                         , file_names[k])
  
  original_file <- read.csv(file_path)
  file_row <- column_to_rownames(original_file, "X")
  final_data <- t(file_row)
  rarefaction_results_rare <- final_data[which(rowSums(final_data)>0),]
  assign(paste0("rarefaction_result_", k), rarefaction_results_rare)
}



# create binary table
for (l in 1:19){
  table_name <- paste0("rarefaction_result_", l)
  tab <- eval(parse(text = table_name))
  data_df <- as.data.frame(tab) %>% rownames_to_column("ASV")
  data_df[-1] <- sapply(data_df[-1], function(x) { as.numeric(x > 0) })
  binary_df <- data_df %>% column_to_rownames("ASV")
  sams <- colnames(binary_df)
  num_samples <- length(sams)
# Create a list of species for each column
  species_lists  <- lapply(binary_df, function(col) rownames(binary_df)[col == 1])
  assign(paste0("species_list_", l), species_lists)
#Combine all species lists
  combined_species_list_all <- unique(unlist(species_lists))
}



# get values from combinations n=4 until n=15 as there are more then 1000 possibilities, but make only 1000 iterations
for (j in 4:15) {
  num_iter <- 1000
  combinations_matrix <- combn(1:num_samples, j)
  size_combinations_matrix <- dim(combinations_matrix)[2]
  res_comb <- matrix(nrow = num_iter, ncol=1)
  for (i in 1:num_iter) {
    # randomly select one of all thousands of combinations  of the combinations matrix
    sams_to_comb <- combinations_matrix[,sample(size_combinations_matrix, 1, replace=TRUE)]
    species_list <- paste0("species_list_", j)
    species_list <- eval(parse(text=species_list))
    species_list_spec <- unique(unlist(species_list[sams[sams_to_comb]]))
    res_comb[i] <- length(species_list_spec)
  write.csv(res_comb, paste0("res_comb", j, ".csv"))
  }
}

# get values from combinations 1,2,3,16,17,18,19 separately (as the have less combinations than 1000 iterations)
nums_smol <- c( 1,2,3,16,17,18,19)
for (j in nums_smol) {
  combinations_matrix <- combn(1:num_samples, j)
  size_combinations_matrix <- dim(combinations_matrix)[2]
  num_iter <- size_combinations_matrix
  res_comb <- matrix(nrow = num_iter, ncol=1)
  for (i in 1:num_iter) {
    sams_to_comb <- combinations_matrix[,i]
    species_list <- paste0("species_list_", j)
    species_list <- eval(parse(text=species_list))
    species_list_spec <- unique(unlist(species_list[sams[sams_to_comb]]))
    res_comb[i] <- length(species_list_spec)
    write.csv(res_comb, paste0("res_comb", j, ".csv"))
  }
}

# read all constructed files to create a bootstrap combination boxplot 
# read all result tables
data_list <- list()
for (i in 1:19) {
  file_name <- paste0("res_comb", i, ".csv")
  data <- read.csv(file_name)
  data_list[[i]] <- data
}

# number of combinations
row_counts <- sapply(data_list, nrow)

# combine all data frames into a single long-format data frame and add leading zeros
combined_data_long <- bind_rows(data_list, .id = "comb")
combined_data_long <- rbind(data.frame(comb = "0", X = 0, V1 = 0), combined_data_long)
combined_data_long$comb <- ifelse(nchar(combined_data_long$comb) == 1,
                                  paste0("0", combined_data_long$comb),
                                  combined_data_long$comb)
combined_data_long_final <- combined_data_long[-1, ]
# boxplot
plot0 <- ggplot(combined_data_long_final, aes(x = comb, y = V1)) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "Distribution of Total OTUs for Different Combination Sizes",
       x = "Combination Size",
       y = "Total Number of OTUs") +
  theme_minimal() +
  theme(legend.position = "none") 
plot0
  ggsave("plot0.pdf", plot0, height=5,width=6)

# 2. fit 3 linear model (linear, mean values)
# Calculate min and max values for each Combination
min_max_values <- combined_data_long %>%
  group_by(comb) %>%
  summarize(min_value = min(V1),
            max_value = max(V1),
            mean_values = mean(V1))

# Calculate linear models for median, min, and max
linear_model_median <-  lm(mean_values ~ as.integer(comb), data = min_max_values)
linear_model_min <- lm(min_value ~ as.integer(comb), data = min_max_values)
linear_model_max <- lm(max_value ~ as.integer(comb), data = min_max_values)



# Create a ggplot with the linear model
plot1 <- ggplot(combined_data_long_final, aes(x = comb, y = V1)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Combination Size",
       y = "Total Number of OTUs") +
  geom_abline(intercept = coef(linear_model_median)[1], slope = coef(linear_model_median)[2], color = "red") +
  geom_abline(intercept = coef(linear_model_min)[1], slope = coef(linear_model_min)[2], color = "grey", lty="dashed") +
  geom_abline(intercept = coef(linear_model_max)[1], slope = coef(linear_model_max)[2], color = "grey", lty="dashed") +
  theme_minimal()+
  coord_cartesian(ylim = c(0, max(combined_data_long$V1) + 100)) +
  geom_smooth(method = 'lm')
plot1
ggsave("plot1.pdf", plot1, height=5,width=6)

# estimate how many samples are needed 
OTU_table_bulk_clean <- OTU_table_bulk[which(rowSums(OTU_table_bulk)>0),]
mean_spec_test <- mean(specnumber(OTU_table_bulk_clean))
target_otus <- mean_spec_test
estimated_samples_needed_med <- (target_otus - coef(linear_model_median)[1]) / coef(linear_model_median)[2]
estimated_samples_needed_max <- (target_otus - coef(linear_model_max)[1]) / coef(linear_model_max)[2]
estimated_samples_needed_min <- (target_otus - coef(linear_model_min)[1]) / coef(linear_model_min)[2]
# one bulk sample equals mean max min x total samples



# calculate a saturation curve model --> is still mostly linear
mean_data <- combined_data_long_final %>% select(-X) %>%
  group_by(comb) %>%
  summarise(mean_Y = mean(V1)) %>%
  mutate(X=1:19) %>% select(-comb)

# fit the Michaelis-Menten model
fit <- nls(mean_Y ~ Vmax * X / (Km + X), data = mean_data, start = list(Vmax = max(mean_data$mean_Y), Km = median(mean_data$X)))
# create a sequence of X values for prediction
x_seq <- seq(min(mean_data$X), max(mean_data$X), length.out = 100)
# predict Y values based on the fitted model
y_pred <- predict(fit, newdata = data.frame(X = x_seq))
mm_line_data <- data.frame(X = x_seq, mean_Y = y_pred)

combined_plot <- plot1 +
  geom_line(data = mm_line_data, aes(x = X, y = mean_Y), col = "blue", size = 0.5) +
  labs(x = "Combination Size",
       y = "Total Number of OTUs", color="s")
combined_plot
ggsave("plot2.pdf", combined_plot, height=5,width=6)

# write a function to estimate the number of samples based on the Michaelis-Menten model
target_otus <- mean_spec_test
estimate_samples <- function(target_otus, model) {
  Vmax <- coef(model)["Vmax"]
  Km <- coef(model)["Km"]
  # Use the rearranged Michaelis-Menten equation
  samples <- (Km * target_otus) / (Vmax - target_otus)
  return(samples)
}

# estimate the number of samples needed 
estimated_samples <- estimate_samples(target_otus, fit)
print(estimated_samples)

# include bulk data to plot
# Create a ggplot with a linear model


# Fit linear model
OTU_tabl_bulk_clean_t <- t(OTU_table_bulk_clean)
specnums_bulk <- as.data.frame(specnumber(t(OTU_tabl_bulk_clean_t))) %>% mutate(comb=25) %>% mutate(X=0) %>% mutate(V1=specnumber(t(OTU_tabl_bulk_clean_t)))
specnums_bulk_c <- specnums_bulk[2:4]
total_bulk <- rbind(combined_data_long_final, specnums_bulk_c)

# Calculate the endpoints of the linear model line so they start at 1 until 19 
x_endpoints <- c(1, 19)
y_endpoints <- coef(linear_model_median)[1] + coef(linear_model_median)[2] * x_endpoints
y_endpoints_min <- coef(linear_model_min)[1] + coef(linear_model_min)[2] * x_endpoints
y_endpoints_max <- coef(linear_model_max)[1] + coef(linear_model_max)[2] * x_endpoints


plot3 <- ggplot(total_bulk, aes(x = comb, y = V1)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Combination Size",
       y = "Total Number of OTUs") +
  geom_line(data = data.frame(X = x_endpoints, V1 = y_endpoints), aes(x = X, y = V1), color = "red")+
  geom_line(data = data.frame(X = x_endpoints, V1 = y_endpoints_min), aes(x = X, y = V1), color = "grey", lty="dashed")+
  geom_line(data = data.frame(X = x_endpoints, V1 = y_endpoints_max), aes(x = X, y = V1), color = "grey", lty="dashed")+
  theme_minimal()+
  coord_cartesian(ylim = c(0, max(combined_data_long$V1) + 100)) +
  scale_x_discrete(labels=paste0(c(1:19, "Bulk")))
plot3

# add michaelis menten model

plot4 <- plot3 +
  geom_line(data = mm_line_data, aes(x = X, y = mean_Y), col = "blue", size = 0.5) +
  labs(x = "Combination Size",
       y = "Total Number of OTUs", color="s")
plot4
ggsave("plot4_50K.pdf", plot4, height=5,width=6)
