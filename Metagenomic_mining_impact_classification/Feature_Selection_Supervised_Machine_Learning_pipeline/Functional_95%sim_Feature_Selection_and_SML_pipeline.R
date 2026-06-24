# M.Sc. Lukas Damm
# RPTU Kaiserslautern-Landau
# 07.03.2025

# curate metagenomic count tables of assembled functional genes and reconstructed/assembled taxonomic genes (16S) 
# eDNA extracted from sediment samples (MUC&PUC) colelcted in the Clarion-Clipperton-Fracture Zone in the Pacific 
# samples are from the BGR area (SO295 cruise), 1.5 years post polymetallic nodule collector test
# count table is merged kallisto output (tpm)

# use CLR-transformed count tables as input for Feature Selection (FS) and Supervised Machine Learning (SML) pipeline. 
# evaluate four different FS methods (Boruta, RFE, PIMP, r2VIM) and use selected features to classify mining impact categories (Collector impact, Plume impact and undisturbed Reference)
# using Random Forest (RF) and Support Vector Machine (SVM) based on microbial signatures (functional genes)



#load packages

my_lib <- "/work/RPTU-MIMeS/R_from_ldamm/R/4.3"
library(tidyverse)
library(dplyr)
library(compositions)
library(vegan)
library(caret)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(openxlsx, lib.loc = my_lib)
library(vegan)
library(pairwiseAdonis, lib.loc = my_lib)
library(cowplot,        lib.loc = my_lib)
library(tibble)
library(e1071)
library(vita, lib.loc = my_lib)
library(kernlab, lib.loc = my_lib)
library(randomForest)
library(MLmetrics, lib.loc = my_lib)
library(Boruta,     lib.loc = my_lib)
library(groupdata2, lib.loc = my_lib)
library(Pomona,     lib.loc = my_lib)

#===================================================================================
# load count tables
count_table_func <- read_tsv("merged_tpm_matrix_BGR_95.tsv")
count_table_func <- count_table_func %>% column_to_rownames("target_id")
count_table_func_t <- t(count_table_func) %>% as.data.frame()

#===================================================================================
# calc overall sparsity
sparsity <- sum(count_table_func_t == 0) / prod(dim(count_table_func_t))
cat("Overall sparsity:", round(sparsity * 100, 2), "%\n")

# per-feature sparsity 
feature_sparsity <- colMeans(count_table_func_t == 0)
summary(feature_sparsity)
hist(feature_sparsity, breaks = 50, 
     main = "Per-feature sparsity", 
     xlab = "Proportion of zero samples")

# per-sample sparsity
sample_sparsity <- rowMeans(count_table_func_t == 0)
summary(sample_sparsity)


#remove rows where rowSums is zero
count_table_func_zero <- count_table_func[rowSums(count_table_func) !=0, ]



#add small pseudocount for following CLR trans
pseudocount <- 0.00001

count_table_func_final <- count_table_func_zero + pseudocount

#===================================================================================
#perform CLR transformation

count_table_func_clr <- clr(count_table_func_final)

write.table(count_table_func_clr, file = "clr_transformed_count_table_func_BGR_95sim.tsv", sep = "\t", quote = F, row.names = T, col.names = NA)

clr_trans_table_func <- read_tsv("clr_transformed_count_table_func_BGR_95sim.tsv")

colnames(clr_trans_table_func)[1] <- "target_id"

clr_trans_table_func <- clr_trans_table_func %>% column_to_rownames("target_id") %>% t() %>% as.data.frame()

#load mapping file and merge with count table
map <- read.csv("mapping_DNA_BGR.csv")

site <- map[, c(1,3,7)]

site$station <- as.character(site$station)
site$station <- as.factor(site$station)

data <- clr_trans_table_func %>% rownames_to_column("sample") %>% inner_join(site) %>% column_to_rownames("sample")
data[, -ncol(data)] <- lapply(data[, -ncol(data)], as.numeric)
data$site <- make.names((data$site))
data$site <- as.factor(data$site)
str(data)
write.csv(data, "data_table_BGR_site_func_95sim.csv")

data[is.na(data)] <- 0.00001

#make sure $site is a factor

data$site <- as.factor(data$site)
data$station <- as.factor(as.character(data$station))
data$sample <- rownames(data)


#===================================================================================
#define functions that are needed throughout the FS+SML 

n_iter <- 50
FALLBACK_N <- 20

my_ctrl <- trainControl(
  method          = "cv",
  number          = 10,
  classProbs      = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary)


extract_metrics <- function(cm, iteration, method, n_features) {
  
  overall    <- cm$overall
  byclass    <- cm$byClass
  f1_scores  <- byclass[, "F1"]
  mean_f1    <- mean(f1_scores, na.rm = TRUE)
  class_names <- rownames(byclass)
  
  df <- data.frame(
    Iteration  = iteration,
    Method     = method,
    n_features = n_features,
    Accuracy   = overall["Accuracy"],
    Kappa      = overall["Kappa"],
    Mean_F1    = mean_f1,
    stringsAsFactors = FALSE
  )
  
  for (cls in class_names) {
    clean <- gsub("Class: ", "", cls)
    df[[paste0("F1_",          clean)]] <- byclass[cls, "F1"]
    df[[paste0("Sensitivity_", clean)]] <- byclass[cls, "Sensitivity"]
    df[[paste0("Specificity_", clean)]] <- byclass[cls, "Specificity"]
  }
  
  return(df)
}

#===================================================================================
# define FS functions for RFE and r2VIM

# --- RFE ---
var.sel.rfe <- function(x, y, prop.rm = 0.2, recalculate = TRUE, tol = 10,
                        ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                        no.threads = 1, method = "ranger", type = "classification",
                        importance = "impurity_corrected", case.weights = NULL) {
  infos    <- NULL
  info.var <- matrix(ncol = 0, nrow = ncol(x), dimnames = list(colnames(x), NULL))
  var      <- colnames(x)
  imp      <- NULL
  
  while (length(var) >= 2) {
    x.sub <- if (length(var) == 1) matrix(x[, var], ncol = 1) else x[, var]
    rf    <- wrapper.rf(x = x.sub, y = y, ntree = ntree, mtry.prop = mtry.prop,
                        nodesize.prop = nodesize.prop, no.threads = no.threads,
                        method = method, type = type, importance = importance,
                        case.weights = case.weights)
    if (is.null(imp) | recalculate) imp <- get.vim(rf)
    imp   <- sort(imp[var])
    error <- calculate.error(rf = rf, true = y)
    infos <- rbind(infos, c(n = length(var), error))
    temp  <- rep(0, nrow(info.var))
    temp[rownames(info.var) %in% var] <- 1
    info.var <- cbind(info.var, temp)
    no.rm    <- round(prop.rm * length(var))
    var      <- names(imp)[-(1:no.rm)]
  }
  
  infos      <- data.frame(infos)
  error.name <- if (type == "regression") "rmse" else "err"
  best       <- min(infos[, error.name])
  ind.min    <- if (best == 0) {
    max(which(infos[, error.name] == best))
  } else {
    error.prop <- (infos[, error.name] - best) / best * 100
    max(which(error.prop <= tol))
  }
  
  sum     <- apply(info.var, 1, sum)
  ind.sel <- as.numeric(info.var[, ind.min] == 1)
  info    <- data.frame(sum, ind.sel)
  colnames(info) <- c("included.until.subset", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1]), info.runs = infos))
}

# --- r2VIM ---
var.sel.r2vim <- function(x, y, no.runs = 10, factor = 0.5, ntree = 500,
                          mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                          method = "ranger", type = "classification",
                          importance = "impurity_corrected", case.weights = NULL) {
  imp.all <- NULL
  for (r in 1:no.runs) {
    rf      <- wrapper.rf(x = x, y = y, ntree = ntree, mtry.prop = mtry.prop,
                          nodesize.prop = nodesize.prop, no.threads = no.threads,
                          method = method, type = type, importance = importance,
                          case.weights = case.weights)
    imp.all <- cbind(imp.all, get.vim(rf))
  }
  
  min.global <- min(imp.all)
  no.neg.min <- 0
  fac        <- matrix(nrow = nrow(imp.all), ncol = ncol(imp.all),
                       dimnames = dimnames(imp.all))
  for (i in 1:ncol(imp.all)) {
    x.col <- imp.all[, i]
    min.i  <- min(x.col)
    if (min.i >= 0) {
      no.neg.min  <- no.neg.min + 1
      fac[, i]    <- x.col / abs(min.global)
    } else {
      fac[, i]    <- x.col / abs(min.i)
    }
  }
  if (no.neg.min > 0) message(no.neg.min, " runs with no negative importance score")
  
  fac.min <- apply(fac, 1, min)
  fac.med <- apply(fac, 1, median)
  ind.sel <- as.numeric(fac.med >= factor)
  info    <- data.frame(imp.all, fac, fac.min, fac.med, ind.sel)
  colnames(info) <- c(paste0("vim.run.", 1:no.runs),
                      paste0("rel.vim.run.", 1:no.runs),
                      "rel.vim.min", "rel.vim.median", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}

#===================================================================================
# predict impact categories based on the full dataset without prior feature selection using RF and SVM
# Due to the high dimensionality of the functional dataset, retain the 15,000 most variable features based on the training set only.


top_n_all <- 15000

all_metrics_all <- list()

for (i in 1:n_iter) {
  
 set.seed(300 + i)
 splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
 train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:244648]
 test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:244648]
 train$site <- as.factor(train$site)
 test$site  <- as.factor(test$site)
  
  # -------------------------
  # Variance filter on train-split only
 x_train <- train[, colnames(train) != "site"]
  
 feature_vars <- apply(x_train, 2, var)
 top_features <- names(sort(feature_vars, decreasing = TRUE))[1:top_n_all]
  
 train_sub <- train[, c(top_features, "site")]
 test_sub  <- test[,  c(top_features, "site")]
  
 n_feat <- length(top_features)
  
  # -------------------------
  # RF
  set.seed(500 + i)
  rf    <- train(site ~ ., data = train_sub, method = "rf",
                 ntree = 500, tuneLength = 3, metric = "Accuracy",
                 trControl = my_ctrl)
  rf_cm <- confusionMatrix(predict(rf, test_sub), test_sub$site)
  
  # -------------------------
  # SVM
  set.seed(600 + i)
  svm    <- train(site ~ ., data = train_sub, method = "svmLinear",
                  tuneLength = 5, metric = "Accuracy",
                  preProcess = c("center", "scale"), trControl = my_ctrl)
  svm_cm <- confusionMatrix(predict(svm, test_sub), test_sub$site)
  
  all_metrics_all[[i]] <- bind_rows(
    extract_metrics(rf_cm,  i, "RF_ALL",  n_feat),
    extract_metrics(svm_cm, i, "SVM_ALL", n_feat)
 )
  
  cat("Iteration", i, "done (ALL) —", n_feat, "features\n")
}

results_all <- bind_rows(all_metrics_all)
saveRDS(results_all, "results_ALL_SVM_RF.rds")


# ========================================================================================
# run BORUTA (FS) and do prediction of impact categories based on Boruta selected features using RF and SVM

all_metrics_boruta  <- list()
all_features_boruta <- list()

for (i in 1:n_iter) {
  
 set.seed(300 + i)
splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:244648]
test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:244648]
train$site <- as.factor(train$site)
test$site  <- as.factor(test$site)
  

# -------------------------
# Boruta
 set.seed(400 + i)
 bor           <- Boruta(site ~ ., data = train, doTrace = 0, maxRuns = 500)
 selected_vars <- getSelectedAttributes(bor, withTentative = FALSE)
  
 if (length(selected_vars) < 2) { cat("Iteration", i, "skipped (Boruta)\n"); next }
  
 all_features_boruta[[i]] <- data.frame(Iteration = i, Feature = selected_vars)
 train_bor <- train[, c(selected_vars, "site")]
 test_bor  <- test[,  c(selected_vars, "site")]
  
 # -------------------------
 # RF
 set.seed(500 + i)
 rf    <- train(site ~ ., data = train_bor, method = "rf",
               ntree = 500, tuneLength = 3, metric = "Accuracy",
              trControl = my_ctrl)
 rf_cm <- confusionMatrix(predict(rf, test_bor), test_bor$site)
  
 
 # -------------------------
 # SVM
 set.seed(600 + i)
 svm    <- train(site ~ ., data = train_bor, method = "svmLinear",
           tuneLength = 5, metric = "Accuracy",
               preProcess = c("center", "scale"), trControl = my_ctrl)
 svm_cm <- confusionMatrix(predict(svm, test_bor), test_bor$site)

 all_metrics_boruta[[i]] <- bind_rows(
 extract_metrics(rf_cm,  i, "RF_Boruta",  length(selected_vars)),
 extract_metrics(svm_cm, i, "SVM_Boruta", length(selected_vars))
 )
 cat("Iteration", i, "done (Boruta)\n")
}

results_boruta  <- bind_rows(all_metrics_boruta)
features_boruta <- bind_rows(all_features_boruta)
saveRDS(results_boruta,  "results_Boruta_SVM_RF.rds")
saveRDS(features_boruta, "features_Boruta_SVM_RF.rds")







# ========================================================================================
## run RFE (FS) and do prediction of impact categories based on RFE selected features using RF and SVM
all_metrics_rfe  <- list()
all_features_rfe <- list()

for (i in 1:n_iter) {
  
set.seed(300 + i)
splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:244648]
test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:244648]
train$site <- as.factor(train$site)
test$site  <- as.factor(test$site)
  
  x <- train[, colnames(train) != "site"]
  y <- train$site
  
  
  # -------------------------
  # RFE
 set.seed(400 + i)
 res           <- var.sel.rfe(x, y)
selected_vars <- res$var

 if (length(selected_vars) < 2) { cat("Iteration", i, "skipped (RFE)\n"); next }
  
 all_features_rfe[[i]] <- data.frame(Iteration = i, Feature = selected_vars)
 train_rfe <- train[, c(selected_vars, "site")]
 test_rfe  <- test[,  c(selected_vars, "site")]
  
 # -------------------------
 # RF
 set.seed(500 + i)
 rf    <- train(site ~ ., data = train_rfe, method = "rf",
             ntree = 500, tuneLength = 3, metric = "Accuracy",
             trControl = my_ctrl)
rf_cm <- confusionMatrix(predict(rf, test_rfe), test_rfe$site)
  
# -------------------------
# SVM
set.seed(600 + i)
svm    <- train(site ~ ., data = train_rfe, method = "svmLinear",
               tuneLength = 5, metric = "Accuracy",
               preProcess = c("center", "scale"), trControl = my_ctrl)
svm_cm <- confusionMatrix(predict(svm, test_rfe), test_rfe$site)

all_metrics_rfe[[i]] <- bind_rows(
extract_metrics(rf_cm,  i, "RF_RFE",  length(selected_vars)),
extract_metrics(svm_cm, i, "SVM_RFE", length(selected_vars))
)
cat("Iteration", i, "done (RFE)\n")
}

results_rfe  <- bind_rows(all_metrics_rfe)
features_rfe <- bind_rows(all_features_rfe)
saveRDS(results_rfe,  "results_RFE_SVM_RF.rds")
saveRDS(features_rfe, "features_RFE_SVM_RF.rds")




# ========================================================================================
## run PIMP (FS) and do prediction of impact categories based on PIMP selected features using RF and SVM
# due to the high computational power PIMP needs, do the FS based on 10.000 features selected prior to FS based on their variance

top_n_features <- 10000
S_permutations <- 100

all_metrics_pimp  <- list()
all_features_pimp <- list()

for (i in 1:n_iter) {
  
set.seed(300 + i)
splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:244648]
test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:244648]
train$site <- as.factor(train$site)
test$site  <- as.factor(test$site)
  
x <- train[, colnames(train) != "site"]
y <- train$site
  
  # -------------------------
  # Pre-filter by variance
  
feature_vars <- apply(x, 2, var)
top_features <- names(sort(feature_vars, decreasing = TRUE))[1:top_n_features]
x_filtered   <- x[, top_features]
  
  # -------------------------
  # PIMP
 
set.seed(400 + i)
rf_model  <- randomForest(x_filtered, y, ntree = 500, importance = TRUE)
pimp_res  <- PIMP(x_filtered, y, rf_model, S = S_permutations)
pimp_test <- PimpTest(pimp_res)
sum_pimp  <- summary(pimp_test, pless = 0.05)
  
selected_vars <- rownames(sum_pimp$cmat2)[
apply(sum_pimp$cmat2, 1, function(x) any(x < 0.05))
]
selected_vars <- intersect(selected_vars, colnames(train))

if (length(selected_vars) < 2) { 
   cat("Iteration", i, "skipped (PIMP)\n")
   next 
 }
  
 all_features_pimp[[i]] <- data.frame(Iteration = i, Feature = selected_vars)
 train_pimp <- train[, c(selected_vars, "site")]
 test_pimp  <- test[,  c(selected_vars, "site")]
  
  # -------------------------
  # RF
 set.seed(500 + i)
 rf    <- train(site ~ ., data = train_pimp, method = "rf",
                ntree = 500, tuneLength = 3, metric = "Accuracy",
               trControl = my_ctrl)
 rf_cm <- confusionMatrix(predict(rf, test_pimp), test_pimp$site)
  
  # -------------------------
  # SVM
 set.seed(600 + i)
 svm    <- train(site ~ ., data = train_pimp, method = "svmLinear",
                 tuneLength = 5, metric = "Accuracy",
                 preProcess = c("center", "scale"), trControl = my_ctrl)
 svm_cm <- confusionMatrix(predict(svm, test_pimp), test_pimp$site)
  
 all_metrics_pimp[[i]] <- bind_rows(
   extract_metrics(rf_cm,  i, "RF_PIMP",  length(selected_vars)),
   extract_metrics(svm_cm, i, "SVM_PIMP", length(selected_vars))
 )

 cat("Iteration", i, "done (PIMP) —", length(selected_vars), "features selected\n")
}

results_pimp  <- bind_rows(all_metrics_pimp)
features_pimp <- bind_rows(all_features_pimp)
saveRDS(results_pimp,  "results_PIMP_SVM_RF.rds")
saveRDS(features_pimp, "features_PIMP_SVM_RF.rds")




# ========================================================================================
## run r2VIM (FS) and do prediction of impact categories based on r2VIM selected features using RF and SVM


all_metrics_r2vim  <- list()
all_features_r2vim <- list()

for (i in 1:n_iter) {
  
 set.seed(300 + i)
 splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
 train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:244648]
 test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:244648]
 train$site <- as.factor(train$site)
 test$site  <- as.factor(test$site)
  
 x <- train[, colnames(train) != "site"]
 y <- train$site
  
 # -------------------------
 # r2VIM
 set.seed(400 + i) 
 res           <- var.sel.r2vim(x, y, no.runs = 10, factor = 0.5)
 selected_vars <- res$var
  
 #fallback: top fallback_n by median relative importance
 if (length(selected_vars) < 2) {
   cat("Iteration", i, ": using fallback (top", FALLBACK_N, "features)\n")
   imp           <- res$info$rel.vim.median
   names(imp)    <- rownames(res$info)
   selected_vars <- names(sort(imp, decreasing = TRUE))[1:FALLBACK_N]
 }
  
 selected_vars <- intersect(selected_vars, colnames(train))
 if (length(selected_vars) < 2) { cat("Iteration", i, "skipped (r2VIM)\n"); next }
  
 all_features_r2vim[[i]] <- data.frame(Iteration = i, Feature = selected_vars)
 train_r2 <- train[, c(selected_vars, "site")]
 test_r2  <- test[,  c(selected_vars, "site")]
  
 
 # -------------------------
 # RF
 set.seed(500 + i)
 rf    <- train(site ~ ., data = train_r2, method = "rf",
                ntree = 500, tuneLength = 3, metric = "Accuracy",
                trControl = my_ctrl)
 rf_cm <- confusionMatrix(predict(rf, test_r2), test_r2$site)
  
 
 # -------------------------
 # SVM
 set.seed(600 + i)
 svm    <- train(site ~ ., data = train_r2, method = "svmLinear",
                tuneLength = 5, metric = "Accuracy",
                 preProcess = c("center", "scale"), trControl = my_ctrl)
 svm_cm <- confusionMatrix(predict(svm, test_r2), test_r2$site)
  
 all_metrics_r2vim[[i]] <- bind_rows(
   extract_metrics(rf_cm,  i, "RF_r2VIM",  length(selected_vars)),
   extract_metrics(svm_cm, i, "SVM_r2VIM", length(selected_vars))
 )
 cat("Iteration", i, "done (r2VIM)\n")
}

results_r2vim  <- bind_rows(all_metrics_r2vim)
features_r2vim <- bind_rows(all_features_r2vim)
saveRDS(results_r2vim,  "results_r2VIM_SVM_RF.rds")
saveRDS(features_r2vim, "features_r2VIM_SVM_RF.rds")


#####load RDS files containing RF_SVM results
results_all_features <- read_rds("results_ALL_SVM_RF.rds")
results_Boruta <- read_rds("results_Boruta_SVM_RF.rds")
results_RFE <- read_rds("results_RFE_SVM_RF.rds")
results_r2VIM <- read_rds("results_r2VIM_SVM_RF.rds")
results_PIMP <- read_rds("results_PIMP_SVM_RF.rds")


features_Boruta <- read_rds("features_Boruta_SVM_RF.rds")
features_RFE <- read_rds("features_RFE_SVM_RF.rds")
features_r2VIM <- read_rds("features_r2VIM_SVM_RF.rds")
features_PIMP <- read_rds("features_PIMP_SVM_RF.rds")







# calculate feature stability

feature_stability <- function(features_df, n_iter) {
  features_df %>%
    count(Feature) %>%
    mutate(Frequency = n / n_iter) %>%
    arrange(desc(Frequency))
}

stability_boruta <- feature_stability(features_Boruta, n_iter)
stability_rfe    <- feature_stability(features_RFE,    n_iter)
stability_pimp   <- feature_stability(features_PIMP,   n_iter)
stability_r2vim  <- feature_stability(features_r2VIM,  n_iter)


# plot the results
all_results <- bind_rows(
  results_all_features, results_Boruta, results_RFE, results_PIMP, results_r2VIM
)

all_results$Method <- factor(all_results$Method, levels = c(
  "RF_ALL",    "SVM_ALL",
  "RF_Boruta", "SVM_Boruta",
  "RF_RFE",    "SVM_RFE",
  "RF_PIMP",   "SVM_PIMP",
  "RF_r2VIM",  "SVM_r2VIM"
))

method_colors <- c(
  "RF_ALL"     = "#1565C0", "SVM_ALL"     = "#42A5F5",
  "RF_Boruta"  = "#1B5E20", "SVM_Boruta"  = "#66BB6A",
  "RF_RFE"     = "#E65100", "SVM_RFE"     = "#FFA726",
  "RF_PIMP"    = "#4A148C", "SVM_PIMP"    = "#AB47BC",
  "RF_r2VIM"   = "#B71C1C", "SVM_r2VIM"   = "#EF5350"
)

plot_metric <- function(data, metric_col, y_label, title) {
  ggplot(data, aes(x = Method, y = .data[[metric_col]], fill = Method)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, width = 0.5, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1.2) +
    scale_fill_manual(values = method_colors) +
    labs(title = title, x = NULL, y = y_label) +
    theme_bw(base_size = 13) +
    theme(legend.position  = "none",
          axis.text.x      = element_text(angle = 40, hjust = 1),
          panel.grid.minor = element_blank()) +
    ylim(0, 1)
}

p_acc   <- plot_metric(all_results, "Accuracy", "Accuracy", "Accuracy 50 iterations")
p_kappa <- plot_metric(all_results, "Kappa",    "Kappa",    "Kappa 50 iterations")
p_f1    <- plot_metric(all_results, "Mean_F1",  "Macro F1", "Mean F 50 iterations")


p_acc
p_kappa
p_f1

# per-class F1
f1_long <- all_results %>%
  select(Iteration, Method, starts_with("F1_")) %>%
  pivot_longer(cols = starts_with("F1_"), names_to = "Class", values_to = "F1") %>%
  mutate(Class = gsub("F1_", "", Class))

p_f1_class <- ggplot(f1_long, aes(x = Method, y = F1, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, width = 0.5, alpha = 0.85) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.2) +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ Class, nrow = 1) +
  labs(title = "Per-class F1-score 50 iterations", x = NULL, y = "F1 Score") +
  theme_bw(base_size = 13) +
  theme(legend.position  = "none",
        axis.text.x      = element_text(angle = 40, hjust = 1),
        panel.grid.minor = element_blank()) +
  ylim(0, 1)
p_f1_class


ggsave("boxplot_accuracy_BGR_func.pdf",     p_acc,      width = 10, height = 5)
ggsave("boxplot_kappa_BGR_func.pdf",        p_kappa,    width = 10, height = 5)
ggsave("boxplot_meanF1_BGR_func.pdf",       p_f1,       width = 10, height = 5)
ggsave("boxplot_perclass_F1_BGR_func.pdf",  p_f1_class, width = 14, height = 5)



##summarize stats for table

summary_table <- all_results %>%
  group_by(Method) %>%
  summarise(
    # Accuracy
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_sd   = sd(Accuracy, na.rm = TRUE),
    Accuracy_min  = min(Accuracy, na.rm = TRUE),
    Accuracy_max  = max(Accuracy, na.rm = TRUE),
    
    # Kappa
    Kappa_mean = mean(Kappa, na.rm = TRUE),
    Kappa_sd   = sd(Kappa, na.rm = TRUE),
    Kappa_min  = min(Kappa, na.rm = TRUE),
    Kappa_max  = max(Kappa, na.rm = TRUE),
    
    # Mean F1
    F1_mean = mean(Mean_F1, na.rm = TRUE),
    F1_sd   = sd(Mean_F1, na.rm = TRUE),
    F1_min  = min(Mean_F1, na.rm = TRUE),
    F1_max  = max(Mean_F1, na.rm = TRUE)
  ) %>%
  arrange(desc(Accuracy_mean))

summary_table


## mean +-std
summary_table_meanstd <- summary_table %>%
  mutate(
    Accuracy = sprintf("%.3f ± %.3f", Accuracy_mean, Accuracy_sd),
    Kappa    = sprintf("%.3f ± %.3f", Kappa_mean, Kappa_sd),
    F1       = sprintf("%.3f ± %.3f", F1_mean, F1_sd)
  ) %>%
  select(Method, Accuracy, Kappa, F1)

summary_table_meanstd


f1_summary <- f1_long %>%
  group_by(Method, Class) %>%
  summarise(
    mean = mean(F1, na.rm = TRUE),
    sd   = sd(F1, na.rm = TRUE),
    min  = min(F1, na.rm = TRUE),
    max  = max(F1, na.rm = TRUE),
    .groups = "drop"
  )

f1_summary
print



results_all    <- readRDS("results_ALL_SVM_RF.rds")
results_boruta <- readRDS("results_Boruta_SVM_RF.rds")
results_rfe    <- readRDS("results_RFE_SVM_RF.rds")
results_pimp   <- readRDS("results_PIMP_SVM_RF.rds")
results_r2vim  <- readRDS("results_r2VIM_SVM_RF.rds")

all_results <- bind_rows(
  results_all, results_boruta, results_rfe,
  results_pimp, results_r2vim
)

# ============================================================================================
# pairwise wilcoxon rank-sum test to compare all RF FS methods against each other and all SVM FS methods against each other


compare_methods <- function(results_df, metric = "Accuracy",
                            model = "RF") {
  
  df <- results_df %>%
    # RF_ methods for RF, SVM_ methods for SVM
    filter(startsWith(Method, paste0(model, "_"))) %>%
    select(Iteration, Method, all_of(metric))
  
  methods <- unique(df$Method)
  
  # all pairwise combinations
  pairs <- combn(methods, 2, simplify = FALSE)
  
  results <- lapply(pairs, function(pair) {
    
    x <- df %>% filter(Method == pair[1]) %>% pull(!!sym(metric))
    y <- df %>% filter(Method == pair[2]) %>% pull(!!sym(metric))
    
    # match lengths
    min_n <- min(length(x), length(y))
    x <- x[1:min_n]
    y <- y[1:min_n]
    
    test <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
    
    data.frame(
      Method_A   = pair[1],
      Method_B   = pair[2],
      Metric     = metric,
      Model      = model,
      median_A   = round(median(x), 4),
      median_B   = round(median(y), 4),
      W          = test$statistic,
      p_value    = test$p.value
    )
  })
  
  bind_rows(results)
}

# ============================================================================================
# run for all metrics and both models


metrics <- c("Accuracy", "Kappa", "Mean_F1")
models  <- c("RF", "SVM")

stat_results <- bind_rows(lapply(metrics, function(m) {
  bind_rows(lapply(models, function(mod) {
    compare_methods(all_results, metric = m, model = mod)
  }))
}))

# ============================================================================================
# multiple test correction
# use Benjamini-Hochberg 


stat_results <- stat_results %>%
  group_by(Metric, Model) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    significance = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )

print(stat_results)

write.csv(stat_results,
          "pairwise_wilcoxon_RF_SVM_FS_methods.csv",
          row.names = FALSE)


# ============================================================================================
# separate RF vs RF and SVM vs SVM comparisons cleanly
# pairwise wilcoxon signed rank test to compare RF FS methods against SVM FS methods (same split per iteration -> paired)

compare_rf_vs_svm <- function(results_df, metric = "Accuracy") {
  
  fs_methods <- c("ALL", "Boruta", "RFE", "PIMP", "r2VIM")
  
  results <- lapply(fs_methods, function(fs) {
    
    rf_vals  <- results_df %>%
      filter(Method == paste0("RF_",  fs)) %>%
      pull(!!sym(metric))
    
    svm_vals <- results_df %>%
      filter(Method == paste0("SVM_", fs)) %>%
      pull(!!sym(metric))
    
    min_n    <- min(length(rf_vals), length(svm_vals))
    test     <- wilcox.test(rf_vals[1:min_n], svm_vals[1:min_n],
                            paired = TRUE, exact = FALSE)
    
    data.frame(
      FS_Method  = fs,
      Metric     = metric,
      median_RF  = round(median(rf_vals),  4),
      median_SVM = round(median(svm_vals), 4),
      W          = test$statistic,
      p_value    = test$p.value
    )
  })
  
  bind_rows(results) %>%
    mutate(p_adjusted    = p.adjust(p_value, method = "BH"),
           significance  = case_when(
             p_adjusted < 0.001 ~ "***",
             p_adjusted < 0.01  ~ "**",
             p_adjusted < 0.05  ~ "*",
             TRUE               ~ "ns"
           ))
}

# run for all metrics
rf_vs_svm <- bind_rows(lapply(
  c("Accuracy", "Kappa", "Mean_F1"),
  function(m) compare_rf_vs_svm(all_results, metric = m)
))

print(rf_vs_svm)
write.csv(rf_vs_svm, "RF_vs_SVM_wilcoxon_func95sim.csv", row.names = FALSE)


## Plot results



# ============================================================================================
# setup


project_label <- "func_95sim"  # change per project (func_95sim, func_99sim, tax_97sim)
n_iter        <- 50

rf_order  <- c("RF_ALL", "RF_Boruta", "RF_RFE", "RF_PIMP", "RF_r2VIM")
svm_order <- c("SVM_ALL", "SVM_Boruta", "SVM_RFE", "SVM_PIMP", "SVM_r2VIM")
all_order <- c(rf_order, svm_order)

# ============================================================================================
# labels

method_labels_full <- c(
  "RF_ALL"     = "RF x All features",
  "RF_Boruta"  = "RF x Boruta",
  "RF_RFE"     = "RF x RFE",
  "RF_PIMP"    = "RF x PIMP",
  "RF_r2VIM"   = "RF x r2VIM",
  "SVM_ALL"    = "SVM x All features",
  "SVM_Boruta" = "SVM x Boruta",
  "SVM_RFE"    = "SVM x RFE",
  "SVM_PIMP"   = "SVM x PIMP",
  "SVM_r2VIM"  = "SVM x r2VIM"
)

x_labels_rf <- c(
  "RF_ALL"    = "All",
  "RF_Boruta" = "Boruta",
  "RF_RFE"    = "RFE",
  "RF_PIMP"   = "PIMP",
  "RF_r2VIM"  = "r2VIM"
)

x_labels_svm <- c(
  "SVM_ALL"    = "All",
  "SVM_Boruta" = "Boruta",
  "SVM_RFE"    = "RFE",
  "SVM_PIMP"   = "PIMP",
  "SVM_r2VIM"  = "r2VIM"
)

x_labels_all <- c(
  "RF_ALL"     = "RF\nAll",
  "RF_Boruta"  = "RF\nBoruta",
  "RF_RFE"     = "RF\nRFE",
  "RF_PIMP"    = "RF\nPIMP",
  "RF_r2VIM"   = "RF\nr2VIM",
  "SVM_ALL"    = "SVM\nAll",
  "SVM_Boruta" = "SVM\nBoruta",
  "SVM_RFE"    = "SVM\nRFE",
  "SVM_PIMP"   = "SVM\nPIMP",
  "SVM_r2VIM"  = "SVM\nr2VIM"
)

# ============================================================================================
# colors — colorblind-safe palette


rf_colors <- c(
  "RF_ALL"    = "#D3D3D3",
  "RF_Boruta" = "#0072B2",
  "RF_RFE"    = "#56B4E9",
  "RF_PIMP"   = "#009E73",
  "RF_r2VIM"  = "#F0E442"
)

svm_colors <- c(
  "SVM_ALL"    = "#BEBEBE",
  "SVM_Boruta" = "#D55E00",
  "SVM_RFE"    = "#E69F00",
  "SVM_PIMP"   = "#CC79A7",
  "SVM_r2VIM"  = "#000000"
)

all_colors <- c(rf_colors, svm_colors)

# ============================================================================================
# get significant pairs vs Boruta only (best performing method)


get_sig_pairs <- function(stat_results, model, metric,
                          ref_method = "Boruta") {
  stat_results %>%
    filter(Model        == model,
           Metric       == metric,
           significance != "ns") %>%
    filter(grepl(ref_method, Method_A) |
             grepl(ref_method, Method_B)) %>%
    mutate(pair = pmap(list(Method_A, Method_B), c)) %>%
    pull(pair)
}

# ============================================================================================
# boxplot function


make_boxplot <- function(data, methods, colors, labels_x,
                         metric, title, y_label,
                         sig_pairs, n_iter = 50) {
  
  df <- data %>%
    filter(Method %in% methods) %>%
    mutate(Method = factor(Method, levels = methods))
  
  p <- ggplot(df, aes(x = Method, y = .data[[metric]],
                      fill = Method)) +
    geom_boxplot(outlier.shape = 21, outlier.size  = 1.5,
                 outlier.alpha = 0.5, width         = 0.55,
                 alpha         = 0.9, linewidth     = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
    scale_fill_manual(
      values = colors,
      labels = method_labels_full[methods],
      name   = "Model x Feature selection"
    ) +
    scale_x_discrete(labels = labels_x) +
    labs(
      title   = title,
      x       = "Feature selection method",
      y       = y_label,
      caption = paste0("n = ", n_iter, " independent iterations")
    ) +
    theme_bw(base_size = 12) +
    theme( legend.position = "none",
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10),
      axis.title         = element_text(size = 11),
      plot.title         = element_text(size = 12, face = "bold"),
      plot.caption       = element_text(size = 8,  color = "grey50"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    ) +
    ylim(0, 1.15)
  
  if (length(sig_pairs) > 0) {
    p <- p + stat_compare_means(
      comparisons   = sig_pairs,
      method        = "wilcox.test",
      label         = "p.signif",
      tip.length    = 0.01,
      step.increase = 0.08,
      size          = 3.5,
      vjust         = 0.5
    )
  }
  
  return(p)
}

# ============================================================================================
# accuracy boxplot

metrics_info <- list(
  list(metric = "Accuracy", y_label = "Accuracy"),
  list(metric = "Kappa",    y_label = "Cohen's Kappa"),
  list(metric = "Mean_F1",  y_label = "Mean F1 score")
)

plot_list <- list()

for (m in metrics_info) {
  
  sig_rf  <- get_sig_pairs(stat_results, "RF",  m$metric)
  sig_svm <- get_sig_pairs(stat_results, "SVM", m$metric)
  
  p_rf <- make_boxplot(
    data      = all_results,
    methods   = rf_order,
    colors    = rf_colors,
    labels_x  = x_labels_rf,
    metric    = m$metric,
    title     = paste("RF", m$y_label),
    y_label   = m$y_label,
    sig_pairs = sig_rf
  )
  
  p_svm <- make_boxplot(
    data      = all_results,
    methods   = svm_order,
    colors    = svm_colors,
    labels_x  = x_labels_svm,
    metric    = m$metric,
    title     = paste("SVM", m$y_label),
    y_label   = m$y_label,
    sig_pairs = sig_svm
  )
  
  plot_list[[m$metric]] <- list(rf = p_rf, svm = p_svm)
}

combined <- (plot_list$Accuracy$rf | plot_list$Accuracy$svm) /
  (plot_list$Kappa$rf    | plot_list$Kappa$svm)    /
  (plot_list$Mean_F1$rf  | plot_list$Mean_F1$svm)

combined_annotated <- combined +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 13, face = "bold")
  )

combined_annotated
ggsave(paste0("Figure_performance_", project_label, ".pdf"),
       combined_annotated, height = 14, width = 14)

ggsave(paste0("Figure_performance_", project_label, ".png"),
       combined_annotated, height = 14, width = 14, dpi = 300)


# ============================================================================================
# per-class F1 figure


f1_long <- all_results %>%
  select(Iteration, Method, starts_with("F1_")) %>%
  pivot_longer(
    cols      = starts_with("F1_"),
    names_to  = "Class",
    values_to = "F1"
  ) %>%
  mutate(
    Class  = gsub("F1_", "", Class),
    Class  = recode(Class,
                    "Collector.track" = "Collector track",
                    "Plume.impact"    = "Plume impact",
                    "Reference.site"  = "Reference site"),
    Method = factor(Method, levels = all_order)
  )

p_f1_class <- ggplot(f1_long,
                     aes(x = Method, y = F1, fill = Method)) +
  geom_boxplot(
    outlier.shape = 21, outlier.size = 1.2,
    width         = 0.55, alpha      = 0.9,
    linewidth     = 0.4
  ) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.7) +
  scale_fill_manual(
    values = all_colors,
    labels = method_labels_full[all_order],
    name   = "Model x Feature selection"
  ) +
  scale_x_discrete(labels = x_labels_all) +
  facet_wrap(~ Class, nrow = 1) +
  labs(
    title   = "Per-class F1 score across feature selection methods",
    x       = NULL,
    y       = "F1 score",
    caption = paste0("n = ", n_iter, " independent iterations")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9),
    axis.text.x        = element_blank(),
    axis.text.y        = element_text(size = 10),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 11, face = "bold"),
    strip.background   = element_rect(fill = "grey95", color = "grey70")
  ) +
  ylim(0, 1)
p_f1_class
ggsave(paste0("Figure_perclass_F1_", project_label, ".pdf"),
       p_f1_class, height = 6, width = 16)

cat("All figures saved.\n")

saveRDS(
  p_f1_class,
  file = paste0("Figure_perclass_F1_", project_label, ".rds")
)





# ============================================================================================
# setup 


project_label <- "func_95" # change per project (func_95sim, func_99sim, tax_97sim)
n_iter        <- 50

site_colors <- c(
  "Collector.track" = "#D32F2F",
  "Plume.impact"    = "#1976D2",
  "Reference.site"  = "#388E3C"
)

site_labels <- c(
  "Collector.track" = "Collector track",
  "Plume.impact"    = "Plume impact",
  "Reference.site"  = "Reference site"
)

site_shapes <- c(
  "Collector.track" = 16,
  "Plume.impact"    = 17,
  "Reference.site"  = 15
)

# ============================================================================================
# prepare data


clr_full <- data[, !colnames(data) %in% c("site", "station", "sample")]
clr_full <- as.data.frame(lapply(clr_full, as.numeric))
rownames(clr_full) <- rownames(data)

map_sub <- data.frame(
  sample    = rownames(data),
  site      = data$site,
  row.names = rownames(data)
)
map_sub$site <- as.factor(map_sub$site)

# sanity checks
cat("Samples:            ", nrow(clr_full), "\n")
cat("Features:           ", ncol(clr_full), "\n")
cat("Site levels:        ", levels(map_sub$site), "\n")
cat("Any NA in clr_full:", anyNA(clr_full), "\n")
stopifnot(all(rownames(clr_full) == rownames(map_sub)))

# ============================================================================================
# PCA plots
# no legend (for combined figure panels)


make_pca_plot <- function(clr_mat, metadata, label,
                          site_col = "site",
                          colors   = site_colors,
                          shapes   = site_shapes,
                          labels   = site_labels) {
  
  cat("\nRunning PCA:", label, "\n")
  
  pca     <- prcomp(clr_mat, center = FALSE, scale. = FALSE)
  scores  <- as.data.frame(pca$x)
  scores$site <- metadata[[site_col]]
  pct_var <- summary(pca)$importance[2, ] * 100
  
  p <- ggplot(scores, aes(x = PC1, y = PC2,
                          color = site,
                          shape = site,
                          fill  = site)) +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(level = 0.95, geom = "polygon",
                 alpha = 0.15, linetype = "dashed") +
    scale_color_manual(values = colors, labels = labels,
                       name   = "Impact category") +
    scale_fill_manual( values = colors, labels = labels,
                       name   = "Impact category") +
    scale_shape_manual(values = shapes, labels = labels,
                       name   = "Impact category") +
    labs(
      title = label,
      x     = paste0("PC1 (", round(pct_var[1], 1), "%)"),
      y     = paste0("PC2 (", round(pct_var[2], 1), "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      legend.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 11, face = "bold")
    )
  
  return(list(plot = p, pca = pca, scores = scores))
}

# ============================================================================================
# PERMANOVA (global + pairwise (BH adjusted p-value)) + betadisper


run_permanova <- function(clr_mat, metadata, label,
                          site_col = "site",
                          n_perm   = 9999) {
  
  cat("\n=== PERMANOVA:", label, "===\n")
  dist_mat <- dist(clr_mat, method = "euclidean")
  
  set.seed(400)
  adonis_res <- adonis2(dist_mat ~ site, data = metadata,
                        permutations = n_perm)
  cat("\nGlobal PERMANOVA:\n"); print(adonis_res)
  write.csv(as.data.frame(adonis_res),
            paste0("adonis_global_", label, ".csv"))
  
  set.seed(401)
  pairwise_res <- pairwise.adonis2(dist_mat ~ site, data = metadata,
                                   permutations = n_perm)
  cat("\nPairwise PERMANOVA:\n"); print(pairwise_res)
  pairwise_df <- do.call(rbind, lapply(names(pairwise_res), function(nm) {
    if (nm == "parent_call") return(NULL)
    x <- pairwise_res[[nm]]
    if (!is.data.frame(x)) return(NULL)
    df <- as.data.frame(x); df$Comparison <- nm
    df$Metric <- rownames(df); rownames(df) <- NULL; df
  }))
  pairwise_df <- pairwise_df %>%
    mutate(p_adjusted = case_when(
      Metric == "Pr(>F)" ~ p.adjust(`Pr(>F)`, method = "BH"),
      TRUE               ~ NA_real_
    ))
  
  write.csv(pairwise_df, paste0("adonis_pairwise_", label, ".csv"),
            row.names = FALSE)
  
  set.seed(402)
  bd      <- betadisper(dist_mat, metadata[[site_col]])
  bd_test <- permutest(bd, permutations = n_perm)
  cat("\nBetadisper:\n"); print(bd_test)
  
  #save betadisper with p-value
  bd_p   <- bd_test$tab["Groups", "Pr(>F)"]
  bd_F   <- bd_test$tab["Groups", "F"]
  bd_df1 <- bd_test$tab["Groups",    "Df"]
  bd_df2 <- bd_test$tab["Residuals", "Df"]
  
  bd_df <- data.frame(
    Group        = names(bd$group.distances),
    MeanDistance = bd$group.distances
  ) %>%
    mutate(
      F_statistic     = round(bd_F,  3),
      Df_groups       = bd_df1,
      Df_residuals    = bd_df2,
      p_value         = round(bd_p, 4),
      p_value_fmt     = case_when(
        bd_p < 0.001 ~ "< 0.001",
        bd_p < 0.01  ~ "< 0.01",
        bd_p < 0.05  ~ "< 0.05",
        TRUE         ~ as.character(round(bd_p, 3))
      )
    )
  
  write.csv(bd_df, paste0("betadisper_", label, ".csv"),
            row.names = FALSE)
  cat("Results saved:", label, "\n")
  
  return(list(adonis     = adonis_res,
              pairwise   = pairwise_res,
              betadisper = bd_test,
              dist       = dist_mat))
}

# ============================================================================================
# extract consensus features


get_consensus <- function(features_df, n_iter = 50, threshold = 0.5) {
  features_df %>%
    count(Feature) %>%
    mutate(Frequency = n / n_iter) %>%
    filter(Frequency >= threshold) %>%
    arrange(desc(Frequency)) %>%
    pull(Feature)
}

# ============================================================================================
# subset CLR table to consensus features


subset_clr <- function(clr_mat, features, label) {
  available <- intersect(features, colnames(clr_mat))
  cat(label, "— requested:", length(features),
      "| matched:", length(available), "\n")
  if (length(available) < 2) {
    stop(paste(label, "— too few features matched. Check feature ID format."))
  }
  clr_mat[, available, drop = FALSE]
}

# ============================================================================================
# feature stability table


make_stability_table <- function(features_df, method_name, n_iter) {
  features_df %>%
    count(Feature) %>%
    mutate(Method    = method_name,
           Frequency = n / n_iter) %>%
    arrange(desc(Frequency)) %>%
    select(Method, Feature, n, Frequency)
}

# ============================================================================================
# PERMANOVA on full dataset


perm_full <- run_permanova(
  clr_mat  = clr_full,
  metadata = map_sub,
  label    = paste0("All features")
)

# ============================================================================================
# PCA on full dataset


pca_full <- make_pca_plot(
  clr_mat  = clr_full,
  metadata = map_sub,
  label    = "All features"
)
pca_full


# ============================================================================================
# check selection frequency of features to define thresholds for consenus feature sets
plot_selection_freq <- function(features_df, method_name, threshold) {
  features_df %>%
    group_by(Feature) %>%
    summarise(freq = n() / 50) %>%
    ggplot(aes(x = freq)) +
    geom_histogram(bins = 30, fill = "#1F4E79", color = "white", alpha = 0.8) +
    geom_vline(xintercept = threshold, color = "red", 
               linetype = "dashed", linewidth = 1) +
    annotate("text", x = threshold + 0.03, y = Inf, 
             label = paste0("Threshold: ", threshold),
             color = "red", vjust = 2, hjust = 0, size = 3.5) +
    labs(
      title = method_name,
      x = "Selection frequency across 50 iterations",
      y = "Number of features"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

# Generate for each method
p1 <- plot_selection_freq(features_boruta, "Boruta (threshold = 0.30)", 0.30)
p2 <- plot_selection_freq(features_r2vim,  "r2VIM (threshold = 0.50)",  0.50)
p3 <- plot_selection_freq(features_rfe,    "RFE (threshold = 0.10)",    0.10)
p4 <- plot_selection_freq(features_pimp,   "PIMP (threshold = 0.80)",   0.80)


combined_threshold <- (p1 | p2) / (p3 | p4)
combined_threshold
combined_annotated_threshold <- combined_threshold +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 16, face = "bold")
  )

combined_annotated_threshold
ggsave("selection_frequency_distributions.pdf", combined_annotated_threshold, 
       width = 10, height = 7, dpi = 300)





# ============================================================================================
# load and prepare consensus features


features_boruta <- readRDS("features_Boruta_SVM_RF.rds")
features_rfe    <- readRDS("features_RFE_SVM_RF.rds")
features_pimp   <- readRDS("features_PIMP_SVM_RF.rds")
features_r2vim  <- readRDS("features_r2VIM_SVM_RF.rds")

# bind if saved as list
for (obj in c("features_boruta", "features_rfe",
              "features_pimp",   "features_r2vim")) {
  x <- get(obj)
  if (is.list(x) & !is.data.frame(x)) assign(obj, bind_rows(x))
}
# ============================================================================================
# extract consensus features — adjust thresholds per method as needed
consensus_boruta <- get_consensus(features_boruta, n_iter, threshold = 0.3)
consensus_rfe    <- get_consensus(features_rfe,    n_iter, threshold = 0.1)
consensus_pimp   <- get_consensus(features_pimp,   n_iter, threshold = 0.8)
consensus_r2vim  <- get_consensus(features_r2vim,  n_iter, threshold = 0.5)

cat("\n--- Consensus feature counts ---\n")
cat("Boruta:", length(consensus_boruta), "\n")
cat("RFE:   ", length(consensus_rfe),    "\n")
cat("PIMP:  ", length(consensus_pimp),   "\n")
cat("r2VIM: ", length(consensus_r2vim),  "\n")

# verify ID format
cat("\nFirst 3 CLR column names:\n");  print(head(colnames(clr_full),    3))
cat("First 3 Boruta feature IDs:\n"); print(head(consensus_boruta, 3))

# subset CLR tables
clr_boruta <- subset_clr(clr_full, consensus_boruta, "Boruta")
clr_rfe    <- subset_clr(clr_full, consensus_rfe,    "RFE")
clr_pimp   <- subset_clr(clr_full, consensus_pimp,   "PIMP")
clr_r2vim  <- subset_clr(clr_full, consensus_r2vim,  "r2VIM")


# ============================================================================================
# PCA on consenus features subset

pca_boruta <- make_pca_plot(clr_boruta, map_sub, "Boruta consensus features")
pca_rfe    <- make_pca_plot(clr_rfe,    map_sub, "RFE consensus features")
pca_pimp   <- make_pca_plot(clr_pimp,   map_sub, "PIMP consensus features")
pca_r2vim  <- make_pca_plot(clr_r2vim,  map_sub, "r2VIM consensus features")

pca_boruta
pca_rfe
pca_pimp
pca_r2vim
# ============================================================================================
# PERMANOVA on consenus features subset


perm_boruta <- run_permanova(clr_boruta, map_sub, paste0("Boruta consensus features"))
perm_rfe    <- run_permanova(clr_rfe,    map_sub, paste0("RFE consensus features"))
perm_pimp   <- run_permanova(clr_pimp,   map_sub, paste0("PIMP consensus features"))
perm_r2vim  <- run_permanova(clr_r2vim,  map_sub, paste0("r2VIM consensus features"))

# ============================================================================================
# combined PCA figure

legend <- get_legend(
  pca_full$plot + theme(legend.position = "right")
)
legend_plot <- as_ggplot(legend)



comb_pca_plots_legend <- ggarrange(pca_full$plot, pca_boruta$plot, pca_rfe$plot, pca_r2vim$plot, pca_pimp$plot, legend_plot, ncol=3, nrow=2, legend = F, labels = c("a", "b", "c", "d", "e", ""))
comb_pca_plots_legend


ggsave(paste0("PCA_combined_", project_label, ".pdf"),
       comb_pca_plots_legend, height = 10, width = 15)

ggsave(paste0("PCA_combined_", project_label, ".png"),
       comb_pca_plots_legend, height = 10, width = 15, dpi = 300)

cat("Combined PCA saved.\n")

# ============================================================================================
# feature stability table

stability_all <- bind_rows(
  make_stability_table(features_boruta, "Boruta", n_iter),
  make_stability_table(features_rfe,    "RFE",    n_iter),
  make_stability_table(features_pimp,   "PIMP",   n_iter),
  make_stability_table(features_r2vim,  "r2VIM",  n_iter)
)

write.csv(stability_all,
          paste0("feature_stability_all_", project_label, ".csv"),
          row.names = FALSE)
cat("Feature stability table saved.\n")







# ============================================================================================
#  annotations of core consensus features
# load your saved feature lists
features_boruta <- readRDS("features_Boruta_SVM_RF.rds")
features_rfe    <- readRDS("features_RFE_SVM_RF.rds")
features_pimp   <- readRDS("features_PIMP_SVM_RF.rds")
features_r2vim  <- readRDS("features_r2VIM_SVM_RF.rds")

# bind if list
for (obj in c("features_boruta", "features_rfe",
              "features_pimp",   "features_r2vim")) {
  x <- get(obj)
  if (is.list(x) & !is.data.frame(x)) assign(obj, bind_rows(x))
}

n_iter <- 50

# stability table per method
make_stability <- function(features_df, method_name, n_iter) {
  features_df %>%
    count(Feature) %>%
    mutate(
      Method        = method_name,
      Frequency     = n / n_iter,
      Stability_pct = round(Frequency * 100, 1)
    ) %>%
    arrange(desc(Frequency)) %>%
    select(Method, Feature, n, Frequency, Stability_pct) %>%
    rename(
      "FS_method"          = Method,
      "Feature_ID"         = Feature,
      "N_selected"         = n,
      "Selection_frequency"= Frequency,
      "Stability_pct"      = Stability_pct
    )
}

stab_boruta <- make_stability(features_boruta, "Boruta", n_iter)
stab_rfe    <- make_stability(features_rfe,    "RFE",    n_iter)
stab_pimp   <- make_stability(features_pimp,   "PIMP",   n_iter)
stab_r2vim  <- make_stability(features_r2vim,  "r2VIM",  n_iter)

# caps per method
annotation_caps <- list(
  Boruta = Inf,   
  r2VIM  = Inf,   
  RFE    = 50,   
  PIMP   = 100     
)

extract_for_annotation <- function(stab_df, method_name,
                                   threshold, cap) {
  df <- stab_df %>%
    filter(Selection_frequency >= threshold) %>%
    arrange(desc(Selection_frequency))
  
  # apply cap
  if (nrow(df) > cap) {
    cat(method_name, ": capping at top", cap,
        "features (", nrow(df), "in consensus)\n")
    df <- df[1:cap, ]
  } else {
    cat(method_name, ":", nrow(df), "features for annotation\n")
  }
  
  return(df)
}

# extract with caps
annot_boruta <- extract_for_annotation(
  stab_boruta, "Boruta", threshold = 0.3, 
  cap = annotation_caps$Boruta)

annot_r2vim  <- extract_for_annotation(
  stab_r2vim,  "r2VIM",  threshold = 0.5,
  cap = annotation_caps$r2VIM)

annot_rfe    <- extract_for_annotation(
  stab_rfe,    "RFE",    threshold = 0.1,
  cap = annotation_caps$RFE)

annot_pimp   <- extract_for_annotation(
  stab_pimp,   "PIMP",   threshold = 0.8,
  cap = annotation_caps$PIMP)

# union of all IDs for eggNOG input
all_annot_ids <- bind_rows(
  annot_boruta, annot_r2vim,
  annot_rfe,    annot_pimp
) %>%
  distinct(Feature_ID)

cat("\nTotal unique features to annotate:",
    nrow(all_annot_ids), "\n")

# save ID list for eggNOG-mapper
write.table(
  all_annot_ids$Feature_ID,
  paste0("features_for_eggnog_", project_label, ".txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

# save full annotation input table with stability scores
annotation_input <- bind_rows(
  annot_boruta, annot_r2vim,
  annot_rfe,    annot_pimp
)

write.csv(
  annotation_input,
  paste0("annotation_input_", project_label, ".csv"),
  row.names = FALSE
)

cat("Files saved — ready for eggNOG-mapper.\n")




# load eggnog annotations
# skip the ## comment lines, keep #query header
eggnog_lines <- readLines("func_95sim_eggnog.emapper.annotations")

# header line  (#query)
header_line <- which(str_detect(eggnog_lines, "^#query"))[1]
cat("Header at line:", header_line, "\n")

# from header line onwards, skipping ## comment lines above it
eggnog <- read_tsv(
  I(eggnog_lines[header_line:length(eggnog_lines)]),
  col_names = TRUE,
  na = c("", "-", "NA"),
  show_col_types = FALSE
) %>%
  rename(Feature = `#query`)

cat("Rows:", nrow(eggnog), "\n")
cat("Columns:", ncol(eggnog), "\n")
print(colnames(eggnog))
head(eggnog$Feature)


# Quick check
glimpse(eggnog)
cat("Total annotated features:", nrow(eggnog), "\n")
cat("Features with COG category:", sum(!is.na(eggnog$COG_category)), "\n")
cat("Features with KEGG_ko:", sum(!is.na(eggnog$KEGG_ko)), "\n")
cat("Features with KEGG_Pathway:", sum(!is.na(eggnog$KEGG_Pathway)), "\n")

# load consenus feature sets

# Func 95% — Boruta consensus (≥30%)
consensus_boruta_95 <- stability_boruta %>%
  filter(Frequency >= 0.30) %>%
  pull(Feature)
cat("Boruta 95% consensus features:", length(consensus_boruta_95), "\n") # expect 13

# Func 95% — r2VIM consensus (≥50%)
consensus_r2vim_95 <- stability_r2vim %>%
  filter(Frequency >= 0.50) %>%
  pull(Feature)
cat("r2VIM 95% consensus features:", length(consensus_r2vim_95), "\n") # expect 7


# join consenus features and annotations

annotate_consensus <- function(consensus_features, eggnog_df, label) {
  df <- data.frame(Feature = consensus_features) %>%
    left_join(eggnog_df, by = "Feature") %>%
    mutate(Dataset = label)
  cat("\n===", label, "===\n")
  cat("Total features:", nrow(df), "\n")
  cat("Annotated:", sum(!is.na(df$COG_category)), "\n")
  cat("Unannotated:", sum(is.na(df$COG_category)), "\n")
  df
}

ann_boruta_95 <- annotate_consensus(consensus_boruta_95, eggnog, "Func_95_Boruta")
ann_r2vim_95  <- annotate_consensus(consensus_r2vim_95,  eggnog, "Func_95_r2VIM")


# shared Boruta r2VIM features
overlap_95 <- intersect(consensus_boruta_95, consensus_r2vim_95)
cat("\nOverlap Boruta ∩ r2VIM (func 95%):", length(overlap_95), "features\n")
cat("Features:", paste(overlap_95, collapse = ", "), "\n")

overlap_95_annotated <- data.frame(Feature = overlap_95) %>%
  left_join(eggnog, by = "Feature")
print(overlap_95_annotated %>% select(Feature, COG_category, Description, Preferred_name, KEGG_ko, KEGG_Pathway))


# COG category summary
cog_descriptions <- c(
  "A" = "RNA processing and modification",
  "B" = "Chromatin structure and dynamics",
  "C" = "Energy production and conversion",
  "D" = "Cell cycle control and mitosis",
  "E" = "Amino acid metabolism and transport",
  "F" = "Nucleotide metabolism and transport",
  "G" = "Carbohydrate metabolism and transport",
  "H" = "Coenzyme metabolism",
  "I" = "Lipid metabolism",
  "J" = "Translation",
  "K" = "Transcription",
  "L" = "Replication and repair",
  "M" = "Cell wall/membrane/envelope biogenesis",
  "N" = "Cell motility",
  "O" = "Post-translational modification, protein turnover, chaperones",
  "P" = "Inorganic ion transport and metabolism",
  "Q" = "Secondary metabolite biosynthesis, transport, catabolism",
  "R" = "General function prediction only",
  "S" = "Function unknown",
  "T" = "Signal transduction",
  "U" = "Intracellular trafficking and secretion",
  "V" = "Defense mechanisms",
  "W" = "Extracellular structures",
  "Z" = "Cytoskeleton"
)

summarise_cog <- function(ann_df, label) {
  ann_df %>%
    filter(!is.na(COG_category)) %>%
    # Split multi-category assignments (e.g. "PQ")
    mutate(COG_category = strsplit(COG_category, "")) %>%
    unnest(COG_category) %>%
    count(COG_category, sort = TRUE) %>%
    mutate(
      Description = cog_descriptions[COG_category],
      Dataset = label
    )
}

cog_boruta_95 <- summarise_cog(ann_boruta_95, "Func_95_Boruta")
cog_r2vim_95  <- summarise_cog(ann_r2vim_95,  "Func_95_r2VIM")
cog_all_95    <- bind_rows(cog_boruta_95, cog_r2vim_95)

cat("\n=== COG categories — Func 95% Boruta consensus ===\n")
print(cog_boruta_95)

cat("\n=== COG categories — Func 95% r2VIM consensus ===\n")
print(cog_r2vim_95)

# KEGG pathway summary
summarise_kegg <- function(ann_df, label) {
  ann_df %>%
    filter(!is.na(KEGG_Pathway)) %>%
    separate_rows(KEGG_Pathway, sep = ",") %>%
    filter(str_detect(KEGG_Pathway, "^map")) %>% # keep map IDs only
    count(KEGG_Pathway, sort = TRUE) %>%
    mutate(Dataset = label)
}

kegg_boruta_95 <- summarise_kegg(ann_boruta_95, "Func_95_Boruta")
kegg_r2vim_95  <- summarise_kegg(ann_r2vim_95,  "Func_95_r2VIM")

cat("\n=== KEGG Pathways — Func 95% Boruta consensus ===\n")
print(kegg_boruta_95)


# Combine all consensus sets into one clean table
full_annotation_table <- bind_rows(
  ann_boruta_95, ann_r2vim_95
) %>%
  select(Dataset, Feature, COG_category, Description,
         Preferred_name, KEGG_ko, KEGG_Pathway, PFAMs) %>%
  arrange(Dataset, COG_category)

write_csv(full_annotation_table, "consensus_features_annotated.csv")
cat("\nFull annotation table saved to consensus_features_annotated.csv\n")





# Boruta + r2VIM shared core (8 features)
core_features <- c(
  "k141_1576221_4",  # queG
  "k141_438091_1",   # phage tail
  "k141_424361_3",   # carA
  "k141_95205_1",    # mutB
  "k141_718818_2",   # Co/Zn/Cd efflux
  "k141_205392_4",   # KaiC
  "k141_1318015_4",  # pilus
  "k141_402012_4_1"  # unannotated
)

# calculate mean clr abundance per imapct class
clr_long <- clr_trans_table_func %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(site %>% select(sample, site), by = "sample") %>%
  pivot_longer(cols = -c(sample, site),
               names_to = "Feature", values_to = "CLR") %>%
  filter(Feature %in% core_features)

# Mean CLR per feature per impact class
mean_clr <- clr_long %>%
  group_by(Feature, site) %>%
  summarise(
    mean_CLR = mean(CLR),
    sd_CLR   = sd(CLR),
    n        = n(),
    .groups  = "drop"
  ) %>%
  left_join(
    data.frame(
      Feature = core_features,
      Gene    = c("queG","phage tail","carA","mutB",
                  "Co/Zn/Cd efflux","KaiC","pilus","unannotated")
    ), by = "Feature"
  )

print(mean_clr %>% arrange(Feature, site))

# kruskal-wallis and pairwise wilcoxon per feature
kw_results <- clr_long %>%
  group_by(Feature) %>%
  summarise(
    kw_p = kruskal.test(CLR ~ site)$p.value,
    .groups = "drop"
  )

# Pairwise Wilcoxon per feature
pw_results <- clr_long %>%
  group_by(Feature) %>%
  do({
    pw <- pairwise.wilcox.test(.$CLR, .$site,
                               p.adjust.method = "BH")
    as.data.frame(pw$p.value) %>%
      rownames_to_column("Group1") %>%
      pivot_longer(-Group1, names_to = "Group2", values_to = "p_adj") %>%
      filter(!is.na(p_adj))
  }) %>%
  ungroup()

print(kw_results)
print(pw_results)

# apply the same to the Boruta features
all_boruta <- stability_boruta %>%
  filter(Frequency >= 0.30) %>%
  pull(Feature)

mean_clr_all_boruta <- clr_trans_table_func %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(site %>% select(sample, site), by = "sample") %>%
  pivot_longer(cols = -c(sample, site),
               names_to = "Feature", values_to = "CLR") %>%
  filter(Feature %in% all_boruta) %>%
  group_by(Feature, site) %>%
  summarise(mean_CLR = mean(CLR), .groups = "drop") %>%
  pivot_wider(names_from = site, values_from = mean_CLR)

print(mean_clr_all_boruta)

# save
write_csv(mean_clr, "core_features_mean_CLR_by_impact.csv")
write_csv(mean_clr_all_boruta, "boruta_consensus_mean_CLR_by_impact.csv")
