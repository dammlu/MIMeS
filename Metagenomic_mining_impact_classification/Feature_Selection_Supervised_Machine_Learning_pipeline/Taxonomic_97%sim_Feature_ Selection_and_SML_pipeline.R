# M.Sc. Lukas Damm
# RPTU Kaiserslautern-Landau
# 07.03.2025

# curate metagenomic count tables of assembled functional genes and reconstructed/assembled taxonomic genes (16S) 
# eDNA extracted from sediment samples (MUC&PUC) colelcted in the Clarion-Clipperton-Fracture Zone in the Pacific 
# samples are from the BGR area (SO295 cruise), 1.5 years post polymmetalic nodule collector test
# count table is merged kallisto output (tpm)

# use CLR-transformed count tables as input for Feature Selection (FS) and Supervised Machine Learning (SML) pipeline. 
# evaluate four different FS methods (Boruta, RFE, PIMP, r2VIM) and use selected features to classify mining impact categories (Collector impact, Plume impact and undisturbed Reference)
# using Random Forest (RF) and Support Vector Machine (SVM) based on microbial signatures (taxonomic marker genes (16S))



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


# Load data and perform CLR trans

count_table_tax <- read_tsv("merged_tpm_matrix_BGR_97sim.tsv") %>%
  column_to_rownames("target_id")

count_table_tax_t <- count_table_tax %>% t() %>% as.data.frame()

map <- read.csv("mapping_DNA_BGR.csv")

count_table_tax_fin <- count_table_tax_t %>%
  rownames_to_column("sample") %>%
  inner_join(map) %>%
  column_to_rownames("sample")

count_table_tax_fin <- count_table_tax_fin[, c(1:929)] %>%
  t() %>% as.data.frame()



count_table_tax_fin_t <- t(count_table_tax_fin) %>% as.data.frame()

#===================================================================================
# calc overall sparsity
sparsity <- sum(count_table_tax_fin_t == 0) / prod(dim(count_table_tax_fin_t))
cat("Overall sparsity:", round(sparsity * 100, 2), "%\n")

# per-feature sparsity 
feature_sparsity <- colMeans(count_table_tax_fin_t == 0)
summary(feature_sparsity)
hist(feature_sparsity, breaks = 50, 
     main = "Per-feature sparsity", 
     xlab = "Proportion of zero samples")

# per-sample sparsity
sample_sparsity <- rowMeans(count_table_tax_fin_t == 0)
summary(sample_sparsity)

#remove rows where rowSums is zero
count_table_tax_zero  <- count_table_tax_fin[rowSums(count_table_tax_fin) != 0, ]

#add small pseudocount for following CLR trans
count_table_tax_count <- count_table_tax_zero + 0.00001


#===================================================================================
#perform CLR transformation
count_table_tax_clr   <- clr(count_table_tax_count)

write.table(count_table_tax_clr, file = "clr_transformed_count_table_tax_97.tsv",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

clr_trans_table_tax <- read_tsv("clr_transformed_count_table_tax_97.tsv") %>%
  rename(target_id = 1) %>%
  column_to_rownames("target_id") %>%
  t() %>% as.data.frame()

site <- map[, c(1, 3, 7)]
site$station <- as.factor(as.character(site$station))

#make sure $site is a factor
data <- clr_trans_table_tax %>%
  rownames_to_column("sample") %>%
  left_join(site) %>%
  column_to_rownames("sample")

data[, -ncol(data)] <- lapply(data[, -ncol(data)], as.numeric)
data$site <- make.names((data$site))
data$site    <- as.factor(data$site)
data$station <- as.factor(as.character(data$station))
data$sample  <- rownames(data)
data[is.na(data)] <- 0.00001

write.csv(data, "data_table_BGR_tax_97_site.csv")

#===================================================================================
#define functions that are needed throughout the FS+SML 

n_iter <- 50
FALLBACK_N <- 20   

my_ctrl <- trainControl(
  method          = "cv",
  number          = 10,
  classProbs      = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary
)


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


all_metrics_all <- list()

for (i in 1:n_iter) {
  
  set.seed(300 + i)
  splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
  train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:926]
  test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:926]
  train$site <- as.factor(train$site)
  test$site  <- as.factor(test$site)
  n_feat <- ncol(train) - 1
  
  # -------------------------
  # RF
  set.seed(500 + i)
  rf    <- train(site ~ ., data = train, method = "rf",
                 ntree = 500, tuneLength = 3, metric = "Accuracy",
                 trControl = my_ctrl)
  rf_cm <- confusionMatrix(predict(rf, test), test$site)
  
  # -------------------------
  # SVM
  set.seed(600 + i)
  svm    <- train(site ~ ., data = train, method = "svmLinear",
                  tuneLength = 5, metric = "Accuracy",
                  preProcess = c("center", "scale"), trControl = my_ctrl)
  svm_cm <- confusionMatrix(predict(svm, test), test$site)
  
  all_metrics_all[[i]] <- bind_rows(
    extract_metrics(rf_cm,  i, "RF_ALL",  n_feat),
    extract_metrics(svm_cm, i, "SVM_ALL", n_feat)
  )
  cat("Iteration", i, "done (ALL)\n")
}

results_all <- bind_rows(all_metrics_all)
saveRDS(results_all, "results_ALL.rds")


# ========================================================================================
# run BORUTA (FS) and do prediction of impact categories based on Boruta selected features using RF and SVM


all_metrics_boruta  <- list()
all_features_boruta <- list()

for (i in 1:n_iter) {
  
  set.seed(300 + i)
  splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
  train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:926]
  test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:926]
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
saveRDS(results_boruta,  "results_Boruta.rds")
saveRDS(features_boruta, "features_Boruta.rds")


# ========================================================================================
## run RFE (FS) and do prediction of impact categories based on RFE selected features using RF and SVM

all_metrics_rfe  <- list()
all_features_rfe <- list()

for (i in 1:n_iter) {
  
  set.seed(300 + i)
  splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
  train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:926]
  test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:926]
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
saveRDS(results_rfe,  "results_RFE.rds")
saveRDS(features_rfe, "features_RFE.rds")


# ========================================================================================
## run PIMP (FS) and do prediction of impact categories based on PIMP selected features using RF and SVM


all_metrics_pimp  <- list()
all_features_pimp <- list()

for (i in 1:n_iter) {
  
  set.seed(300 + i)
  splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
  train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:926]
  test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:926]
  train$site <- as.factor(train$site)
  test$site  <- as.factor(test$site)
  
  x <- train[, colnames(train) != "site"]
  y <- train$site
  
  # -------------------------
  # PIMP
  set.seed(400 + i)
  rf_model  <- randomForest(x, y, ntree = 500, importance = TRUE)
  pimp_res  <- PIMP(x, y, rf_model, S = 50)
  pimp_test <- PimpTest(pimp_res)
  sum_pimp  <- summary(pimp_test, pless = 0.1)
  
  # robust extraction across all classes
  selected_vars <- rownames(sum_pimp$cmat2)[
    apply(sum_pimp$cmat2, 1, function(x) any(x < 0.1))
  ]
  selected_vars <- intersect(selected_vars, colnames(train))
  
  if (length(selected_vars) < 2) { cat("Iteration", i, "skipped (PIMP)\n"); next }
  
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
  cat("Iteration", i, "done (PIMP)\n")
}

results_pimp  <- bind_rows(all_metrics_pimp)
features_pimp <- bind_rows(all_features_pimp)
saveRDS(results_pimp,  "results_PIMP.rds")
saveRDS(features_pimp, "features_PIMP.rds")


# ========================================================================================
## run r2VIM (FS) and do prediction of impact categories based on r2VIM selected features using RF and SVM


all_metrics_r2vim  <- list()
all_features_r2vim <- list()

for (i in 1:n_iter) {
  
  set.seed(300 + i)
  splits <- partition(data, p = 0.6, id_col = "station", cat_col = "site")
  train  <- splits[[1]] %>% column_to_rownames("sample") %>% .[, 2:926]
  test   <- splits[[2]] %>% column_to_rownames("sample") %>% .[, 2:926]
  train$site <- as.factor(train$site)
  test$site  <- as.factor(test$site)
  
  x <- train[, colnames(train) != "site"]
  y <- train$site
  
  # -------------------------
  # r2VIM
  set.seed(400 + i) 
  res           <- var.sel.r2vim(x, y, no.runs = 10, factor = 0.5)
  selected_vars <- res$var
  
  # fallback: top fallback_n by median relative importance
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
saveRDS(results_r2vim,  "results_r2VIM.rds")
saveRDS(features_r2vim, "features_r2VIM.rds")


#####load RDS files containing RF_SVM results
results_all    <- readRDS("results_ALL.rds")
results_boruta <- readRDS("results_Boruta.rds")
results_rfe    <- readRDS("results_RFE.rds")
results_pimp   <- readRDS("results_PIMP.rds")
results_r2vim  <- readRDS("results_r2VIM.rds")


features_Boruta <- read_rds("features_Boruta.rds")
features_RFE <- read_rds("features_RFEF.rds")
features_r2VIM <- read_rds("features_r2VIM.rds")
features_PIMP <- read_rds("features_PIMP.rds")



# calculate feature stability

feature_stability <- function(features_df, n_iter) {
  features_df %>%
    count(Feature) %>%
    mutate(Frequency = n / n_iter) %>%
    arrange(desc(Frequency))
}

stability_boruta <- feature_stability(features_boruta, n_iter)
stability_rfe    <- feature_stability(features_rfe,    n_iter)
stability_pimp   <- feature_stability(features_pimp,   n_iter)
stability_r2vim  <- feature_stability(features_r2vim,  n_iter)


# plot the results

all_results <- bind_rows(
  results_all, results_boruta, results_rfe, results_pimp, results_r2vim
)

# per-class F1
f1_long <- all_results %>%
  select(Iteration, Method, starts_with("F1_")) %>%
  pivot_longer(cols = starts_with("F1_"), names_to = "Class", values_to = "F1") %>%
  mutate(Class = gsub("F1_", "", Class))


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



#statistics 

results_all    <- readRDS("results_ALL.rds")
results_boruta <- readRDS("results_Boruta.rds")
results_rfe    <- readRDS("results_RFE.rds")
results_pimp   <- readRDS("results_PIMP.rds")
results_r2vim  <- readRDS("results_r2VIM.rds")

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
    
    # match lengths in case some iterations were skipped
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
write.csv(rf_vs_svm, "RF_vs_SVM_wilcoxon_tax97.csv", row.names = FALSE)


## Plot results



# ============================================================================================
# setup

project_label <- "tax_97"  # change per project
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
# colors colorblind-safe palette

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
# combined accuracy figure

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


project_label <- "tax_97"
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

# verify all numeric
clr_full <- as.data.frame(lapply(clr_full, as.numeric))
rownames(clr_full) <- rownames(data)

# metadata aligned to clr_full row order
map_sub <- data.frame(
  sample = rownames(data),
  site   = data$site,
  row.names = rownames(data)
)
map_sub$site <- as.factor(map_sub$site)

# sanity checks
cat("Samples:     ", nrow(clr_full), "\n")
cat("Features:    ", ncol(clr_full), "\n")
cat("Site levels: ", levels(map_sub$site), "\n")
cat("Any NA in clr_full:", anyNA(clr_full), "\n")
stopifnot(all(rownames(clr_full) == rownames(map_sub)))

# ============================================================================================
# PCA plot


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
  cat(label, "— features requested:", length(features),
      "| matched:", length(available), "\n")
  if (length(available) < 2) {
    stop(paste(label, "— too few features matched. Check feature ID format."))
  }
  clr_mat[, available, drop = FALSE]
}

# ============================================================================================
# PERMANOVA on full dataset

perm_full <- run_permanova(
  clr_mat  = clr_full,
  metadata = map_sub,
  label    = paste0(project_label, "_full")
)

# ============================================================================================
# PCA on full dataset

pca_full <- make_pca_plot(
  clr_mat  = clr_full,
  metadata = map_sub,
  label    = paste0("All features")
)
pca_full

# ============================================================================================
# load consensus features

features_boruta <- readRDS("features_Boruta.rds")
features_rfe    <- readRDS("features_RFE.rds")
features_pimp   <- readRDS("features_PIMP.rds")
features_r2vim  <- readRDS("features_r2VIM.rds")

# bind in case saved as list of dataframes
if (is.list(features_boruta) & !is.data.frame(features_boruta)) {
  features_boruta <- bind_rows(features_boruta)
}
if (is.list(features_rfe) & !is.data.frame(features_rfe)) {
  features_rfe <- bind_rows(features_rfe)
}
if (is.list(features_pimp) & !is.data.frame(features_pimp)) {
  features_pimp <- bind_rows(features_pimp)
}
if (is.list(features_r2vim) & !is.data.frame(features_r2vim)) {
  features_r2vim <- bind_rows(features_r2vim)
}

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
# extract consensus at thresholds
consensus_boruta <- get_consensus(features_boruta, n_iter, threshold = 0.3)
consensus_rfe    <- get_consensus(features_rfe,    n_iter, threshold = 0.1)
consensus_pimp   <- get_consensus(features_pimp,   n_iter, threshold = 0.8)
consensus_r2vim  <- get_consensus(features_r2vim,  n_iter, threshold = 0.5)

cat("\n--- Consensus feature counts (>=50% selection frequency) ---\n")
cat("Boruta:", length(consensus_boruta), "\n")
cat("RFE:   ", length(consensus_rfe),    "\n")
cat("PIMP:  ", length(consensus_pimp),   "\n")
cat("r2VIM: ", length(consensus_r2vim),  "\n")

# verify ID format before subsetting
cat("\nFirst 3 CLR column names:\n"); print(head(colnames(clr_full), 3))
cat("First 3 Boruta feature IDs:\n"); print(head(consensus_boruta,  3))

# subset CLR tables
clr_boruta <- subset_clr(clr_full, consensus_boruta, "Boruta")
clr_rfe    <- subset_clr(clr_full, consensus_rfe,    "RFE")
clr_pimp   <- subset_clr(clr_full, consensus_pimp,   "PIMP")
clr_r2vim  <- subset_clr(clr_full, consensus_r2vim,  "r2VIM")

# ============================================================================================
# PCA on consensus subset

pca_boruta <- make_pca_plot(clr_boruta, map_sub, paste0("Boruta consensus features"))
pca_rfe    <- make_pca_plot(clr_rfe,    map_sub, paste0("RFE consensus features"))
pca_pimp   <- make_pca_plot(clr_pimp,   map_sub, paste0("PIMP consensus features"))
pca_r2vim  <- make_pca_plot(clr_r2vim,  map_sub, paste0("r2VIM consensus features"))

pca_boruta
pca_r2vim
pca_full
pca_rfe
pca_pimp

# ============================================================================================
# PERMANOVA on consensus subset

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

make_stability_table <- function(features_df, method_name, n_iter) {
  features_df %>%
    count(Feature) %>%
    mutate(Method    = method_name,
           Frequency = n / n_iter) %>%
    arrange(desc(Frequency)) %>%
    select(Method, Feature, n, Frequency)
}

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
# load SILVA BLAST output
# Standard BLAST tabular format columns:
# qseqid sseqid stitle pident length evalue bitscore
silva_blast <- read_tsv(
  "blast_tax_vs_silva.txt",  
  col_names = c("qseqid", "sseqid", "stitle", 
                "pident", "length", "evalue", "bitscore"),
  show_col_types = FALSE
)

# Keep only best hit per query (highest bitscore)
silva_best <- silva_blast %>%
  group_by(qseqid) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

# parse taxonomy from stitle
# stitle format: "accession Bacteria;Phylum;Class;Order;Family;Genus;description"
silva_tax <- silva_best %>%
  mutate(
    # Extract taxonomy string — everything after first space up to last semicolon block
    tax_string = str_extract(stitle, "(?<=\\d\\s)[Bacteria|Archaea].*"),
    # Split into levels
    tax_levels = str_split(tax_string, ";")
  ) %>%
  mutate(
    Domain  = map_chr(tax_levels, ~ .x[1] %>% str_trim(), .default = NA),
    Phylum  = map_chr(tax_levels, ~ .x[2] %>% str_trim(), .default = NA),
    Class   = map_chr(tax_levels, ~ .x[3] %>% str_trim(), .default = NA),
    Order   = map_chr(tax_levels, ~ .x[4] %>% str_trim(), .default = NA),
    Family  = map_chr(tax_levels, ~ .x[5] %>% str_trim(), .default = NA),
    Genus   = map_chr(tax_levels, ~ .x[6] %>% str_trim(), .default = NA)
  ) %>%
  select(Feature = qseqid, pident, bitscore, 
         Domain, Phylum, Class, Order, Family, Genus)

cat("Total OTUs with taxonomy:", nrow(silva_tax), "\n")

# ============================================================================================
# join tax to consenus features

consensus_boruta_tax <- stability_boruta %>%
  filter(Frequency >= 0.30) %>%
  pull(Feature)

consensus_r2vim_tax <- stability_r2vim %>%
  filter(Frequency >= 0.50) %>%
  pull(Feature)

cat("Tax Boruta consensus:", length(consensus_boruta_tax), "features\n")
cat("Tax r2VIM consensus:", length(consensus_r2vim_tax), "features\n")   

# ============================================================================================
# annotate
tax_boruta_ann <- data.frame(Feature = consensus_boruta_tax) %>%
  left_join(silva_tax, by = "Feature") %>%
  mutate(Method = "Boruta")

tax_r2vim_ann <- data.frame(Feature = consensus_r2vim_tax) %>%
  left_join(silva_tax, by = "Feature") %>%
  mutate(Method = "r2VIM")

# ============================================================================================
# shjared tax between Boruta and r2VM consenus features
overlap_tax <- intersect(consensus_boruta_tax, consensus_r2vim_tax)
cat("\nOverlap Boruta ∩ r2VIM (tax 97%):", length(overlap_tax), "features\n")

# ============================================================================================
# Phylum summary
phylum_summary <- bind_rows(tax_boruta_ann, tax_r2vim_ann) %>%
  filter(!is.na(Phylum)) %>%
  count(Method, Phylum, sort = TRUE)

cat("\n=== Phylum distribution — tax consensus features ===\n")
print(phylum_summary)

# ============================================================================================
# complete table
options(na.print = "—")

tax_r2vim_ann %>%
  select(Feature, Phylum, Class, Order, pident) %>%
  arrange(Phylum) %>%
  as.data.frame() %>%
  print()

cat("\n=== Boruta tax consensus — taxonomy ===\n")
tax_boruta_ann %>%
  select(Feature, Phylum, Class, Order, pident) %>%
  arrange(Phylum) %>%
  as.data.frame() %>%
  print()

cat("\n=== r2VIM tax consensus — taxonomy ===\n")
tax_r2vim_ann %>%
  select(Feature, Phylum, Class, Order, pident, Method) %>%
  mutate(across(everything(), ~replace_na(as.character(.), "—"))) %>%
  arrange(Phylum) %>%
  print(n = 20)

#overlap taxonomy
cat("\n=== Shared features (Boruta ∩ r2VIM) — taxonomy ===\n")
data.frame(Feature = overlap_tax) %>%
  left_join(silva_tax, by = "Feature") %>%
  select(Feature, Phylum, Class, Order, pident) %>%
  mutate(across(everything(), ~replace_na(as.character(.), "—"))) %>%
  arrange(Phylum) %>%
  print(n = 20)
