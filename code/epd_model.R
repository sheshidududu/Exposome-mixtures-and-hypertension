# =============================================================================
# E-P-D model (baseline covariates + clinical indicators + 70 pollutants
# + Gene (protein) biological feature matrix)
# Stacked generalisation pipeline with elastic-net gene feature screening
# =============================================================================

rm(list = ls())

# Load packages
pacman::p_load(
  tidyverse, readxl, naniar, R6, rlist, mlr3verse, mlr3tuning, mlr3pipelines,
  DALEX, DALEXtra, ddpcr, caret, xgboost, ggplot2, patchwork, ggsci, vip, future,
  furrr, writexl, mice, pROC, glmnet
)

seed <- 20260215
set.seed(seed)

# -----------------------------------------------------------------------------
# Input (load before running):
#   exposure data with 70 creatinine-adjusted pollutants, baseline covariates,
#   seven fasting clinical indicators and the binary outcome "Hypertension"
#   gene_matrix_raw: pollutant-gene matrix from the ExposomeX platform
#   (rows: exposure IDs "EX:E...", columns: GO terms, entries 0/1)
# -----------------------------------------------------------------------------

eSet <- R6::R6Class("eSet", lock_class = FALSE, lock_objects = FALSE)
# eSet$Data <- ...   # load data here
# gene_matrix_raw <- ...   # load Gene matrix here

VarY <- "Hypertension"

# Pollutant abbreviation -> ExposomeX exposure ID
pollutant_mapping <- c(
  "SFM" = "EX:E000005038", "4AN6C2B" = "EX:E000060919", "SFSX" = "EX:E000005055",
  "TMTH" = "EX:E000005278", "SFTHI" = "EX:E000005052", "NOR" = "EX:E000004290",
  "ENOR" = "EX:E000065012", "CIP" = "EX:E000002599", "OFLX" = "EX:E000004331",
  "TET" = "EX:E028250140", "OTC" = "EX:E028250143", "AM" = "EX:E028250141",
  "DOXY" = "EX:E028245645", "IMI" = "EX:E051947117", "NIT" = "EX:E002558306",
  "ACE" = "EX:E000191299", "THM" = "EX:E000098621", "CLO" = "EX:E051947118",
  "THI" = "EX:E000105427", "IPP" = "EX:E023920745", "CYC" = "EX:E024269009",
  "DIN" = "EX:E090574179", "IMID" = "EX:E000166282", "DM-THM" = "EX:E008359588",
  "DN-IMI" = "EX:E008578608", "5-OH-IMI" = "EX:E083973165", "DM-ACE" = "EX:E009759597",
  "FP" = "EX:E000003163", "DMM" = "EX:E000078878", "CBDZ" = "EX:E000023751",
  "CPS-O" = "EX:E000020455", "DZN" = "EX:E000002567", "DDVP" = "EX:E000002860",
  "MALT" = "EX:E000003780", "3OH-CBF" = "EX:E000026053", "ATZ" = "EX:E000002109",
  "ATZ-DE-DSI-2HD" = "EX:E090564678", "BPF" = "EX:E000011447", "BPS" = "EX:E000006270",
  "BPC" = "EX:E000006264", "BPE" = "EX:E000538902", "BPP" = "EX:E000558126",
  "2-OHN" = "EX:E000008218", "1-OHN" = "EX:E000006630", "3-OHF" = "EX:E000088028",
  "2-OHF" = "EX:E000068871", "2-OHP" = "EX:E000062971", "1-OHP" = "EX:E000090236",
  "4-OHP" = "EX:E000075021", "1-OHPYR" = "EX:E000020066", "3-PBA" = "EX:E000018334",
  "4-F-3-PBA" = "EX:E000142364", "2-MB-3CA" = "EX:E013156095",
  "trans-DCCA" = "EX:E000015745", "cis-DCCA" = "EX:E000019493",
  "IMI-olefin" = "EX:E011214727", "3,5,6-TCP" = "EX:E000021564",
  "V" = "EX:E000022460", "Fe" = "EX:E000022401", "Co" = "EX:E000096046",
  "Cu" = "EX:E000022448", "Zn" = "EX:E000022464", "As" = "EX:E004815849",
  "Se" = "EX:E005270968", "Sr" = "EX:E004815779", "Mo" = "EX:E000022407",
  "Cd" = "EX:E000022444", "Te" = "EX:E005271068", "Tl" = "EX:E004815817",
  "Pb" = "EX:E004811303"
)

# 70 creatinine-adjusted pollutants (12 chemical classes)
pollutants <- c(
  "SFM", "4AN6C2B", "SFSX", "TMTH", "SFTHI",
  "NOR", "ENOR", "CIP", "OFLX", "TET",
  "OTC", "AM", "DOXY", "IMI", "NIT",
  "ACE", "THM", "CLO", "THI", "IPP",
  "CYC", "DIN", "IMID", "DM-THM", "DN-IMI",
  "5-OH-IMI", "DM-ACE", "FP", "DMM", "CBDZ",
  "CPS-O", "DZN", "DDVP", "MALT", "3OH-CBF",
  "ATZ", "ATZ-DE-DSI-2HD", "BPF", "BPS",
  "BPC", "BPE", "BPP", "2-OHN", "1-OHN",
  "3-OHF", "2-OHF", "2-OHP", "1-OHP", "4-OHP",
  "1-OHPYR", "3-PBA", "4-F-3-PBA", "2-MB-3CA",
  "trans-DCCA", "cis-DCCA", "IMI-olefin", "3,5,6-TCP",
  "V", "Fe", "Co", "Cu", "Zn", "As", "Se",
  "Sr", "Mo", "Cd", "Te", "Tl", "Pb"
)
pollutants_data_cols <- pollutants

# Baseline covariates
covariates_cat <- c("Gender", "Exercisefrequency", "Eatinghabits", "Smoke", "Drink")
covariates_cont <- c("Leftfastingbloodglucose2mmol", "BMI", "Age2")
baseline_vars <- c(covariates_cont, covariates_cat)

# Seven fasting clinical indicators
clinical_vars <- c("Platelet2gl", "WBC2gl", "Hemoglobin2gL", "ALB", "GLO", "UA", "HCY")

stopifnot(all(c(VarY, baseline_vars) %in% colnames(eSet$Data)))
clinical_vars <- intersect(clinical_vars, colnames(eSet$Data))

# -----------------------------------------------------------------------------
# Preprocessing
# -----------------------------------------------------------------------------

# Keep variables with < 33% missingness
eSet$Data <- eSet$Data %>%
  select(all_of(c(VarY, pollutants_data_cols, baseline_vars, clinical_vars))) %>%
  naniar::miss_var_summary() %>%
  filter(pct_miss < 66.67) %>%
  pull(variable) %>%
  { select(eSet$Data, all_of(.)) }

# Outcome as factor
eSet$Data <- eSet$Data %>%
  mutate(!!VarY := as.factor(!!sym(VarY)))

# One-hot encoding of categorical covariates
cat_vars <- intersect(covariates_cat, colnames(eSet$Data))
if (length(cat_vars) > 0) {
  dummy <- caret::dummyVars(~ ., data = eSet$Data[, cat_vars, drop = FALSE], fullRank = TRUE)
  dummy_df <- predict(dummy, newdata = eSet$Data[, cat_vars, drop = FALSE]) %>% as.data.frame()
  colnames(dummy_df) <- make.names(colnames(dummy_df), unique = TRUE)
  colnames(dummy_df) <- gsub("\\.", "_", colnames(dummy_df))
  eSet$Data <- eSet$Data %>%
    select(-all_of(cat_vars)) %>%
    bind_cols(dummy_df)
  baseline_vars <- c(covariates_cont, colnames(dummy_df))
}

# Scale continuous variables to [0, 1]
scale_to_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

pollutants_data_cols <- intersect(pollutants_data_cols, colnames(eSet$Data))

cont_vars <- intersect(c(covariates_cont, pollutants_data_cols, clinical_vars), colnames(eSet$Data))
eSet$Data <- eSet$Data %>%
  mutate(across(all_of(cont_vars), ~ scale_to_01(.x)))

# -----------------------------------------------------------------------------
# Feature engineering: Gene feature activity (exposure matrix %*% Gene matrix)
# -----------------------------------------------------------------------------
existing_pollutant_cols <- intersect(pollutants_data_cols, colnames(eSet$Data))
existing_exids <- pollutant_mapping[existing_pollutant_cols]
valid_mapping <- !is.na(existing_exids)
existing_pollutant_cols <- existing_pollutant_cols[valid_mapping]
existing_exids <- existing_exids[valid_mapping]

available_exids <- intersect(existing_exids, rownames(gene_matrix_raw))

gene_features <- character(0)
if (length(available_exids) > 0) {
  cols_to_use <- existing_pollutant_cols[existing_exids %in% available_exids]

  P <- eSet$Data %>% select(all_of(cols_to_use)) %>% as.matrix()
  B <- gene_matrix_raw[available_exids, , drop = FALSE] %>% as.matrix()
  B[is.na(B)] <- 0
  T_mat <- P %*% B
  colnames(T_mat) <- colnames(B)
  T_scaled <- apply(T_mat, 2, scale_to_01) %>%
    as_tibble() %>%
    rename_with(~ paste0("Gene_", .x))
  eSet$Data <- bind_cols(eSet$Data, T_scaled)
  gene_features <- colnames(T_scaled)
} else {
  warning("No pollutants could be matched to the Gene matrix; Gene features skipped.")
}

# -----------------------------------------------------------------------------
# Elastic-net screening of Gene features
# (alpha = 0.5, 10-fold CV, lambda.1se scaled by lambda_factor)
# -----------------------------------------------------------------------------
select_features_elasticnet <- function(features, data, target,
                                       alpha = 0.5, nfolds = 10,
                                       use_1se = TRUE,
                                       lambda_factor = 0.01,
                                       seed = 20240821) {
  cat("\n===== Elastic-net screening (alpha=0.5, 10-fold CV, lambda.1se, factor=",
      lambda_factor, ") =====\n")
  df <- data %>% select(all_of(c(target, features))) %>% drop_na()
  if (ncol(df) < 2 || nrow(df) == 0) {
    warning("Insufficient data; no Gene features returned.")
    return(character(0))
  }

  x <- as.matrix(df[, features])
  y <- as.numeric(df[[target]]) - 1
  set.seed(seed)
  cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = alpha, nfolds = nfolds,
                      type.measure = "class", parallel = FALSE)
  lambda_base <- if (use_1se) cv_fit$lambda.1se else cv_fit$lambda.min
  lambda_choice <- lambda_base * lambda_factor
  coef_vec <- as.matrix(coef(cv_fit, s = lambda_choice))
  selected <- rownames(coef_vec)[coef_vec[, 1] != 0 & rownames(coef_vec) != "(Intercept)"]
  cat("Non-zero feature count =", length(selected), "\n")
  if (length(selected) == 0) {
    warning("No Gene feature selected at the current lambda_factor.")
  }
  return(selected)
}

if (length(gene_features) > 0) {
  gene_features <- select_features_elasticnet(
    features = gene_features,
    data = eSet$Data,
    target = VarY,
    alpha = 0.5,
    nfolds = 10,
    use_1se = TRUE,
    lambda_factor = 0.01,
    seed = seed
  )
  cat("Gene feature count after screening:", length(gene_features), "\n")
} else {
  cat("No Gene features; elastic-net screening skipped.\n")
}

# -----------------------------------------------------------------------------
# Prespecified sensitivity analysis (2026-09-17):
# the 153 deep-pathway protein list replaces the elastic-net-selected genes.
# fixed_153 <- read.csv("protein_selection_153.csv")
# gene_features <- intersect(paste0("Gene_", fixed_153$feature), colnames(eSet$Data))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Feature groups
# -----------------------------------------------------------------------------
groups <- list(
  Baseline = baseline_vars,
  Clinical = clinical_vars,
  Exposome = c(existing_pollutant_cols, gene_features)
)
groups <- groups[lengths(groups) > 0]

# -----------------------------------------------------------------------------
# Stage 1: recursive feature screening per group (XGBoost, 5-fold CV)
# -----------------------------------------------------------------------------
select_features <- function(group_name, features, data, target, seed = 20240821) {
  cat("\n===== Screening group:", group_name, "=====\n")

  df <- data %>% select(all_of(c(target, features))) %>% drop_na()
  if (ncol(df) < 2) {
    warning("Group ", group_name, " has too few features; skipped.")
    return(NULL)
  }

  task <- as_task_classif(df, target = target, id = group_name)

  set.seed(seed)
  rsmp <- rsmp("cv", folds = 5)
  rsmp$instantiate(task)

  learner <- lrn("classif.xgboost",
                 predict_type = "prob",
                 predict_sets = c("train", "test"),
                 eta = 0.1)

  at <- AutoTuner$new(
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("classif.acc"),
    search_space = ps(
      subsample = p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = p_dbl(lower = 0.5, upper = 1),
      max_depth = p_int(lower = 3, upper = 10),
      nrounds = p_int(lower = 2, upper = 100)
    ),
    terminator = trm("evals", n_evals = 50),
    tuner = tnr("random_search")
  )
  ddpcr::quiet(at$train(task))
  learner$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  learner$train(task)
  imp <- learner$importance()
  if (length(imp) == 0) {
    warning("Group ", group_name, ": importance empty; all features kept.")
    return(features)
  }

  imp_sorted <- sort(imp, decreasing = TRUE)
  feature_order <- names(imp_sorted)

  results <- list()
  for (k in seq_along(feature_order)) {
    feat_sub <- feature_order[1:k]
    task_sub <- as_task_classif(df[, c(target, feat_sub)], target = target)

    design <- benchmark_grid(
      tasks = task_sub,
      learners = learner$clone(),
      resamplings = rsmp
    )
    bmr <- benchmark(design, store_models = FALSE)

    measures <- list(
      msr("classif.acc", id = "acc"),
      msr("classif.precision", id = "precision"),
      msr("classif.sensitivity", id = "sensitivity"),
      msr("classif.specificity", id = "specificity"),
      msr("classif.npv", id = "npv")
    )
    scores <- bmr$score(measures)

    avg_score <- scores %>%
      summarise(across(c(acc, precision, sensitivity, specificity, npv),
                       ~ mean(.x, na.rm = TRUE))) %>%
      mutate(avg = (acc + precision + sensitivity + specificity + npv) / 5) %>%
      pull(avg)

    if (!is.na(avg_score)) {
      results[[k]] <- tibble(n_features = k, avg_score = avg_score)
    }
  }

  if (length(results) == 0) {
    warning("Group ", group_name, ": all combinations NA; all features kept.")
    return(features)
  }

  best <- bind_rows(results) %>%
    filter(avg_score == max(avg_score, na.rm = TRUE)) %>%
    slice(1)
  best_n <- best$n_features
  final_features <- feature_order[1:best_n]
  cat("Best feature count:", best_n, "\n")
  cat("Final features:", final_features, "\n")

  return(final_features)
}

# -----------------------------------------------------------------------------
# Stage 1 base models: extract out-of-fold predicted probabilities
# -----------------------------------------------------------------------------
train_base_model <- function(group_name, features, data, target, seed) {
  df <- data %>% select(all_of(c(target, features))) %>% drop_na()
  task <- as_task_classif(df, target = target)

  set.seed(seed)
  rsmp <- rsmp("cv", folds = 5)
  rsmp$instantiate(task)

  learner <- lrn("classif.xgboost",
                 predict_type = "prob",
                 predict_sets = c("train", "test"),
                 eta = 0.1)

  at <- AutoTuner$new(
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("classif.acc"),
    search_space = ps(
      subsample = p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = p_dbl(lower = 0.5, upper = 1),
      max_depth = p_int(lower = 3, upper = 10),
      nrounds = p_int(lower = 2, upper = 100)
    ),
    terminator = trm("evals", n_evals = 50),
    tuner = tnr("random_search")
  )
  ddpcr::quiet(at$train(task))
  learner$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  design <- benchmark_grid(
    tasks = task,
    learners = learner,
    resamplings = rsmp
  )
  bmr <- benchmark(design)

  preds <- map_dfr(1:5, function(i) {
    as.data.table(bmr, reassemble_learners = TRUE, convert_predictions = TRUE,
                  predict_sets = "test")$prediction[[i]] %>%
      as.data.table() %>%
      select(row_ids, prob.1) %>%
      set_names(c("row_id", "prob"))
  }) %>%
    arrange(row_id) %>%
    distinct(row_id, .keep_all = TRUE)

  prob_vec <- rep(NA, nrow(data))
  prob_vec[preds$row_id] <- preds$prob
  prob_vec
}

# Run stage-1 screening
eSet$FinalFeatures <- list()
for (g in names(groups)) {
  feats <- select_features(g, groups[[g]], eSet$Data, VarY, seed)
  eSet$FinalFeatures[[g]] <- feats
}
for (g in names(groups)) {
  if (is.null(eSet$FinalFeatures[[g]])) {
    warning("Group ", g, ": no features selected; all features used.")
    eSet$FinalFeatures[[g]] <- groups[[g]]
  }
}

# Run stage-1 base models
base_probs <- list()
for (g in names(groups)) {
  probs <- train_base_model(g, eSet$FinalFeatures[[g]], eSet$Data, VarY, seed)
  base_probs[[paste0("pred.1.", g)]] <- probs
}
base_features <- as_tibble(base_probs)

# -----------------------------------------------------------------------------
# Build SG model input
# -----------------------------------------------------------------------------
selected_raw <- list()
for (g in names(groups)) {
  feats <- eSet$FinalFeatures[[g]]
  if (length(feats) > 0) {
    df_temp <- eSet$Data %>% select(all_of(feats))
    colnames(df_temp) <- paste0(g, "_", colnames(df_temp))
    selected_raw[[g]] <- df_temp
  }
}
selected_raw <- bind_cols(selected_raw)

df_sg <- bind_cols(
  eSet$Data %>% select(all_of(VarY)),
  selected_raw,
  base_features
)

voca <- tibble(
  VarName = colnames(df_sg),
  VarLabel = colnames(df_sg),
  VarType = case_when(
    VarName == VarY ~ "Outcome",
    str_detect(VarName, "^pred\\.1\\.") ~ "BaseModel",
    str_detect(VarName, "^Baseline_") ~ "Baseline",
    str_detect(VarName, "^Clinical_") ~ "Clinical",
    str_detect(VarName, "^Exposome_") ~ "Exposome",
    TRUE ~ "Other"
  )
)

# -----------------------------------------------------------------------------
# Stage 2: SG meta-model (feature screening + final model)
# -----------------------------------------------------------------------------
run_sg_model <- function(df, voca, model_name) {
  cat("\n===== Running SG model:", model_name, "=====\n")
  df_clean <- df %>% drop_na()
  task <- as_task_classif(df_clean, target = VarY, id = model_name)

  set.seed(seed)
  rsmp <- rsmp("cv", folds = 5)
  rsmp$instantiate(task)

  learner <- lrn("classif.xgboost",
                 predict_type = "prob",
                 predict_sets = c("train", "test"),
                 eta = 0.1)
  at <- AutoTuner$new(
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = msr("classif.acc"),
    search_space = ps(
      subsample = p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = p_dbl(lower = 0.5, upper = 1),
      max_depth = p_int(lower = 3, upper = 10),
      nrounds = p_int(lower = 2, upper = 100)
    ),
    terminator = trm("evals", n_evals = 50),
    tuner = tnr("random_search")
  )
  ddpcr::quiet(at$train(task))
  learner$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  learner$train(task)
  imp <- learner$importance()
  feature_order <- names(sort(imp, decreasing = TRUE))

  results <- list()
  for (k in seq_along(feature_order)) {
    feat_sub <- feature_order[1:k]
    task_sub <- as_task_classif(df_clean[, c(VarY, feat_sub)], target = VarY)

    design <- benchmark_grid(
      tasks = task_sub,
      learners = learner$clone(),
      resamplings = rsmp
    )
    bmr <- benchmark(design, store_models = FALSE)

    measures <- list(
      msr("classif.acc", id = "acc"),
      msr("classif.precision", id = "precision"),
      msr("classif.sensitivity", id = "sensitivity"),
      msr("classif.specificity", id = "specificity"),
      msr("classif.npv", id = "npv")
    )
    scores <- bmr$score(measures)
    avg_scores <- scores %>%
      summarise(across(c(acc, precision, sensitivity, specificity, npv),
                       ~ mean(.x, na.rm = TRUE))) %>%
      mutate(avg = (acc + precision + sensitivity + specificity + npv) / 5)

    results[[k]] <- tibble(n_features = k, avg_score = avg_scores$avg)
  }

  best <- bind_rows(results) %>%
    filter(avg_score == max(avg_score, na.rm = TRUE)) %>%
    slice(1)
  best_n <- best$n_features
  final_features <- feature_order[1:best_n]
  cat("Final feature count:", best_n, "\n")
  cat("Final features:", final_features, "\n")

  task_final <- as_task_classif(df_clean[, c(VarY, final_features)], target = VarY)
  set.seed(seed)
  rsmp_final <- rsmp("cv", folds = 5)
  rsmp_final$instantiate(task_final)

  learner_final <- lrn("classif.xgboost",
                       predict_type = "prob",
                       predict_sets = c("train", "test"),
                       eta = 0.1)
  at_final <- AutoTuner$new(
    learner = learner_final,
    resampling = rsmp("cv", folds = 5),
    measure = msr("classif.acc"),
    search_space = ps(
      subsample = p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = p_dbl(lower = 0.5, upper = 1),
      max_depth = p_int(lower = 3, upper = 10),
      nrounds = p_int(lower = 2, upper = 100)
    ),
    terminator = trm("evals", n_evals = 50),
    tuner = tnr("random_search")
  )
  ddpcr::quiet(at_final$train(task_final))
  learner_final$param_set$values <- at_final$model$tuning_instance$result_learner_param_vals

  design_final <- benchmark_grid(
    tasks = task_final,
    learners = learner_final,
    resamplings = rsmp_final
  )
  bmr_final <- benchmark(design_final, store_models = TRUE)

  measures_final <- list(
    msr("classif.acc", id = "acc"),
    msr("classif.precision", id = "precision"),
    msr("classif.sensitivity", id = "sensitivity"),
    msr("classif.specificity", id = "specificity"),
    msr("classif.npv", id = "npv"),
    msr("classif.auc", id = "auc")
  )
  scores <- bmr_final$score(measures_final) %>%
    summarise(across(c(acc, precision, sensitivity, specificity, npv, auc),
                     ~ mean(.x, na.rm = TRUE))) %>%
    mutate(across(everything(), ~ round(.x, 4)))

  learner_final$train(task_final)
  imp_final <- learner_final$importance()
  imp_df <- tibble(
    Features = names(imp_final),
    Importance = imp_final,
    Importance_pct = Importance / sum(Importance),
    Importance_cum = cumsum(Importance_pct)
  ) %>%
    left_join(voca, by = c("Features" = "VarName"))

  return(list(
    scores = scores,
    imp_df = imp_df,
    final_features = final_features,
    bmr_final = bmr_final,
    task_final = task_final,
    rsmp_final = rsmp_final
  ))
}

res <- run_sg_model(df_sg, voca, "Final_Model_3Groups")

# -----------------------------------------------------------------------------
# Robustness grid: 24 fixed hyperparameter settings (max_depth x nrounds),
# pooled AUC across the same CV folds (paired comparisons across models)
# -----------------------------------------------------------------------------
rob_grid <- expand.grid(max_depth = c(2, 3, 4, 6, 8, 10), nrounds = c(10, 25, 50, 100))
rob_out <- vector("list", nrow(rob_grid))
for (ri in seq_len(nrow(rob_grid))) {
  learner_fixed <- lrn("classif.xgboost", predict_type = "prob",
                       eta = 0.1, subsample = 1, colsample_bytree = 1,
                       max_depth = rob_grid$max_depth[ri], nrounds = rob_grid$nrounds[ri])
  design_fixed <- benchmark_grid(tasks = res$task_final, learners = learner_fixed,
                                 resamplings = res$rsmp_final)
  bmr_fixed <- benchmark(design_fixed, store_models = FALSE)
  rrs <- bmr_fixed$resample_results$resample_result
  pooled <- do.call(rbind, lapply(rrs, function(rr) as.data.table(rr$prediction())))
  roc_obj <- pROC::roc(pooled$truth, pooled$prob.1, quiet = TRUE)
  rob_out[[ri]] <- data.frame(max_depth = rob_grid$max_depth[ri], nrounds = rob_grid$nrounds[ri],
                              pooled_auc = as.numeric(roc_obj$auc))
}
rob_df <- bind_rows(rob_out)
cat("Robustness grid completed:", nrow(rob_df), "settings\n")

# -----------------------------------------------------------------------------
# Model performance
# -----------------------------------------------------------------------------
print(res$scores)

# -----------------------------------------------------------------------------
# Base-model importances (used for network figures)
# -----------------------------------------------------------------------------
base_imp <- list()
for (g in names(groups)) {
  df_temp <- eSet$Data %>% select(all_of(c(VarY, eSet$FinalFeatures[[g]]))) %>% drop_na()
  task_temp <- as_task_classif(df_temp, target = VarY)

  set.seed(seed)
  rsmp_temp <- rsmp("cv", folds = 5)
  rsmp_temp$instantiate(task_temp)

  learner_temp <- lrn("classif.xgboost",
                      predict_type = "prob",
                      predict_sets = c("train", "test"),
                      eta = 0.1)
  at_temp <- AutoTuner$new(
    learner = learner_temp,
    resampling = rsmp_temp,
    measure = msr("classif.acc"),
    search_space = ps(
      subsample = p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = p_dbl(lower = 0.5, upper = 1),
      max_depth = p_int(lower = 3, upper = 10),
      nrounds = p_int(lower = 2, upper = 100)
    ),
    terminator = trm("evals", n_evals = 50),
    tuner = tnr("random_search")
  )
  ddpcr::quiet(at_temp$train(task_temp))
  learner_temp$param_set$values <- at_temp$model$tuning_instance$result_learner_param_vals
  learner_temp$train(task_temp)

  imp_temp <- learner_temp$importance()
  imp_df_temp <- tibble(
    Feature = names(imp_temp),
    Importance = imp_temp,
    Importance_pct = Importance / sum(Importance)
  ) %>% arrange(desc(Importance_pct))

  base_imp[[g]] <- imp_df_temp
}

# -----------------------------------------------------------------------------
# Feature importance plot
# -----------------------------------------------------------------------------
imp_sg <- res$imp_df %>%
  filter(!is.na(Importance)) %>%
  arrange(desc(Importance_pct)) %>%
  mutate(Features = fct_reorder(Features, Importance_pct))

imp_sg <- imp_sg %>%
  mutate(GroupDetail = case_when(
    VarType == "BaseModel" ~ "BaseModel",
    VarType == "Baseline" ~ "Baseline",
    VarType == "Clinical" ~ "Clinical",
    VarType == "Exposome" & str_detect(Features, "^Exposome_Gene_") ~ "GO",
    VarType == "Exposome" & !str_detect(Features, "^Exposome_Gene_") ~ "Pollutant",
    TRUE ~ "Other"
  ))

color_palette <- c(
  "BaseModel" = "#FF4500",
  "Baseline"  = "#A9A9A9",
  "Clinical"  = "#FF8C00",
  "Gene"      = "#F8D582",
  "Pollutant" = "#2E8B57",
  "Other"     = "#000000"
)

p_imp <- ggplot(imp_sg, aes(x = Importance_pct, y = Features, fill = GroupDetail)) +
  geom_col() +
  scale_fill_manual(values = color_palette, name = "Variable Type") +
  labs(x = "Relative Importance", y = NULL,
       title = "Feature Importance (SG Model with 3 Groups - Gene matrix)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("importance_epd.png", p_imp,
       width = 8, height = max(5, nrow(imp_sg) * 0.3))
