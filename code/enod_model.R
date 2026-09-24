# =============================================================================
# Eno-D model (baseline covariates + clinical indicators only, no pollutants)
# Stacked generalisation pipeline for hypertension classification
# =============================================================================

rm(list = ls())

# Load packages
pacman::p_load(
  tidyverse, readxl, naniar, R6, rlist, mlr3verse, mlr3tuning, mlr3pipelines,
  DALEX, DALEXtra, ddpcr, caret, xgboost, ggplot2, patchwork, ggsci, vip, future,
  furrr, writexl, mice, pROC
)

seed <- 20260215
set.seed(seed)

# -----------------------------------------------------------------------------
# Input (load before running):
#   exposure data with baseline covariates, seven fasting clinical indicators
#   and the binary outcome "Hypertension"
# -----------------------------------------------------------------------------

eSet <- R6::R6Class("eSet", lock_class = FALSE, lock_objects = FALSE)
# eSet$Data <- ...   # load data here
# colnames(eSet$Data) <- make.names(colnames(eSet$Data), unique = TRUE)
# colnames(eSet$Data) <- gsub("\\.", "_", colnames(eSet$Data))

VarY <- "Hypertension"

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
  select(all_of(c(VarY, baseline_vars, clinical_vars))) %>%
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

cont_vars <- intersect(c(covariates_cont, clinical_vars), colnames(eSet$Data))
eSet$Data <- eSet$Data %>%
  mutate(across(all_of(cont_vars), ~ scale_to_01(.x)))

# Feature groups
groups <- list(
  Baseline = baseline_vars,
  Clinical = clinical_vars
)

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

res_BD <- run_sg_model(df_sg, voca, "Baseline_Clinical_Only")

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
  design_fixed <- benchmark_grid(tasks = res_BD$task_final, learners = learner_fixed,
                                 resamplings = res_BD$rsmp_final)
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
print(res_BD$scores)

# -----------------------------------------------------------------------------
# Base-model importances (used for network figures)
# -----------------------------------------------------------------------------
base_imp_BD <- list()
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

  base_imp_BD[[g]] <- imp_df_temp
}

# -----------------------------------------------------------------------------
# Feature importance plot
# -----------------------------------------------------------------------------
imp_bd <- res_BD$imp_df %>%
  filter(!is.na(Importance)) %>%
  arrange(desc(Importance_pct)) %>%
  mutate(Features = fct_reorder(Features, Importance_pct))

imp_bd <- imp_bd %>%
  mutate(GroupDetail = case_when(
    VarType == "BaseModel" ~ "BaseModel",
    VarType == "Baseline" ~ "Baseline",
    VarType == "Clinical" ~ "Clinical",
    TRUE ~ "Other"
  ))

color_palette_bd <- c(
  "BaseModel" = "#FF4500",
  "Baseline"  = "#A9A9A9",
  "Clinical"  = "#FF8C00",
  "Other"     = "#000000"
)

p_imp_bd <- ggplot(imp_bd, aes(x = Importance_pct, y = Features, fill = GroupDetail)) +
  geom_col() +
  scale_fill_manual(values = color_palette_bd, name = "Variable Type") +
  labs(x = "Relative Importance", y = NULL,
       title = "Feature Importance (Baseline + Clinical Only)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("importance_enod.png", p_imp_bd,
       width = 8, height = max(5, nrow(imp_bd) * 0.3))
