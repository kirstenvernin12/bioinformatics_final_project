#North American Breeding Bird Survey ML workflow - regression (predict slope)

setwd("C:/Users/kirst/Desktop/UTA/Fall 2025/Bioinformatics/Final Project/bioinformatics_final_project")


library(tidymodels)
library(workflows)
library(dplyr)
library(ggplot2)
library(vip)
library(furrr)
library(doFuture)
library(caret)

run_slope_regression <- function(
    data,
    drop_state_number = TRUE,
    save_prefix = "slope_model",
    parallel = TRUE,
    n_cores = max(1, parallel::detectCores() - 1)
) {
  
  # --------------------------------------------------------------------
  # 1. CLEANING
  # --------------------------------------------------------------------
  
  df <- data
  
  # Remove categorical ID fields
  drop_cols <- c(
    "scientific_name", "English_Common_Name", "Order", "Family",
    "Genus", "Species", "Route", "RouteName", "Stratum",
    "State", "X", "Unnamed..0", "Unnamed.0"
  )
  df <- df %>% select(-any_of(drop_cols))
  
  # Optionally drop StateNum
  if (drop_state_number) {
    df <- df %>% select(-any_of("StateNum"))
  }
  
  # Convert slope to numeric (required)
  df$slope <- suppressWarnings(as.numeric(df$slope))
  df <- df %>% filter(!is.na(slope))
  
  # Keep numeric columns only
  num_cols <- names(df)[sapply(df, is.numeric)]
  df <- df %>% select(all_of(num_cols))
  
  # Remove zero-variance predictors
  nzv <- nearZeroVar(df, names = TRUE)
  df <- df %>% select(-any_of(nzv))
  
  # --------------------------------------------------------------------
  # 2. TRAIN/TEST SPLIT
  # --------------------------------------------------------------------
  set.seed(123)
  split <- initial_split(df, prop = 0.8)
  train <- training(split)
  test  <- testing(split)
  
  # predictor list
  predictor_cols <- setdiff(names(train), "slope")
  
  # --------------------------------------------------------------------
  # 3. RECIPE
  # --------------------------------------------------------------------
  rec <- recipe(slope ~ ., data = train) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors())
  
  rec_prep <- prep(rec)
  
  # --------------------------------------------------------------------
  # 4. CROSS-VALIDATION
  # --------------------------------------------------------------------
  folds <- vfold_cv(train, v = 5)
  
  # --------------------------------------------------------------------
  # 5. MODEL SPECIFICATIONS
  # --------------------------------------------------------------------
  
  # Random Forest (regression)
  rf_mod <- rand_forest(
    mtry = tune(),
    min_n = tune(),
    trees = 500
  ) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("regression")
  
  # SVM RBF
  svm_mod <- svm_rbf(
    cost = tune(),
    rbf_sigma = tune()
  ) %>%
    set_engine("kernlab") %>%
    set_mode("regression")
  
  # Workflows
  rf_wf  <- workflow() %>% add_recipe(rec) %>% add_model(rf_mod)
  svm_wf <- workflow() %>% add_recipe(rec) %>% add_model(svm_mod)
  
  # --------------------------------------------------------------------
  # 6. GRID SETUP (*** FIX: Proper mtry bounds ***)
  # --------------------------------------------------------------------
  
  rf_params <- parameters(
    mtry(range = c(1L, length(predictor_cols))),  # HARD LIMIT
    min_n()
  )
  
  rf_grid <- grid_space_filling(
    rf_params,
    size = 20
  )
  
  svm_grid <- grid_space_filling(
    parameters(cost(), rbf_sigma()),
    size = 20
  )
  
  # --------------------------------------------------------------------
  # 7. PARALLEL BACKEND
  # --------------------------------------------------------------------
  if (parallel) {
    plan(multisession, workers = n_cores)
    registerDoFuture()
  }
  
  # --------------------------------------------------------------------
  # 8. TUNING
  # --------------------------------------------------------------------
  metrics_reg <- metric_set(rmse, rsq)
  
  ctrl <- control_grid(save_pred = TRUE, parallel_over = "resamples")
  
  rf_tuned <- tune_grid(
    rf_wf,
    resamples = folds,
    grid = rf_grid,
    metrics = metrics_reg,
    control = ctrl
  )
  
  svm_tuned <- tune_grid(
    svm_wf,
    resamples = folds,
    grid = svm_grid,
    metrics = metrics_reg,
    control = ctrl
  )
  
  # --------------------------------------------------------------------
  # 9. SELECT BEST MODELS
  # --------------------------------------------------------------------
  best_rf  <- select_best(rf_tuned, metric = "rmse")
  best_svm <- select_best(svm_tuned, metric = "rmse")
  
  rf_final  <- finalize_workflow(rf_wf, best_rf)
  svm_final <- finalize_workflow(svm_wf, best_svm)
  
  # --------------------------------------------------------------------
  # 10. FIT FINAL MODELS
  # --------------------------------------------------------------------
  rf_fit  <- fit(rf_final, train)
  svm_fit <- fit(svm_final, train)
  
  # --------------------------------------------------------------------
  # 11. PREDICTIONS + METRICS
  # --------------------------------------------------------------------
  rf_preds <- predict(rf_fit, test) %>%
    bind_cols(test %>% select(slope))
  
  svm_preds <- predict(svm_fit, test) %>%
    bind_cols(test %>% select(slope))
  
  rf_metrics  <- metrics_reg(rf_preds, truth = slope, estimate = .pred)
  svm_metrics <- metrics_reg(svm_preds, truth = slope, estimate = .pred)
  
  # --------------------------------------------------------------------
  # 12. VARIABLE IMPORTANCE
  # --------------------------------------------------------------------
  rf_imp <- vi(extract_fit_parsnip(rf_fit))
  
  # --------------------------------------------------------------------
  # 13. SAVE OUTPUT
  # --------------------------------------------------------------------
  write.csv(rf_metrics,  paste0(save_prefix, "_rf_metrics.csv"),  row.names = FALSE)
  write.csv(svm_metrics, paste0(save_prefix, "_svm_metrics.csv"), row.names = FALSE)
  write.csv(rf_imp,      paste0(save_prefix, "_rf_importance.csv"), row.names = FALSE)
  
  # --------------------------------------------------------------------
  # 14. RETURN RESULTS
  # --------------------------------------------------------------------
  return(list(
    rf_fit = rf_fit,
    svm_fit = svm_fit,
    rf_metrics = rf_metrics,
    svm_metrics = svm_metrics,
    rf_importance = rf_imp,
    rf_tuned = rf_tuned,
    svm_tuned = svm_tuned
  ))
}

res_State <- run_slope_regression(
  data = training_TX_clean,
  save_prefix = "State_slope",
  parallel = TRUE,
  n_cores = 6
)

res_BCR <- run_slope_regression(
  data = training_BCR_clean,
  save_prefix = "BCR_slope",
  parallel = TRUE,
  n_cores = 6
)
