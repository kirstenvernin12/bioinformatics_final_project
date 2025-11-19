
# Bioinformatics_Project_converted.R
# Converted from Python Colab notebook to R (tidymodels + tidyverse)
# Assumes user will upload a CSV file and specify its path below.
# Author: ChatGPT (converted)
# Date: 2025-11-13

# --- Setup --------------------------------------------------------------------
# Install packages if needed (uncomment to install)
# install.packages(c("tidyverse", "tidymodels", "vip", "doParallel", "themis"))

library(tidyverse)
library(tidymodels)  # includes parsnip, recipes, workflows, yardstick, rsample, tune
library(vip)         # variable importance plotting (optional)
library(doParallel)  # parallel processing (optional)
# library(themis)    # for dealing with class imbalance (optional)

# Set a seed for reproducibility
set.seed(1234)

# --- User inputs --------------------------------------------------------------
# Replace with your uploaded CSV filename or path
data_path <- "data.csv"  # <-- change this to your CSV file

# The column name of the class/target variable. Change as needed.
target_col <- "Condition"  # e.g., "Condition" or whatever your outcome column is
# If target is numeric but represents classes, set this to TRUE
target_is_factor <- TRUE

# --- Load data ----------------------------------------------------------------
df <- readr::read_csv(data_path)

# Quick look
glimpse(df)
summary(df)

# Ensure target is a factor (classification use-case)
if (target_is_factor) {
  df <- df %>% mutate(!!target_col := as.factor(.data[[target_col]]))
}

# If there are obvious ID columns, you can remove them here (uncomment and edit)
# df <- df %>% select(-sample_id, -another_id)

# --- Train/Test Split --------------------------------------------------------
# Stratified split if classification
split <- initial_split(df, prop = 0.8, strata = all_of(target_col))
train_df <- training(split)
test_df  <- testing(split)

# --- Preprocessing recipe ----------------------------------------------------
# Build a recipe: impute (median for numeric), center/scale, handle categorical, and create any
# feature engineering steps you had in Python.
rec <- recipe(as.formula(paste(target_col, "~ .")), data = train_df) %>%
  # Remove near-zero variance predictors
  step_nzv(all_predictors()) %>%
  # Impute numeric predictors with median
  step_impute_median(all_numeric_predictors()) %>%
  # Impute categorical predictors with most frequent
  step_impute_mode(all_nominal_predictors()) %>%
  # Convert categorical to dummy variables
  step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
  # Center and scale numeric predictors (helps models like SVM and logistic)
  step_normalize(all_numeric_predictors())

# If you want to remove correlated predictors:
# rec <- rec %>% step_corr(all_numeric_predictors(), threshold = 0.9)

# Prep recipe and visualize (optional)
prep_rec <- prep(rec, training = train_df)
juice(prep_rec) %>% glimpse()

# --- Resampling / Cross-validation -------------------------------------------
# 5-fold cross-validation on training set (stratified)
cv_folds <- vfold_cv(train_df, v = 5, strata = all_of(target_col))

# --- Model specifications ----------------------------------------------------
# 1) Logistic regression
logreg_spec <- logistic_reg(mode = "classification", penalty = 0, mixture = 0) %>%
  set_engine("glm")

# 2) Random forest
rf_spec <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

# 3) Support Vector Machine (RBF)
svm_spec <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# --- Workflows ---------------------------------------------------------------
logreg_wf <- workflow() %>% add_model(logreg_spec) %>% add_recipe(rec)
rf_wf     <- workflow() %>% add_model(rf_spec)     %>% add_recipe(rec)
svm_wf    <- workflow() %>% add_model(svm_spec)    %>% add_recipe(rec)

# --- Tuning grids ------------------------------------------------------------
# Use reasonable ranges or let tune do the work
rf_grid <- grid_latin_hypercube(
  mtry(range = c(1L, ncol(train_df) - 1L)),
  min_n(range = c(2L, 10L)),
  size = 10
)

svm_grid <- grid_latin_hypercube(
  cost(range = c(-5, 2)),  # log2-scale in kernlab - tune will transform as needed
  rbf_sigma(range = c(-10, -1)),
  size = 10
)

# --- Parallel backend (optional) ---------------------------------------------
# Use all but one core
cores <- parallel::detectCores(logical = TRUE) - 1
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)

# --- Tuning (random forest and SVM) -----------------------------------------
# For speed, tune rf and svm; logistic regression no tuning here.
set.seed(234)
rf_tune_res <- tune_grid(
  rf_wf,
  resamples = cv_folds,
  grid = rf_grid,
  metrics = metric_set(accuracy, roc_auc, sens, spec),
  control = control_grid(save_pred = TRUE)
)

set.seed(345)
svm_tune_res <- tune_grid(
  svm_wf,
  resamples = cv_folds,
  grid = svm_grid,
  metrics = metric_set(accuracy, roc_auc, sens, spec),
  control = control_grid(save_pred = TRUE)
)

# Choose best hyperparameters by ROC AUC (change as desired)
best_rf <- select_best(rf_tune_res, "roc_auc")
best_svm <- select_best(svm_tune_res, "roc_auc")

# Finalize workflows with best parameters
rf_final_wf  <- finalize_workflow(rf_wf, best_rf)
svm_final_wf <- finalize_workflow(svm_wf, best_svm)

# Fit the final logistic (no tuning) on the training data
logreg_fit <- fit(logreg_wf, data = train_df)

# Fit finalized models on full training data
rf_fit  <- fit(rf_final_wf, data = train_df)
svm_fit <- fit(svm_final_wf, data = train_df)

# Stop parallel backend
stopCluster(cl)
registerDoSEQ()

# --- Evaluation on test set --------------------------------------------------
# Helper to get predictions + metrics for a fitted workflow
eval_model <- function(fit_obj, test_data, model_name) {
  preds <- predict(fit_obj, test_data, type = "prob") %>%
    bind_cols(predict(fit_obj, test_data, type = "class")) %>%
    bind_cols(test_data %>% select(all_of(target_col)))
  colnames(preds)[ncol(preds)] <- target_col
  preds <- preds %>% rename(.pred_class = .pred_class)
  preds <- preds %>% mutate(.pred_positive = .data[[paste0(".pred_", levels(test_data[[target_col]])[2])]])
  # Note: the above assumes a two-class problem and that the second level is the "positive" class.
  metrics <- metric_set(accuracy, roc_auc, sens, spec)
  res <- metrics(preds, truth = !!sym(target_col), estimate = .pred_class, .pred = .pred_positive)
  return(list(preds = preds, metrics = res))
}

logreg_res <- eval_model(logreg_fit, test_df, "LogReg")
rf_res     <- eval_model(rf_fit, test_df, "RF")
svm_res    <- eval_model(svm_fit, test_df, "SVM")

# Combine metrics into a tibble for easy comparison
all_metrics <- bind_rows(
  logreg_res$metrics %>% mutate(model = "LogReg"),
  rf_res$metrics     %>% mutate(model = "RF"),
  svm_res$metrics    %>% mutate(model = "SVM")
) %>% select(model, .metric, .estimate)

print(all_metrics)

# --- Confusion matrix and plots ----------------------------------------------
# Confusion matrix function
plot_conf_mat <- function(preds_df, model_name) {
  cm <- conf_mat(preds_df, truth = !!sym(target_col), estimate = .pred_class)
  autoplot(cm, type = "heatmap") + ggtitle(paste("Confusion Matrix:", model_name))
}

# Plot for each model
plot_conf_mat(logreg_res$preds, "LogReg")
plot_conf_mat(rf_res$preds, "RF")
plot_conf_mat(svm_res$preds, "SVM")

# --- ROC curves --------------------------------------------------------------
# ROC plotting helper (assumes binary classification)
plot_roc_curve <- function(fit_obj, test_data, model_name) {
  probs <- predict(fit_obj, test_data, type = "prob") %>% bind_cols(test_data %>% select(all_of(target_col)))
  # rename positive probability column
  positive_class <- levels(test_data[[target_col]])[2]
  prob_col <- paste0(".pred_", positive_class)
  roc_data <- roc_curve(probs, truth = !!sym(target_col), !!sym(prob_col))
  autoplot(roc_data) + ggtitle(paste("ROC curve -", model_name))
}

plot_roc_curve(logreg_fit, test_df, "LogReg")
plot_roc_curve(rf_fit, test_df, "RF")
plot_roc_curve(svm_fit, test_df, "SVM")

# --- Variable importance (for RF) -------------------------------------------
if ("ranger" %in% .packages(all.available = TRUE)) {
  try({
    rf_vip <- vip(rf_fit$fit$fit, num_features = 20) + ggtitle("RF Variable Importance")
    print(rf_vip)
  }, silent = TRUE)
}

# --- Save fitted models and results -----------------------------------------
# Save models as RDS if you want to reload later
saveRDS(logreg_fit, "logreg_fit.rds")
saveRDS(rf_fit, "rf_fit.rds")
saveRDS(svm_fit, "svm_fit.rds")
write_csv(all_metrics, "model_metrics.csv")

# --- Notes / To adapt --------------------------------------------------------
# - If you have more than two classes, adjust ROC calculations accordingly (one-vs-all or macro-averaging).
# - If your data needs more specialized preprocessing (e.g., rare OTU filtering, log transforms), insert steps into the recipe.
# - For huge datasets, consider using vroom::vroom() for faster CSV reading or simplifying the recipe.
# - If your positive class is not the second level, either relevel the factor: df[[target_col]] <- relevel(df[[target_col]], ref = "negative") OR adjust how .pred_positive is selected.

# End of script
