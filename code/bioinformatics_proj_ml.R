#North American Breeding Bird Survey Data Exploration (2024 Release) + NOAA Weather DATA (POR:1997-2023)

setwd("C:/Users/kirst/Desktop/UTA/Fall 2025/Bioinformatics/Final Project/bioinformatics_final_project")


library(tidyverse)
library(tidymodels)
library(vip)
library(future)
library(doFuture)

# Will model each trend separately (at the BCR scale and at the State scale). Because of this, may want to add more species
# to the training data. This is the most rigorous approach for abundance trends at different spatial scales. By treating
# the labels separately:
#     1. Because the trends at different scales often differ in magnitude and direction. Which we know is the case for 
#     these data.
#     2. Each responds to distinct ecological and sampling processes - and our hypothesis is potentially management policies?
#     3. Model interpretation is more clear. 

# The column name of the class/target variable. Change as needed.
# e.g., "Condition" or whatever your outcome column is. This tells the model which column in 
# the training dataset that is the outcome you want to be able to predict. This is the label. 
# This line points to the label column, not the classifier. 

# The full modeling pipeline is wrapped inside a function so that it can be called twice - once for each target variable (state and BCR trend)
#This avoids code duplication and ensures identical preprocessing and tuning steps across both modeling tasks. 




run_classification_models <- function(
    data,
    outcome,
    label_columns = c("State_Trend", "BCR_Trend"),
    save_prefix = "model_output",
    save_csv = TRUE,
    save_rds = TRUE,
    parallel = TRUE,
    n_cores = max(1, parallel::detectCores() - 1)
) {
  # ------------------------- Prevent leakage ------------------------------
  cols_to_drop <- setdiff(label_columns, outcome)
  data <- data %>% select(-any_of(cols_to_drop))
  
  # ------------------------- Clean character fields -----------------------
  data <- data %>%
    mutate(
      across(
        where(is.character),
        ~ .x |>
          na_if("NULL") |>
          na_if("") |>
          stringr::str_trim() |>
          na_if("")
      )
    )
  
  # ------------------------- Clean outcome safely -------------------------
  # Work on a local vector to avoid scoping issues
  out_vec <- data[[outcome]] |> as.character()
  out_vec <- out_vec |> stringr::str_trim()
  out_vec <- na_if(out_vec, "NULL")
  out_vec <- na_if(out_vec, "")
  
  keep <- !is.na(out_vec)
  data <- data[keep, , drop = FALSE]
  out_vec <- out_vec[keep]
  
  data[[outcome]] <- factor(out_vec)
  
  if (length(levels(data[[outcome]])) != 3) {
    stop(
      "Outcome '", outcome,
      "' must have exactly 3 classes after cleaning. Found: ",
      paste(levels(data[[outcome]]), collapse = ", ")
    )
  }
  
  # ------------------------- Predictors -----------------------------------
  allowed_predictors <- c(
    "SpeciesTotal", "RecordedCar", "TotalCarObs", "NoiseDetected",
    "PRCP", "TMAX", "TMIN", "Avg_Temp"
  )
  
  missing <- setdiff(allowed_predictors, names(data))
  if (length(missing) > 0) {
    stop("Missing predictor columns: ", paste(missing, collapse = ", "))
  }
  
  data <- data %>%
    mutate(
      SpeciesTotal  = as.numeric(SpeciesTotal),
      TotalCarObs   = as.numeric(TotalCarObs),
      PRCP          = as.numeric(PRCP),
      TMAX          = as.numeric(TMAX),
      TMIN          = as.numeric(TMIN),
      Avg_Temp      = as.numeric(Avg_Temp),
      RecordedCar   = as.factor(RecordedCar),
      NoiseDetected = as.factor(NoiseDetected)
    )
  
  # ------------------------- Parallel backend -----------------------------
  if (parallel) {
    plan(multisession, workers = n_cores)
    registerDoFuture()
  }
  
  # ------------------------- Train/Test split -----------------------------
  set.seed(123)
  split <- initial_split(data, prop = 0.8, strata = all_of(outcome))
  train <- training(split)
  test  <- testing(split)
  
  # ------------------------- Recipe ---------------------------------------
  rec <- recipe(
    as.formula(paste(outcome, "~", paste(allowed_predictors, collapse = "+"))),
    data = train
  ) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors())
  
  prep(rec, training = train)
  
  # ------------------------- CV folds -------------------------------------
  folds <- vfold_cv(train, v = 5, strata = all_of(outcome))
  
  # ------------------------- Model specs ----------------------------------
  rf_mod <- rand_forest(
    mtry = tune(),
    min_n = tune(),
    trees = 500
  ) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")
  
  svm_mod <- svm_rbf(
    cost = tune(),
    rbf_sigma = tune()
  ) %>%
    set_engine("kernlab") %>%
    set_mode("classification")
  
  rf_wf  <- workflow() %>% add_recipe(rec) %>% add_model(rf_mod)
  svm_wf <- workflow() %>% add_recipe(rec) %>% add_model(svm_mod)
  
  # ------------------------- Grids ----------------------------------------
  final_mtry <- min(ncol(train) - 1, length(allowed_predictors))
  
  rf_grid <- grid_space_filling(
    parameters(
      mtry(range = c(1L, final_mtry)),
      min_n()
    ),
    size = 10
  )
  
  svm_grid <- grid_space_filling(
    parameters(cost(), rbf_sigma()),
    size = 10
  )
  
  # ------------------------- Metrics --------------------------------------
  metrics_class <- metric_set(accuracy, kap, mcc)
  ctrl <- control_grid(save_pred = TRUE)
  
  # ------------------------- Tuning ---------------------------------------
  rf_tuned <- tune_grid(
    rf_wf,
    resamples = folds,
    grid = rf_grid,
    metrics = metrics_class,
    control = ctrl
  )
  
  svm_tuned <- tune_grid(
    svm_wf,
    resamples = folds,
    grid = svm_grid,
    metrics = metrics_class,
    control = ctrl
  )
  
  # ------------------------- Select best (named arg!) ---------------------
  best_rf  <- select_best(rf_tuned, metric = "accuracy")
  best_svm <- select_best(svm_tuned, metric = "accuracy")
  
  rf_final  <- finalize_workflow(rf_wf, best_rf)
  svm_final <- finalize_workflow(svm_wf, best_svm)
  
  # ------------------------- Fit final models -----------------------------
  rf_fit  <- fit(rf_final, train)
  svm_fit <- fit(svm_final, train)
  
  # ------------------------- Predictions ----------------------------------
  rf_preds <- predict(rf_fit, test, type = "prob") %>%
    bind_cols(predict(rf_fit, test)) %>%
    bind_cols(test %>% select(all_of(outcome)))
  
  svm_preds <- predict(svm_fit, test, type = "prob") %>%
    bind_cols(predict(svm_fit, test)) %>%
    bind_cols(test %>% select(all_of(outcome)))
  
  # ------------------------- Class metrics --------------------------------
  rf_results  <- metrics_class(rf_preds, truth = .data[[outcome]], estimate = .pred_class)
  svm_results <- metrics_class(svm_preds, truth = .data[[outcome]], estimate = .pred_class)
  
  # ------------------------- Macro-weighted ROC-AUC -----------------------
  # Only probability columns (exclude .pred_class)
  rf_prob_cols <- grep("^\\.pred_", names(rf_preds), value = TRUE)
  rf_prob_cols <- setdiff(rf_prob_cols, ".pred_class")
  
  svm_prob_cols <- grep("^\\.pred_", names(svm_preds), value = TRUE)
  svm_prob_cols <- setdiff(svm_prob_cols, ".pred_class")
  
  rf_auc <- roc_auc(
    rf_preds,
    truth = !!sym(outcome),
    !!!rlang::syms(rf_prob_cols),
    estimator = "macro_weighted"
  )
  
  svm_auc <- roc_auc(
    svm_preds,
    truth = !!sym(outcome),
    !!!rlang::syms(svm_prob_cols),
    estimator = "macro_weighted"
  )
  
  comparison <- bind_rows(
    rf_results %>% mutate(model = "RandomForest"),
    svm_results %>% mutate(model = "SVM"),
    tibble(model = "RandomForest", .metric = "roc_auc_macro", .estimate = rf_auc$.estimate),
    tibble(model = "SVM",         .metric = "roc_auc_macro", .estimate = svm_auc$.estimate)
  )
  
  # ------------------------- Confusion matrices ---------------------------
  rf_cm  <- conf_mat(rf_preds, truth = .data[[outcome]], estimate = .pred_class)
  svm_cm <- conf_mat(svm_preds, truth = .data[[outcome]], estimate = .pred_class)
  
  # ------------------------- ROC curves -----------------------------------
  rf_roc_plot  <- autoplot(roc_curve(rf_preds, truth = .data[[outcome]], starts_with(".pred_"))) +
    ggtitle("Random Forest ROC") +
    theme_minimal(base_size = 14)
  
  svm_roc_plot <- autoplot(roc_curve(svm_preds, truth = .data[[outcome]], starts_with(".pred_"))) +
    ggtitle("SVM ROC") +
    theme_minimal(base_size = 14)
  
  # ------------------------- RF variable importance -----------------------
  rf_importance <- vi(extract_fit_parsnip(rf_fit))
  
  rf_imp_plot <- rf_importance %>%
    ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col() +
    coord_flip() +
    xlab("Predictor") +
    ylab("Importance") +
    ggtitle("Random Forest Variable Importance") +
    theme_minimal(base_size = 14)
  
  # ------------------------- Confusion matrix plots -----------------------
  rf_cm_plot  <- autoplot(rf_cm, type = "heatmap") +
    ggtitle("Random Forest Confusion Matrix") +
    theme_minimal(base_size = 14)
  
  svm_cm_plot <- autoplot(svm_cm, type = "heatmap") +
    ggtitle("SVM Confusion Matrix") +
    theme_minimal(base_size = 14)
  
  # ------------------------- Save numeric outputs -------------------------
  if (save_csv) {
    write_csv(comparison, paste0(save_prefix, "_comparison.csv"))
    write_csv(rf_importance, paste0(save_prefix, "_rf_importance.csv"))
  }
  
  if (save_rds) {
    saveRDS(list(
      rf_fit = rf_fit,
      svm_fit = svm_fit,
      comparison = comparison,
      rf_cm = rf_cm,
      svm_cm = svm_cm,
      rf_importance = rf_importance,
      rf_tuned = rf_tuned,
      svm_tuned = svm_tuned
    ), paste0(save_prefix, "_models.rds"))
  }
  
  # ------------------------- Optionally save plots ------------------------
  if (save_csv || save_rds) {
    ggsave(paste0(save_prefix, "_rf_ROC.png"), rf_roc_plot,  width = 7, height = 6)
    ggsave(paste0(save_prefix, "_svm_ROC.png"), svm_roc_plot, width = 7, height = 6)
    ggsave(paste0(save_prefix, "_rf_importance_plot.png"), rf_imp_plot, width = 7, height = 6)
    ggsave(paste0(save_prefix, "_rf_confusion_matrix.png"), rf_cm_plot,  width = 7, height = 6)
    ggsave(paste0(save_prefix, "_svm_confusion_matrix.png"), svm_cm_plot, width = 7, height = 6)
  }
  
  # ------------------------- Return everything ----------------------------
  return(list(
    rf_fit = rf_fit,
    svm_fit = svm_fit,
    comparison = comparison,
    rf_cm = rf_cm,
    svm_cm = svm_cm,
    rf_importance = rf_importance,
    rf_roc_plot = rf_roc_plot,
    svm_roc_plot = svm_roc_plot,
    rf_imp_plot = rf_imp_plot,
    rf_cm_plot = rf_cm_plot,
    svm_cm_plot = svm_cm_plot,
    rf_tuned = rf_tuned,
    svm_tuned = svm_tuned
  ))
}

# =========================================================================
# Load the training datasets
# =========================================================================
training_TX <- read.csv("training/training_TX.csv")
training_BCR <- read.csv("training/training_BCR.csv")

#----------------------------------------------------------------------------
# Call the function to run the model at different geographic scales
#----------------------------------------------------------------------------

#BCR 
res_BCR <- run_classification_models(
  data = training_BCR,
  outcome = "BCR_Trend",
  label_columns = c("State_Trend", "BCR_Trend"),
  save_prefix = "BCR_results",
  parallel = TRUE,
  n_cores = 6
)

#State
res_State <- run_classification_models(
  data = training_TX,
  outcome = "State_Trend",
  label_columns = c("State_Trend", "BCR_Trend"),
  save_prefix = "State_results",
  parallel = TRUE,
  n_cores = 6
)


#----------------------------------------------------------------------------
# Debug section - DO NOT RUN ONCE THE MODEL IS WORKING!
#----------------------------------------------------------------------------
run_classification_models_debug <- function(data,
                                            outcome,
                                            label_columns = c("State_Trend", "BCR_Trend"),
                                            save_prefix = "model_output",
                                            parallel = FALSE,
                                            n_cores = 1) {
  
  # Remove other trend labels
  cols_to_drop <- setdiff(label_columns, outcome)
  data <- data %>% select(-any_of(cols_to_drop))
  
  # Clean NULL
  data <- data %>% mutate(across(everything(), ~ ifelse(. %in% c("NULL", ""), NA, .)))
  
  # Factor outcome
  data[[outcome]] <- factor(data[[outcome]])
  
  # Predictors
  allowed_predictors <- c(
    "SpeciesTotal","RecordedCar","TotalCarObs","NoiseDetected",
    "PRCP","TMAX","TMIN","Avg_Temp"
  )
  
  data <- data %>% mutate(
    SpeciesTotal = as.numeric(SpeciesTotal),
    TotalCarObs  = as.numeric(TotalCarObs),
    PRCP         = as.numeric(PRCP),
    TMAX         = as.numeric(TMAX),
    TMIN         = as.numeric(TMIN),
    Avg_Temp     = as.numeric(Avg_Temp),
    RecordedCar  = as.factor(RecordedCar),
    NoiseDetected = as.factor(NoiseDetected)
  )
  
  # Split
  set.seed(123)
  split <- initial_split(data, prop = 0.8, strata = all_of(outcome))
  train <- training(split)
  test  <- testing(split)
  
  # Recipe
  rec <- recipe(as.formula(paste(outcome, "~", paste(allowed_predictors, collapse="+"))),
                data=train) %>%
    step_zv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors())
  
  prep(rec, training=train)
  
  # CV folds
  folds <- vfold_cv(train, v=5, strata = all_of(outcome))
  
  # Models
  rf_mod <- rand_forest(mtry=tune(), min_n=tune(), trees=200) %>%
    set_engine("ranger", importance="impurity") %>%
    set_mode("classification")
  
  svm_mod <- svm_rbf(cost=tune(), rbf_sigma=tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")
  
  rf_wf  <- workflow() %>% add_recipe(rec) %>% add_model(rf_mod)
  svm_wf <- workflow() %>% add_recipe(rec) %>% add_model(svm_mod)
  
  # Manual mtry range
  num_pred <- length(allowed_predictors)
  
  rf_grid <- grid_latin_hypercube(
    mtry(range=c(1L, num_pred)),
    min_n(),
    size = 6
  )
  
  svm_grid <- grid_latin_hypercube(
    cost(),
    rbf_sigma(),
    size = 6
  )
  
  metrics_class <- metric_set(accuracy, kap, mcc)
  
  # Tune BOTH models
  rf_tuned <- tune_grid(
    rf_wf, resamples = folds,
    grid = rf_grid,
    metrics = metrics_class,
    control = control_grid(save_pred = TRUE)
  )
  
  svm_tuned <- tune_grid(
    svm_wf, resamples = folds,
    grid = svm_grid,
    metrics = metrics_class,
    control = control_grid(save_pred = TRUE)
  )
  
  # Try selecting best
  best_rf <- tryCatch(select_best(rf_tuned, "accuracy"),
                      error=function(e) e)
  
  best_svm <- tryCatch(select_best(svm_tuned, "accuracy"),
                       error=function(e) e)
  
  return(list(
    rf_tuned = rf_tuned,
    svm_tuned = svm_tuned,
    best_rf = best_rf,
    best_svm = best_svm
  ))
}
debug_BCR <- run_classification_models_debug(
  training_BCR,
  outcome="BCR_Trend"
)
debug_TX <- run_classification_models_debug(
  training_TX,
  outcome="State_Trend"
)

show_notes(debug_BCR$rf_tuned)
show_notes(debug_BCR$svm_tuned)
collect_metrics(debug_BCR$rf_tuned)
collect_metrics(debug_BCR$rf_tuned)

