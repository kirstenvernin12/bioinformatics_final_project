#North American Breeding Bird Survey Data Exploration (2024 Release) + NOAA Weather DATA (POR:1997-2023)

setwd("C:/Users/kirst/Desktop/data/aves/na_bbs/data")
getwd()

library(tidyverse)
library(tidymodels)  # includes parsnip, recipes, workflows, yardstick, rsample, tune
library(vip)         # variable importance plotting (optional)
library(doParallel)  # parallel processing (optional)
library(themis)    # for dealing with class imbalance (optional)
library(ranger)
library(kernlab)

# Set a seed for reproducibility (ensures that the random component behaves the same way each time the code is run, which
# is what guarantees reproducibility. This is necessary for models like Random Forest.)
set.seed(1234)

# Load the training dataset
training <- read.csv("training_final.csv") #NOTE: ctrl+shift+alt+m to replace all objects with the same name
#Convert NULL to NA
training[training == "NULL"] <- NA

#Remove unneeded fields. This will open memory and increase speed because large datasets with many unused columns take more
#memory and longer processing times. Fewer columns also make preprocessing and feature selection steps easier to follow. 
#This also reduces the risk of accidental leakage. Columns not intended as predictors may inadvertently leak information
#about the target (if they contain IDs, dates, or future information - basically, better safe than sorry!) It is best
#practice to keep only the columns you intend to use as predictors or target when feeding data into preprocessing/recipe 
# steps. 

#state training dataset
training_TX <- training %>% filter(StateNum==83) %>% 
  select(RouteName, Active, Stratum, BCR, Month, Day, Year, SpeciesTotal, 
  RecordedCar, TotalCarObs, NoiseDetected, English_Common_Name, Order, Family, 
  scientific_name, State_Trend, BCR_trend, PRCP, TMAX, TMIN)
  

#BCR training dataset
training_BCR <- training %>% 
  select(RouteName, Active, Stratum, BCR, Month, Day, Year, SpeciesTotal, 
         RecordedCar, TotalCarObs, NoiseDetected, English_Common_Name, Order, Family, 
         scientific_name, State_Trend, BCR_trend, PRCP, TMAX, TMIN)

#Remove the one BCR that is NA
training_BCR <- training_BCR  %>% filter(!is.na(BCR_trend))

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

run_model_for_target <- function(df, target_col, target_is_factor = TRUE) {
  
  # Ensure target is factor for classification
  if (target_is_factor) {
    df <- df %>% mutate(!!target_col := as.factor(.data[[target_col]]))
  }
  
  # Train/Test split (stratified)
  split <- initial_split(df, prop = 0.8, strata = all_of(target_col))
  train_df <- training(split)
  test_df  <- testing(split)
  
  # ------------------ Preprocessing recipe --------------------------------
  rec <- recipe(as.formula(paste(target_col, "~ .")), data = train_df) %>%
    step_nzv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_novel(all_nominal_predictors())%>% #handle unseen factor levels
    step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
    step_zv(all_predictors()) %>%           # remove zero-variance columns
    step_normalize(all_numeric_predictors())
  
  prep_rec <- prep(rec, training = train_df)
  
  # ------------------ Cross-validation ------------------------------------
  cv_folds <- vfold_cv(train_df, v = 5, strata = all_of(target_col))
  
  # ------------------ Model specifications --------------------------------
  logreg_spec <- multinom_reg(mode = "classification") %>%
    set_engine("nnet")
  
  rf_spec <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")
  
  svm_spec <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")
  
  # ------------------ Workflows ------------------------------------------
  logreg_wf <- workflow() %>% add_model(logreg_spec) %>% add_recipe(rec)
  rf_wf     <- workflow() %>% add_model(rf_spec)     %>% add_recipe(rec)
  svm_wf    <- workflow() %>% add_model(svm_spec)    %>% add_recipe(rec)
  
  # ------------------ Tuning grids ----------------------------------------
  rf_grid <- grid_latin_hypercube(
    mtry(range = c(1L, ncol(train_df) - 1L)),
    min_n(range = c(2L, 10L)),
    size = 10
  )
  
  svm_grid <- grid_latin_hypercube(
    cost(range = c(-5, 2)),
    rbf_sigma(range = c(-10, -1)),
    size = 10
  )
  
  # ------------------ Parallel backend ------------------------------------
  cores <- parallel::detectCores(logical = TRUE) - 1
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  
  # ------------------ Hyperparameter tuning -------------------------------
  set.seed(234)
  rf_tune_res <- tune_grid(
    rf_wf,
    resamples = cv_folds,
    grid = rf_grid,
    metrics = metric_set(accuracy, roc_auc,kap,mcc),
    control = control_grid(save_pred = TRUE)
  )
  
  set.seed(345)
  svm_tune_res <- tune_grid(
    svm_wf,
    resamples = cv_folds,
    grid = svm_grid,
    metrics = metric_set(accuracy, roc_auc,kap,mcc),
    control = control_grid(save_pred = TRUE)
  )
  
  best_rf  <- select_best(rf_tune_res, metric = "roc_auc")
  best_svm <- select_best(svm_tune_res, metric = "roc_auc")
  
  rf_final_wf  <- finalize_workflow(rf_wf, best_rf)
  svm_final_wf <- finalize_workflow(svm_wf, best_svm)
  
  # ------------------ Fit tuned models ------------------------------------
  logreg_fit <- fit(logreg_wf, data = train_df)
  rf_fit     <- fit(rf_final_wf, data = train_df)
  svm_fit    <- fit(svm_final_wf, data = train_df)
  
  stopCluster(cl)
  registerDoSEQ()
  
  # ------------------ Evaluation ------------------------------------------
  eval_model <- function(fit_obj, test_data, model_name) {
    preds <- predict(fit_obj, test_data, type = "prob") %>%
      bind_cols(predict(fit_obj, test_data, type = "class")) %>%
      bind_cols(test_data %>% select(all_of(target_col)))
    colnames(preds)[ncol(preds)] <- target_col
    preds <- preds %>% rename(.pred_class = .pred_class)
    pos <- levels(test_data[[target_col]])[2]
    preds <- preds %>% mutate(.pred_positive = .data[[paste0(".pred_", pos)]])
    metrics <- metric_set(accuracy, roc_auc, sens, spec)
    res <- metrics(preds, truth = !!sym(target_col), estimate = .pred_class, .pred = .pred_positive)
    return(list(preds = preds, metrics = res))
  }
  
  logreg_res <- eval_model(logreg_fit, test_df, "LogReg")
  rf_res     <- eval_model(rf_fit, test_df, "RF")
  svm_res    <- eval_model(svm_fit, test_df, "SVM")
  
  all_metrics <- bind_rows(
    logreg_res$metrics %>% mutate(model = "LogReg"),
    rf_res$metrics     %>% mutate(model = "RF"),
    svm_res$metrics    %>% mutate(model = "SVM")
  ) %>% select(model, .metric, .estimate)
  
  return(list(
    target = target_col,
    logreg = logreg_res,
    rf     = rf_res,
    svm    = svm_res,
    metrics = all_metrics
  ))
}


# Run model for BCR trend
res_BCR <- run_model_for_target(training_BCR, target_col = "BCR_trend")

# Run model for state trend
res_State <- run_model_for_target(training_TX, target_col = "State_Trend")

# Inspect metrics
res_BCR$metrics
res_State$metrics