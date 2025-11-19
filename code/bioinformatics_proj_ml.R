#North American Breeding Bird Survey Data Exploration (2024 Release) + NOAA Weather DATA (POR:1997-2023)

setwd("C:/Users/kirst/Desktop/data/aves/na_bbs/data")


library(tidymodels)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)

#Data processing and splitting of state and BCR training dataset

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




run_classification_models <- function(data,
                                      outcome = "status",
                                      save_prefix = "model_output",
                                      save_csv = TRUE,
                                      save_rds = TRUE) {
  
  # Convert outcome to factor
  data[[outcome]] <- factor(data[[outcome]])
  
  # -------------------------------------------------------------------------
  # Split
  # -------------------------------------------------------------------------
  set.seed(123)
  split <- initial_split(data, prop = 0.8, strata = outcome)
  train <- training(split)
  test  <- testing(split)
  
  # -------------------------------------------------------------------------
  # Recipe
  # -------------------------------------------------------------------------
  rec <- recipe(as.formula(paste(outcome, "~ .")), data = train) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_predictors())
  
  # -------------------------------------------------------------------------
  # Models
  # -------------------------------------------------------------------------
  rf_mod <- rand_forest(
    mtry = tune(),
    trees = 500,
    min_n = tune()
  ) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")
  
  svm_mod <- svm_rbf(
    cost = tune(),
    rbf_sigma = tune()
  ) %>%
    set_engine("kernlab") %>%
    set_mode("classification")
  
  # -------------------------------------------------------------------------
  # Workflows
  # -------------------------------------------------------------------------
  rf_wf  <- workflow() %>% add_model(rf_mod)  %>% add_recipe(rec)
  svm_wf <- workflow() %>% add_model(svm_mod) %>% add_recipe(rec)
  
  # -------------------------------------------------------------------------
  # CV folds
  # -------------------------------------------------------------------------
  folds <- vfold_cv(train, v = 5, strata = outcome)
  
  # -------------------------------------------------------------------------
  # Grids
  # -------------------------------------------------------------------------
  rf_grid <- grid_space_filling(
    parameters(
      finalize(mtry(), train),
      min_n()
    ),
    size = 10
  )
  
  svm_grid <- grid_space_filling(
    parameters(
      cost(),
      rbf_sigma()
    ),
    size = 12
  )
  
  # -------------------------------------------------------------------------
  # Metrics
  # -------------------------------------------------------------------------
  metrics <- metric_set(roc_auc, accuracy, sensitivity, specificity)
  
  # -------------------------------------------------------------------------
  # Tuning
  # -------------------------------------------------------------------------
  rf_tuned <- tune_grid(
    rf_wf, resamples = folds, grid = rf_grid,
    metrics = metrics, control = control_grid(save_pred = TRUE)
  )
  
  svm_tuned <- tune_grid(
    svm_wf, resamples = folds, grid = svm_grid,
    metrics = metrics, control = control_grid(save_pred = TRUE)
  )
  
  # -------------------------------------------------------------------------
  # Select best by ROC-AUC
  # -------------------------------------------------------------------------
  rf_best  <- select_best(rf_tuned, "roc_auc")
  svm_best <- select_best(svm_tuned, "roc_auc")
  
  rf_final_wf  <- finalize_workflow(rf_wf, rf_best)
  svm_final_wf <- finalize_workflow(svm_wf, svm_best)
  
  # -------------------------------------------------------------------------
  # Fit final
  # -------------------------------------------------------------------------
  rf_fit  <- fit(rf_final_wf, data = train)
  svm_fit <- fit(svm_final_wf, data = train)
  
  # -------------------------------------------------------------------------
  # Predictions
  # -------------------------------------------------------------------------
  rf_preds <- predict(rf_fit, test, type = "prob") %>%
    bind_cols(predict(rf_fit, test)) %>%
    bind_cols(test %>% select(all_of(outcome)))
  
  svm_preds <- predict(svm_fit, test, type = "prob") %>%
    bind_cols(predict(svm_fit, test)) %>%
    bind_cols(test %>% select(all_of(outcome)))
  
  # -------------------------------------------------------------------------
  # Test metrics
  # -------------------------------------------------------------------------
  rf_results  <- metrics(rf_preds, truth = .data[[outcome]], .pred_class)
  svm_results <- metrics(svm_preds, truth = .data[[outcome]], .pred_class)
  
  comparison <- bind_rows(
    rf_results %>% mutate(model = "Random Forest"),
    svm_results %>% mutate(model = "SVM")
  ) %>%
    select(model, .metric, .estimate)
  
  # -------------------------------------------------------------------------
  # ROC curves
  # -------------------------------------------------------------------------
  rf_roc <- roc_curve(rf_preds, truth = .data[[outcome]], starts_with(".pred_"))
  svm_roc <- roc_curve(svm_preds, truth = .data[[outcome]], starts_with(".pred_"))
  
  p_rf_roc <- autoplot(rf_roc) + ggtitle("Random Forest ROC Curve")
  p_svm_roc <- autoplot(svm_roc) + ggtitle("SVM ROC Curve")
  
  # -------------------------------------------------------------------------
  # Variable importance
  # -------------------------------------------------------------------------
  rf_imp <- rf_fit %>%
    extract_fit_parsnip() %>%
    vip::vi()
  
  # -------------------------------------------------------------------------
  # Save outputs
  # -------------------------------------------------------------------------
  if (save_csv) {
    write_csv(comparison, paste0(save_prefix, "_comparison.csv"))
    write_csv(rf_imp, paste0(save_prefix, "_rf_variable_importance.csv"))
  }
  
  if (save_rds) {
    saveRDS(list(
      rf_fit = rf_fit,
      svm_fit = svm_fit,
      rf_preds = rf_preds,
      svm_preds = svm_preds,
      rf_importance = rf_imp,
      comparison = comparison
    ), paste0(save_prefix, "_models.rds"))
  }
  
  # -------------------------------------------------------------------------
  # Return everything
  # -------------------------------------------------------------------------
  return(list(
    rf_fit = rf_fit,
    svm_fit = svm_fit,
    rf_preds = rf_preds,
    svm_preds = svm_preds,
    comparison = comparison,
    rf_importance = rf_imp,
    rf_roc_plot = p_rf_roc,
    svm_roc_plot = p_svm_roc
  ))
}
