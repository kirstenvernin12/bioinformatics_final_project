# AI Use Log For Final Project

# TASK 1 - Proposal 
- Tool/model & version: ChatGPT-5
- What I asked for: To review idea for machine learning workflow.
- Snippet of prompt(s): "Is this a reasonable approach to train a machine learning model using supervised learning to detect trends (decreasing, increasing or stable) in bird abundance data and identify the environmental variables that are most predictive of these trends: First do a correlation analysis to determine which environmental variables are correlated with the trends in the training bird population data. Then use supervised learning thorugh a decision tree and random forest approach to train on one dataset, rank which environmental variables are the most prediction of the population trends and then test on a different dataset. Also to determine how well the model the is working, how should I determine this?"
- How I verified correctness (tests, sample data): Verfied the output by comparing to class notes and primary literature.

# TASK 2 - Code Translation
- Tool/model & version: ChatGPT-5
- What I asked for: To convert python colab code provided in class (the baseline classifiers with pipelines, confusion matrix, cross-validation & small grid search) to R code.
- Snippet of prompt(s): "convert python colab code to R code"
- How I verified correctness (tests, sample data): Verified it worked by running though R, and installing additional packages that were necesarry to run the code.

# TASK 3 - Data Processing and Creation of Training Datasets
- Tool/model & version: ChatGPT-5
- What I asked for: To write a function to join and filter both the NOAA weather data and the North American Breeding Bird Datasets (raw data consists of multiple csv files)
- Snippet of prompt(s): "is there a faster way to join a 5824 obs df to a 6047800 obs data frame in R"; "the smaller df is the lookup table ad all of the fields should be joined to the matching route in the larger dataset. There will be multiple matching routes in the larger data frame."; "is there a way to take one data frame with latitude and longitude coordinates and match to closest latitude and longitude coordinates in another data frame and then export a list of the associated station number from this second data frame?"; can you explain what the crs is in this code: points <- st_as_sf(df1, coords = c("lon", "lat"), crs = 4326) stations <- st_as_sf(df2, coords = c("lon", "lat"), crs = 4326)"; "Is there a way to read in a large zip file add column names and filter it for a list of stations and then for the variables that you are interested in and then match the data based on station ID and Date"; "will fread work for .csv.gz files?"; "convert from long to wide format using ELEMENT field as the source for the new field names and populate with VALUE field"; "can you modify this code to get all of the zip files that range from the years 1997 to 2025 and then do all of the following processing steps to the weather data:"; "how to pull the Adjusted R-Squared value from summary (lm)"
- How I verified correctness (tests, sample data): I initially wrote an R script that did what I wanted the functions to do and then checked that the function provided the same output. I also verified the number of rows based on the sum of the number of rows in the smaller data frames when a full join was performed.

# TASK 4 - Writing the Function for Data Pre-processing, Training the ML Model, and Generating Outputs for Evaluation
- Tool/model & version: ChatGPT-5
- What I asked for: To understand and write a function that did what the code that was translated from Python did (result from Task #2 above).
- Snippet of prompt(s): "do you need to set.seed(1234) if you are using your own data rather than simulated data"; "using machine learning is this referencing the labele or classifiers? # The column name of the class/target variable. Change as needed. target_col <- "Condition" # e.g., "Condition" or whatever your outcome column is # If target is numeric but represents classes, set this to TRUE target_is_factor <- TRUE"; "Yes, can you please recommend how to structure the dataset, as well as modify this R code for the strategy that treats the labels seaprately (runs two models)."; "is it okay to have extra fields in your csv file that will not be used in training? or should I remove them?"; Where would I put the predictor variables in this code. The fields for the predictor variables are "RecordedCar", "TotalCarObs","NoiseDetected" if some of these are NULL values is that ok": "should I change all NULL to NA?"
- How I verified correctness (tests, sample data): Ran the function and if I got an error asked Chat what the error meant. See Task 5. 

# TASK 5 - Debugging the Function for Data Pre-processing, Training the ML Model, and Generating Outputs for Evaluation
- Tool/model & version: ChatGPT-5
- What I asked for: To help with debugging process for the function for data pre-processing, training the ML model and generating outputs
- Snippet of prompt(s): "what does this mean: → A | error: Error: Missing data in dependent variable. → B | warning: ! There are new levels in RouteName: "COLUMBUS", "FOUR LAKES", "OAK GROVE", "PEDRO CREEK", "PADRE ISLAND PAIS", and "CHAMPION". ℹ Consider using step_novel() before step_dummy() to handle unseen values. There were issues with some computations A: x220 B: x55" [[This was repeated many times with different error codes.]]
- How I verified correctness (tests, sample data): Ran the updated function, and if I got another error, asked Chat what it meant. I also re-examined the structure of the data and revisited the processing and model pre-processing steps to determine if the underlying structure was the cause of the errors. 
  
# TASK 6 - Data Reprocessing to Create New Training Dataset
- Tool/model & version: ChatGPT-5
- What I asked for: To write a function for data processing steps to get the slope of the linear regression line (representative of abundance trend over time)
- Snippet of prompt(s): "Write a function that dies the following for all species (scientific_name) for the attached dataset..."
- How I verified correctness (test, sample data): Ensured that the output of the function was the same that I got from my previously written code that did the same thing.
