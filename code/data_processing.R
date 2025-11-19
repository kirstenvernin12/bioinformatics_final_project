#North American Breeding Bird Survey Data Exploration (2024 Release)
setwd("C:/Users/kirst/Desktop/data/aves/na_bbs")
getwd()

#load libraries
library(dplyr)
library(ggplot2)
library(sf)
library(data.table)
library(R.utils)
library(lubridate)

#Data processing
#join the US BBS data into one file
#join state data into one df

# Define the folder path
path <- "States/US/"

# Get all CSV file names in the folder
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSVs into one dataframe
bbs_us <- files %>%
  lapply(read.csv, stringsAsFactors = FALSE) %>%
  bind_rows()

# Routes data (includes BCRs)
routes <- read.csv("Routes.csv")
#Remove "CountryNum"
routes <- routes %>% select(-CountryNum)
#join bbs_us and routes so that this can be filtered
bbs_us <- merge(bbs_us,routes,by=c("StateNum","Route"))

#to get list of BCRs that enter TX
tx <- read.csv("texas_analysis_09182025/Texas.csv")
tx <- merge(routes,tx,by=c("StateNum","Route"))
BCR_tx <- unique(tx$BCR)

#Texas Workflow
bbs_tx_bcr <- bbs_us %>% filter(BCR %in% BCR_tx)
#filter for just state data for TX (StateNum=83)
#bbs_tx <- bbs_tx_bcr %>% filter(StateNum==83)

#Texas observations with SpeciesList to obtain species names and routes to gather route information
species <- read.csv("SpeciesList.csv")
dist <- read.csv("VehicleData.csv")
weather <- read.csv("Weather.csv") 

#remove some fields 
dist <- dist %>% 
  select(-CountryNum,-StateNum, -Route, -RPID, -Year) 

#Get total cars observed
dist$TotalCarObs <- rowSums(dist[paste0("Car",1:50)],na.rm=TRUE)

#Create NoiseDetected field (1=Y,0=N)
dist$NoiseDetected <- ifelse(
  rowSums(dist[paste0("Noise", 1:50)] == 1, na.rm = TRUE) > 0,
  1,
  0
)

#Remove Car and Noise stop data 1-50
dist <- dist[ , !names(dist) %in% c(paste0("Car", 1:50), paste0("Noise", 1:50))]

#Remove CountryNum, StateNum, Route, RPID, Year
datetime <- weather %>% 
  select(RouteDataID,Month,Day,StartTime,EndTime)

#No filtering required because "AOU" only shared field
bbs_tx_bcr <- merge(bbs_tx_bcr,species, by="AOU")
bbs_tx_bcr <- merge(dist,bbs_tx_bcr, by="RouteDataID")
bbs_tx_bcr <- merge(datetime,bbs_tx_bcr,by="RouteDataID")

#Create summary fields/tables and then clean up species (what to do with hybrids, ect. prior to plotting)
#Step 1 - create a combined Genus species field
bbs_tx_bcr <- bbs_tx_bcr %>% mutate(scientific_name = paste(Genus,Species, sep = " "))

#organize fields
bbs_tx_bcr <- bbs_tx_bcr[,c("RouteDataID", "Route", "RouteName", "RouteTypeID", "RouteTypeDetailID","Active",
                          "AOU", "CountryNum", "StateNum","RPID","Latitude", "Longitude", "Stratum", "BCR",
                          "Month", "Day", "Year","StartTime", "EndTime", "Count10", "Count20", "Count30",            
                          "Count40", "Count50", "StopTotal", "SpeciesTotal", "RecordedCar", "TotalCarObs",
                          "NoiseDetected", "Seq", "English_Common_Name", "French_Common_Name", "Order","Family",             
                          "Genus", "Species","scientific_name")] 

#subset for 1997-2023 and species of interest
training <- bbs_tx_bcr %>%
  filter(scientific_name== "Riparia riparia"|scientific_name=="Chaetura pelagica"|
           scientific_name=="Chordeiles minor") %>%
  filter(Year>=1997)

#Colorado Workflow - copy from TX and update with proper codes

#Additional Processing - training data only
# Get the abundance data for entire state and BCRs (these are the BCRs that intersect the state),
# Adjust the abundance over time for the number of routes over time (group by year)
# Then get the trend over time (slope of the linear regression)
# use that to create label (increasing - slope > 0.00, decreasing - slope < 0.00, stable - slope = 0.00)
#these labels will go in either of two new fields (StateTrend, BCRTrend) - will have to rejoin the tables and then export a new csv file

#from macroecology - "in performing correlation analyses it is important yo use only sites where a species occurs to calculate its average abundance" (exclude 0 values)

#Bank Swallow (All trends decreasing)
#State Trend - Bank Swallow
train_state <- training %>% filter(StateNum==83 & scientific_name== "Riparia riparia") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")
  
# Fit linear regression model
bank_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(bank_lm)["Year"] #-0.4987682 (Decreasing for the state)
summary(bank_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Riparia riparia"
State_Trend <- "Decreasing"

State_bank <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Bank Swallow
train_bcr <- training %>% filter(scientific_name== "Riparia riparia") %>%
  group_by(Year,BCR) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

models_by_bcr <- train_bcr %>%
  group_by(BCR) %>%
  group_modify(~ {
    m <- lm(adj_species_abundance ~ Year, data = .x)
    tibble(
      slope = coef(m)["Year"],
      intercept = coef(m)["(Intercept)"]
    )
  })

models_by_bcr

ggplot(train_bcr, aes(x = Year, y = adj_species_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ BCR, scales = "free_y") +
  labs(
    x = "Year",
    y = "Adjusted Species Abundance",
    title = "Linear Regression of Adjusted Species Abundance Over Time by BCR"
  )

BCR_trend <- c("Increasing", "Increasing", NA, "Decreasing","Increasing","Decreasing","Decreasing") #One with no trend
BCR <- c(18,19,20,21,35,36,37)
scientific_name <- c("Riparia riparia","Riparia riparia","Riparia riparia","Riparia riparia",
                     "Riparia riparia","Riparia riparia","Riparia riparia")

BCR_bank <- data.frame(BCR, BCR_trend,scientific_name)

#Chimney Swift
#State Trend - Chimney Swift
train_state <- training %>% filter(StateNum==83 & scientific_name== "Chaetura pelagica") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
chimney_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(chimney_lm)["Year"] #-0.1377858 (Decreasing for the state)
summary(chimney_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Chaetura pelagica"
State_Trend <- "Decreasing"

State_chimney <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Bank Swallow
train_bcr <- training %>% filter(scientific_name== "Chaetura pelagica") %>%
  group_by(Year,BCR) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

models_by_bcr <- train_bcr %>%
  group_by(BCR) %>%
  group_modify(~ {
    m <- lm(adj_species_abundance ~ Year, data = .x)
    tibble(
      slope = coef(m)["Year"],
      intercept = coef(m)["(Intercept)"]
    )
  })

models_by_bcr

ggplot(train_bcr, aes(x = Year, y = adj_species_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ BCR, scales = "free_y") +
  labs(
    x = "Year",
    y = "Adjusted Species Abundance",
    title = "Linear Regression of Adjusted Species Abundance Over Time by BCR"
  )
#If 0.00 then stable
BCR_trend <- c("Stable", "Decreasing", "Decreasing","Decreasing", "Decreasing", "Increasing", "Decreasing","Decreasing")
BCR <- c(18,19,20,21,25,35,36,37)
scientific_name <- c("Chaetura pelagica","Chaetura pelagica","Chaetura pelagica","Chaetura pelagica","Chaetura pelagica","Chaetura pelagica", "Chaetura pelagica","Chaetura pelagica")

BCR_chimney <- data.frame(BCR, BCR_trend,scientific_name)

#Common Nighthawk
#State Trend - Common Nighthawk
train_state <- training %>% filter(StateNum==83 & scientific_name== "Chordeiles minor") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
night_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(night_lm)["Year"] #-0.191436 (Decreasing for the state)
summary(night_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Chordeiles minor"
State_Trend <- "Decreasing"

State_night <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Common Nighthawk
train_bcr <- training %>% filter(scientific_name== "Chordeiles minor") %>%
  group_by(Year,BCR) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

models_by_bcr <- train_bcr %>%
  group_by(BCR) %>%
  group_modify(~ {
    m <- lm(adj_species_abundance ~ Year, data = .x)
    tibble(
      slope = coef(m)["Year"],
      intercept = coef(m)["(Intercept)"]
    )
  })

models_by_bcr

ggplot(train_bcr, aes(x = Year, y = adj_species_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ BCR, scales = "free_y") +
  labs(
    x = "Year",
    y = "Adjusted Species Abundance",
    title = "Linear Regression of Adjusted Species Abundance Over Time by BCR"
  )

BCR_trend <- c("Decreasing", "Increasing", "Decreasing","Decreasing", "Decreasing", "Decreasing","Decreasing","Decreasing")
BCR <- c(18,19,20,21,25,35,36,37)
scientific_name <- c("Chordeiles minor","Chordeiles minor","Chordeiles minor","Chordeiles minor","Chordeiles minor","Chordeiles minor", "Chordeiles minor","Chordeiles minor")

BCR_night <- data.frame(BCR, BCR_trend,scientific_name)

#join the tables to the full training dataset
State <- merge(State_bank,State_chimney, all=TRUE)
State <- merge(State,State_night, all=TRUE)

BCR <- merge(BCR_bank,BCR_chimney, all=TRUE)
BCR <- merge(BCR,BCR_night,all=TRUE)

training <- training %>% left_join(State, by=c("StateNum","scientific_name"))
training <- training %>% left_join(BCR,by=c("BCR","scientific_name"))

#create csv of the bbs data
write.csv(training,"data/bbs_texas_training.csv", row.names=FALSE)

training <- read.csv("data/bbs_texas_training.csv")

#Match weather stations with closest route point
weatherstations <- read.csv("data/NOAA_weather/ghcnd-stations.csv")
#First convert both data frames to sf point objects
points <- st_as_sf(training,
                   coords = c("Longitude", "Latitude"),
                   crs = 4326)

stations <- st_as_sf(weatherstations,
                     coords = c("LONGITUDE", "LATITUDE"),
                     crs = 4326)

#Find the closest weather station for each starting point for the BBS routes
nearest_index <- st_nearest_feature(points, stations)

#Extract the station numbers
points$NearestStation <- stations$ID[nearest_index]
weatherstations <- unique(points$NearestStation)

#Export the matching weather stations to know which stations to pull data for
write.csv(points, "data/nearest_stations_training.csv", row.names = FALSE)

#Import, filter, and process NOAA weather data

#Write function
process_weather_years <- function(
    path,
    years = 1997:2025,
    stations = NULL
) {
  # Build file list
  files <- sprintf("%s/%d.csv.gz", path, years)
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    stop("No .csv.gz files found for the requested years.")
  }
  
  # NOAA column names
  noaa_cols <- c(
    "ID", "DATE_YYYYMMDD", "ELEMENT", "VALUE",
    "MFLAG", "QFLAG", "SFLAG", "OBSTIME_HRMM"
  )
  
  # Function to process one year's file
  process_one <- function(file) {
    
    dt <- fread(file)
    setnames(dt, noaa_cols)
    
    setDT(dt)
    
    # Filter by station list only if provided
    if (!is.null(stations)) {
      dt <- dt[ID %in% stations]
    }
    
    # Keep only weather variables of interest
    dt <- dt[ELEMENT %in% c("PRCP", "TMAX", "TMIN")]
    
    # Keep essential columns
    dt <- dt[, .(ID, DATE_YYYYMMDD, ELEMENT, VALUE)]
    
    # Long → wide
    dt <- dcast(
      dt,
      ID + DATE_YYYYMMDD ~ ELEMENT,
      value.var = "VALUE"
    )
    
    # Convert units
    dt[, PRCP := PRCP / 10]   # tenths mm → mm
    dt[, TMAX := TMAX / 10]   # tenths °C → °C
    dt[, TMIN := TMIN / 10]   # tenths °C → °C
    
    return(dt)
  }
  
  # Apply to all years and combine
  result <- rbindlist(
    lapply(files, process_one),
    use.names = TRUE,
    fill = TRUE
  )
  
  return(result)
}

#Implement the function
path <- "data/NOAA_WEATHER/"
years <- 1997:2023
stations <- weatherstations

weather_processed <- process_weather_years(path, years, stations)


#BE SURE THAT THE DATE FIELDS ARE THE SAME IN BOTH DATA FRAMES WHEN PROCESSING TESTING DATA FOR CO!!!!! (see below)
#Create date and time fields in training data that are the same format as the weather data
points <- points %>% mutate(DATE_YYYYMMDD=sprintf("%04d%02d%02d",Year,Month,Day))
#make sure the same data type
weather_processed$DATE_YYYYMMDD <- as.character(weather_processed$DATE_YYYYMMDD)


#join the training data with the weather data - nearest weather collection station by match date and match station IDs
merged <- points %>%
  left_join(
    weather_processed,
    by = c("DATE_YYYYMMDD", "NearestStation" = "ID")
  )

#PRCP: mm
#TMIN and TMAX: degrees celsius
#re-write training dataset
write.csv(merged,"data/training_final.csv",row.names=FALSE)

