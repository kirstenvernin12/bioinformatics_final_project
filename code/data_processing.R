#North American Breeding Bird Survey Data Exploration (2024 Release)
setwd("C:/Users/kirst/Desktop/UTA/Fall 2025/Bioinformatics/Final Project/bioinformatics_final_project")

#load libraries
library(dplyr)
library(ggplot2)
library(sf)
library(data.table)
library(R.utils)
library(lubridate)
library(purrr)
library(tidyr)
library(broom)
library(stringr)

#Data processing
#join the US BBS data into one file
#join state data into one df

# Define the folder path
path <- "raw_data/States/US/"

# Get all CSV file names in the folder
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSVs into one dataframe
bbs_us <- files %>%
  lapply(read.csv, stringsAsFactors = FALSE) %>%
  bind_rows()

# Routes data (includes BCRs)
routes <- read.csv("raw_data/Routes.csv")
#Remove "CountryNum"
routes <- routes %>% select(-CountryNum)
#join bbs_us and routes so that this can be filtered
bbs_us <- merge(bbs_us,routes,by=c("StateNum","Route"))

#to get list of BCRs that enter state of interest
get_state_bcrs <- function(data, state_number) {
  data %>%
    filter(StateNum == state_number) %>%
    distinct(BCR) %>%
    pull(BCR)
}

#Texas Workflow
tx_bcrs <- get_state_bcrs(bbs_us, 83)
bbs_tx_bcr <- bbs_us %>% filter(BCR %in% tx_bcrs)

#Texas observations with SpeciesList to obtain species names and routes to gather route information
species <- read.csv("raw_data/SpeciesList.csv")
dist <- read.csv("raw_data/VehicleData.csv")
weather <- read.csv("raw_data/Weather.csv") 

#remove some fields 
dist <- dist %>% 
  select(-CountryNum,-StateNum, -Route, -RPID, -Year) 

#Get total cars observed
dist$TotalCarObs <- rowSums(dist[paste0("Car",1:50)],na.rm=TRUE)

#Remove Car and Noise stop data 1-50
dist <- dist[ , !names(dist) %in% c(paste0("Car", 1:50), paste0("Noise", 1:50))]

#Remove CountryNum, StateNum, Route, RPID, Year
weather <- weather %>% 
  select(RouteDataID,Month,Day,StartTime,EndTime,StartTemp,EndTemp) 

weather$StartTemp <- as.numeric(weather$StartTemp)
weather$EndTemp <- as.numeric(weather$EndTemp)


weather <- weather %>%
  select(RouteDataID,Month,Day,StartTime,EndTime,StartTemp, EndTemp)

#No filtering required because "AOU" only shared field
bbs_tx_bcr <- merge(bbs_tx_bcr,species, by="AOU")
bbs_tx_bcr <- merge(dist,bbs_tx_bcr, by="RouteDataID")
bbs_tx_bcr <- merge(weather,bbs_tx_bcr,by="RouteDataID")

#Create summary fields/tables and then clean up species (what to do with hybrids, ect. prior to plotting)
#Step 1 - create a combined Genus species field
bbs_tx_bcr <- bbs_tx_bcr %>% mutate(scientific_name = paste(Genus,Species, sep = " "))

bbs_tx_bcr <- bbs_tx_bcr %>% select(c(RouteDataID,Route, AOU,CountryNum,StateNum, Latitude,Longitude,Stratum,
                                      BCR,Month,Day,Year,StartTemp, EndTemp,TotalCarObs,SpeciesTotal,English_Common_Name,Order,Family,
                                      scientific_name))%>%
  rename(Abundance=SpeciesTotal)


#organize fields
bbs_tx_bcr <- bbs_tx_bcr[,c("RouteDataID", "Route","AOU", "CountryNum", "StateNum","Latitude", "Longitude", "Stratum", "BCR",
                          "Month", "Day", "Year","StartTemp", "EndTemp", "Abundance", "TotalCarObs",
                          "English_Common_Name", "Order","Family","scientific_name")] 
write.csv(bbs_tx_bcr,"raw_data/bbs_tx_bcr.csv",row.names = FALSE)

#filter for period of interest (1997-2023)
bbs_tx_sub <- bbs_tx_bcr %>% filter(Year>=1997)


#Match weather stations with closest route point
weatherstations <- read.csv("raw_data/NOAA_weather/ghcnd-stations.csv")
#First convert both data frames to sf point objects
points <- st_as_sf(bbs_tx_sub,
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
path <- "raw_data/NOAA_WEATHER/"
years <- 1997:2023
stations <- weatherstations

weather_processed <- process_weather_years(path, years, stations)

#BE SURE THAT THE DATE FIELDS ARE THE SAME IN BOTH DATA FRAMES
#Create date and time fields in training data that are the same format as the weather data
points <- points %>% mutate(DATE_YYYYMMDD=sprintf("%04d%02d%02d",Year,Month,Day))
#make sure the same data type
weather_processed$DATE_YYYYMMDD <- as.character(weather_processed$DATE_YYYYMMDD)


#join the bbs data with the weather data - nearest weather collection station by match date and match station IDs
merged <- points %>%
  left_join(
    weather_processed,
    by = c("DATE_YYYYMMDD", "NearestStation" = "ID")
  )
#PRCP: mm
#TMIN and TMAX: degrees celsius
#StartTemp and EndTemp degrees F
#Convert daily data into annual summaries per route 
env <- as.data.frame(merged)

write.csv(env,"raw_data/bbs_bcr_tx_env.csv")

#Calculate annual summary stats

num_routes_state <- env %>% filter(StateNum==83) %>%
  group_by(Year)%>%
  summarise(num_routes = n_distinct(Route))
  
num_routes_BCR <- env %>%
  group_by(Year,BCR)%>%
  summarise(num_routes = n_distinct(Route))

#Species richness
#State
state_rich <- env %>% filter(StateNum==83) %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year) %>%
  summarise(species_rich=n_distinct(scientific_name)) 

state_numroutes <- env %>% filter(StateNum==83) %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year) %>%
  summarise(num_routes=n_distinct(Route)) 

state_rich_perroute <- env %>% filter(StateNum==83) %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year) %>%
  summarise(species_per_route=round(n_distinct(scientific_name)/n_distinct(Route)))

state_annual_summ <- merge(state_rich,state_numroutes,by="Year")
state_annual_summ <- merge(state_annual_summ,state_rich_perroute,by="Year")

#convert from long to wide format
state_rich <- state_annual_summ %>%
  pivot_longer(
    cols = c(species_rich, num_routes, species_per_route),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value,
    values_fill = 0
  )

#BCR
BCR_rich <- env %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year,BCR) %>%
  summarise(species_rich=n_distinct(scientific_name)) 

BCR_numroutes <- env %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year,BCR) %>%
  summarise(num_routes=n_distinct(Route)) 

BCR_numroutes <- env %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year,BCR,scientific_name) %>%
  summarise(num_routes=n_distinct(Route)) 

BCR_rich_perroute <- env %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year,BCR) %>%
  summarise(species_per_route=round(n_distinct(scientific_name)/n_distinct(Route)))

BCR_annual_summ <- merge(BCR_rich,BCR_numroutes,by=c("Year","BCR"))
BCR_annual_summ <- merge(BCR_annual_summ,BCR_rich_perroute,by=c("Year","BCR"))

#convert from long to wide format
BCR_rich <- BCR_annual_summ %>%
  pivot_longer(
    cols = c(species_rich, num_routes, species_per_route),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(BCR,feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value,
    values_fill = 0
  )

#Abundance
#State
state_abund <- env %>% filter(StateNum==83) %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year,scientific_name) %>%
  summarise(abundance=sum(Abundance)) 

state_numroutes <- state_annual_summ %>% select(c(Year,num_routes))

state_abund <- merge(state_abund,state_numroutes,by="Year")


state_abund_summ <-state_abund %>%
  group_by(Year,scientific_name) %>%
  summarise(abund_per_route=round(abundance/num_routes))

state_annual_abund_summ<- merge(state_abund_summ,state_abund,by=c("Year","scientific_name"))

#convert from long to wide format
state_abundance <- state_annual_abund_summ %>%
  pivot_longer(
    cols = c(abund_per_route, abundance),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(scientific_name, feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value,
    values_fill = 0
  )

#BCR
BCR_abund <- env %>% filter(StateNum==83) %>%
  filter(                                              #this filters out unknown species and hybrids
    !grepl(" / ", scientific_name),
    !grepl(" sp\\.", scientific_name),
    !grepl(" x ", scientific_name, ignore.case = TRUE)
  ) %>%
  group_by(Year,scientific_name,BCR) %>%
  summarise(abundance=sum(Abundance)) 

BCR_numroutes <- BCR_annual_summ %>% select(c(Year,BCR,num_routes))

BCR_abund <- merge(BCR_abund,BCR_numroutes,by=c("Year","BCR"))


BCR_abund_summ <-BCR_abund %>%
  group_by(Year,scientific_name,BCR) %>%
  summarise(abund_per_route=round(abundance/num_routes))

BCR_annual_abund_summ<- merge(BCR_abund_summ,BCR_abund,by=c("Year","scientific_name","BCR"))

#convert from long to wide format
BCR_abundance <- BCR_annual_abund_summ %>%
  pivot_longer(
    cols = c(abund_per_route, abundance),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(scientific_name, BCR,feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value,
    values_fill = 0
  )


#join abundance and richness
#State
state_abundance$StateNum <- 83
state_rich$StateNum <- 83

state_trends <- merge(state_abundance,state_rich,by="StateNum")

#BCR
bcr_trends <- merge(BCR_abundance,BCR_rich,by="BCR")

#calculate median values for environmental data based on spatial area of interest (state or BCR)
#State
env_med <- env %>% filter(StateNum==83) %>%
  group_by(Year)%>%
  summarise(Med_StartTemp=median(StartTemp),Med_EndTemp=median(EndTemp),,
            Med_CarsObs=median(TotalCarObs),Med_PRCP=median(PRCP[PRCP !=0],na.rm=TRUE),
            Med_TMAX=median(TMAX,na.rm=TRUE),Med_TMIN=median(TMIN,na.rm=TRUE))

#convert F to C
fahrenheit_to_celsius <- function(F) {
  ifelse(is.na(F), NA, (F - 32) * 5/9)
}

env_med <- env_med %>%
  mutate(
    across(c(Med_StartTemp, Med_EndTemp),
           ~ fahrenheit_to_celsius(.x),
           .names = "{.col}_C")
  )


env_med <- env_med %>% select(-Med_StartTemp,-Med_EndTemp)
  
env_med <- env_med %>%  
  mutate(Med_StartTemp=round(Med_StartTemp_C,1),Med_EndTemp=round(Med_EndTemp_C,1)) %>%
  select(-Med_StartTemp_C,-Med_EndTemp_C)

env_med <- env_med %>%
  pivot_longer(
    cols = c(Med_CarsObs, Med_PRCP,Med_TMAX,Med_TMIN,Med_StartTemp,Med_EndTemp),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value,
    values_fill = 0
  )

#merge with state data
env_med$StateNum <- 83

state_trends <- merge(state_trends,env_med,by="StateNum")

#BCR
env_med_b <- env %>% 
  group_by(Year,BCR)%>%
  summarise(Med_StartTemp=median(StartTemp),Med_EndTemp=median(EndTemp),,
            Med_CarsObs=median(TotalCarObs),Med_PRCP=median(PRCP[PRCP !=0],na.rm=TRUE),
            Med_TMAX=median(TMAX,na.rm=TRUE),Med_TMIN=median(TMIN,na.rm=TRUE))

env_med_b <- env_med_b %>%
  mutate(
    across(c(Med_StartTemp, Med_EndTemp),
           ~ fahrenheit_to_celsius(.x),
           .names = "{.col}_C")
  )

env_med_b <- env_med_b %>% select(-Med_StartTemp,-Med_EndTemp)

env_med_b <- env_med_b %>%  
  mutate(Med_StartTemp=round(Med_StartTemp_C,1),Med_EndTemp=round(Med_EndTemp_C,1)) %>%
  select(-Med_StartTemp_C,-Med_EndTemp_C)



env_med_b <- env_med_b %>%
  ungroup() %>%    
  pivot_longer(
    cols = c(Med_CarsObs, Med_PRCP, Med_TMAX, Med_TMIN, 
             Med_StartTemp, Med_EndTemp),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(BCR, feature, value) %>%
  pivot_wider(
    names_from  = feature,
    values_from = value,
    values_fill = 0
  )


#merge with bcr data
bcr_trends <- merge(bcr_trends,env_med_b,by="BCR")

# Trend calculations [training dataset only]
# Get the abundance data for entire state and BCRs (these are the BCRs that intersect the state),
# Adjust the abundance over time for the number of routes over time (group by year)
# Then get the trend over time (slope of the linear regression)
# use that to create label (increasing - slope > 0.00, decreasing - slope < 0.00, stable - slope = 0.00)
#these labels will go in either of two new fields (StateTrend, BCRTrend) - will have to rejoin the tables and 
#then export a new csv file

#from macroecology - "in performing correlation analyses it is important you use only sites where a species 
#occurs to calculate its average abundance" (exclude 0 values)

#filtered for the subset of interest for our study (1997-2023)


#slope calculations for 1997-2023
analyze_all_species <- function(data,
                                state_number = 83,
                                output_dir = "trends/texas",
                                slope_threshold = 0.04) {
  
  # Ensure directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ---------------------------------------------------
  # 0. FILTER SCIENTIFIC NAMES
  # ---------------------------------------------------
  data_clean <- data %>%
    filter(
      !grepl(" / ", scientific_name),
      !grepl(" sp\\.", scientific_name),
      !grepl(" x ", scientific_name, ignore.case = TRUE)
    )
  
  message("Filtered species count: ", length(unique(data_clean$scientific_name)))
  
  species_list <- sort(unique(data_clean$scientific_name))
  
  state_results <- list()
  BCR_master <- list()
  
  for (sp in species_list) {
    
    message("Processing: ", sp)
    
    # ---------------------------------------------------
    # 1. STATE-LEVEL TREND (FILTERED BY StateNum)
    # ---------------------------------------------------
    train_state <- data_clean %>%
      filter(StateNum == state_number,
             scientific_name == sp) %>%
      group_by(Year) %>%
      summarise(adj_species_abundance = round(sum(Abundance) /
                                                n_distinct(Route)), .groups = "drop")
    
    slope <- adjR2 <- lower_CI <- upper_CI <- NA
    trend <- NA
    
    if (nrow(train_state) >= 3) {
      mod_state <- lm(adj_species_abundance ~ Year, data = train_state)
      slope <- coef(mod_state)["Year"]
      adjR2 <- summary(mod_state)$adj.r.squared
      ci <- confint(mod_state)["Year", ]
      lower_CI <- ci[1]
      upper_CI <- ci[2]
      
      # Trend classification with slope threshold + CI override
      if (lower_CI <= 0 && upper_CI >= 0) {
        trend <- "Stable"
      } else if (slope >= slope_threshold) {
        trend <- "Increasing"
      } else if (slope <= -slope_threshold) {
        trend <- "Decreasing"
      } else {
        trend <- "Stable"
      }
    }
    
    state_results[[sp]] <- data.frame(
      scientific_name = sp,
      StateNum = state_number,
      slope = slope,
      lower_CI = lower_CI,
      upper_CI = upper_CI,
      Adj_R2 = adjR2,
      Trend = trend
    )
    
    
    # ---------------------------------------------------
    # 1b. SAFE STATE PLOT
    # ---------------------------------------------------
    if (nrow(train_state) > 0) {
      
      title_text <- paste0(
        "Linear Regression of Adjusted Species Abundance Over Time\n",
        "Slope = ", round(slope, 4),
        " | CI = [", round(lower_CI, 4), ", ", round(upper_CI, 4), "]",
        " | Adjusted R² = ", round(adjR2, 4),
        " | Trend: ", trend
      )
      
      png(filename = paste0(output_dir, "/", gsub(" ", "_", sp),
                            "_state_plot.png"),
          width = 1200, height = 900, res = 150)
      
      plot(train_state$Year, train_state$adj_species_abundance,
           xlab = "Year",
           ylab = "Adjusted Species Abundance (#Individuals/#Routes)",
           main = title_text)
      
      if (nrow(train_state) >= 3) abline(mod_state, lwd = 2)
      
      dev.off()
    }
    
    
    # ---------------------------------------------------
    # 2. BCR-LEVEL TRENDS (FULL DATA, NO STATE FILTER)
    # ---------------------------------------------------
    train_bcr <- data_clean %>%
      filter(scientific_name == sp) %>%
      group_by(Year, BCR) %>%
      summarise(adj_species_abundance = round(sum(Abundance) /
                                                n_distinct(Route)), .groups = "drop")
    
    if (nrow(train_bcr) == 0) next
    
    BCR_trends <- train_bcr %>%
      group_by(BCR) %>%
      group_modify(~{
        
        if (nrow(.x) >= 3) {
          m <- lm(adj_species_abundance ~ Year, data = .x)
          
          slope <- coef(m)["Year"]
          ci <- confint(m)["Year", ]
          lower_CI <- ci[1]
          upper_CI <- ci[2]
          
          # CI override + threshold rule
          if (lower_CI <= 0 && upper_CI >= 0) {
            trend <- "Stable"
          } else if (slope >= slope_threshold) {
            trend <- "Increasing"
          } else if (slope <= -slope_threshold) {
            trend <- "Decreasing"
          } else {
            trend <- "Stable"
          }
          
          tibble(
            slope = slope,
            lower_CI = lower_CI,
            upper_CI = upper_CI,
            intercept = coef(m)["(Intercept)"],
            adj_r2 = summary(m)$adj.r.squared,
            Trend = trend
          )
          
        } else {
          tibble(
            slope = NA, lower_CI = NA, upper_CI = NA,
            intercept = NA, adj_r2 = NA, Trend = NA
          )
        }
      }) %>%
      mutate(scientific_name = sp) %>%
      select(scientific_name, BCR, everything())
    
    BCR_master[[sp]] <- BCR_trends
    
    write.csv(
      BCR_trends,
      paste0(output_dir, "/", gsub(" ", "_", sp), "_BCR_trends.csv"),
      row.names = FALSE
    )
    
    
    # ---------------------------------------------------
    # 2b. SAFE BCR PLOT
    # ---------------------------------------------------
    gg <- ggplot(train_bcr, aes(x = Year, y = adj_species_abundance)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
      facet_wrap(~ BCR, scales = "free_y") +
      labs(
        x = "Year",
        y = "Adjusted Species Abundance",
        title = paste("Linear Regression by BCR for:", sp)
      )
    
    ggsave(filename = paste0(output_dir, "/", gsub(" ", "_", sp),
                             "_BCR_plot.png"),
           plot = gg, width = 12, height = 10, dpi = 150)
  }
  
  
  # ---------------------------------------------------
  # 3. MASTER TABLES
  # ---------------------------------------------------
  state_table <- bind_rows(state_results)
  BCR_master_table <- bind_rows(BCR_master)
  
  write.csv(state_table,
            paste0(output_dir, "/STATE_master_trends_with_CI.csv"),
            row.names = FALSE)
  
  write.csv(BCR_master_table,
            paste0(output_dir, "/BCR_master_trends_with_CI.csv"),
            row.names = FALSE)
  
  list(
    State_Trends = state_table,
    BCR_Trends = BCR_master_table
  )
}

#call the function
analyze_all_species(data=bbs_tx_sub)


#Final training dataset creation
#state
state <- read.csv("trends/texas/STATE_master_trends_with_CI.csv")
state <- state %>% select(scientific_name,StateNum,slope)

training_state <- merge(state_trends,state,by=c("scientific_name","StateNum"))

training_state <- training_state %>%
  select(
    scientific_name,StateNum,
    {
      yr_cols <- names(.)[str_detect(names(.), "_\\d{4}$")]
      yr_cols[order(as.numeric(str_extract(yr_cols, "\\d{4}")))]
    },slope
  )

write.csv(training_state,"training/training_TX.csv",row.names = FALSE)

bcr <- read.csv("trends/texas/BCR_master_trends_with_CI.csv")
bcr <- bcr %>% select(scientific_name,BCR,slope)
training_bcr <- merge(bcr_trends,bcr,by=c("scientific_name","BCR"))
training_bcr <- training_bcr %>%
  select(
    scientific_name,BCR,
    {
      yr_cols <- names(.)[str_detect(names(.), "_\\d{4}$")]
      yr_cols[order(as.numeric(str_extract(yr_cols, "\\d{4}")))]
    },slope
  )

write.csv(training_bcr,"training/training_BCR.csv",row.names = FALSE)
 
#Merge the two datasets and include a spatial scale variable for BCR or State
state <- read.csv("training/training_TX.csv")
state$spatial_scale <- "state"
state <- state %>% select(-StateNum)

bcr <- read.csv("training/training_bcr.csv")
bcr$spatial_scale <- "bcr"
bcr <- bcr %>% select(-BCR)

training <- rbind(state,bcr)

write.csv(training, "training/training_texas.csv")
