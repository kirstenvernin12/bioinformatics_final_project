#North American Breeding Bird Survey Data Exploration (2024 Release)
setwd("C:/Users/kirst/Desktop/UTA/Fall 2025/Bioinformatics/Final Project/bioinformatics_final_project")

#used to create training and testing datasets

#load libraries
library(dplyr)
library(data.table)
library(lubridate)
library(sf)

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
  select(-CountryNum,-StateNum, -Route, -RPID, -Year,-RecordedCar) 

#Get total cars observed
dist$TotalCarObs <- rowSums(dist[paste0("Car",1:50)],na.rm=TRUE)

#Remove Car and Noise stop data 1-50
dist <- dist[ , !names(dist) %in% c(paste0("Car", 1:50), paste0("Noise", 1:50))]

#Remove CountryNum, StateNum, Route, RPID, Year
weather <- weather %>% 
  select(RouteDataID,Month,Day,StartTime,EndTime) 

weather <- weather %>%
  select(RouteDataID,Month,Day,StartTime,EndTime)

#No filtering required because "AOU" only shared field
bbs_tx_bcr <- merge(bbs_tx_bcr,species, by="AOU")
bbs_tx_bcr <- merge(dist,bbs_tx_bcr, by="RouteDataID")
bbs_tx_bcr <- merge(weather,bbs_tx_bcr,by="RouteDataID")

#Create summary fields/tables and then clean up species (what to do with hybrids, ect. prior to plotting)
#Step 1 - create a combined Genus species field
bbs_tx_bcr <- bbs_tx_bcr %>% mutate(scientific_name = paste(Genus,Species, sep = " "))

bbs_tx_bcr <- bbs_tx_bcr %>% select(c(RouteDataID,Route,StateNum,Latitude,Longitude,
                                      BCR,Month,Day,Year,TotalCarObs,SpeciesTotal,
                                      scientific_name))%>%
  rename(Abundance=SpeciesTotal)


#organize fields
bbs_tx_bcr <- bbs_tx_bcr[,c("RouteDataID", "Route", "StateNum","Latitude", "Longitude", "BCR",
                          "Month", "Day", "Year", "Abundance", "TotalCarObs",
                          "scientific_name")] 
write.csv(bbs_tx_bcr,"raw_data/bbs_tx_bcr.csv",row.names = FALSE)

#filter for period of interest (1997-2023)
bbs_tx_sub <- bbs_tx_bcr %>% filter(Year>=1997)

#Import, filter, and process NOAA weather data
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


#Match weather stations with closest route point
weatherstations <- read.csv("raw_data/NOAA_weather/ghcnd-stations.csv")

#Get the state values only - these will be the ones added to the final state abundance df
tx_stations <- weatherstations %>% filter(STATE=="TX")
bbs_texas <- bbs_tx_sub %>% filter(StateNum==83)
#Implement the function
path <- "raw_data/NOAA_WEATHER/"
years <- 1997:2023
stations <- tx_stations$ID

weather_state <- process_weather_years(path, years, stations)
weather_state <- weather_state %>% mutate(DATE_YYYYMMDD=ymd(DATE_YYYYMMDD),
                                          Year=year(DATE_YYYYMMDD))

#Get annual PRCP, TMAX, TMIN by year for the state. This will be joined with the annual abundance data by species for the state.
annual_weather_tx <- weather_state %>% group_by(Year) %>% summarize(PRCP_mean=mean(PRCP,na.rm=TRUE), 
                                                                    TMAX_mean=mean(TMAX, na.rm=TRUE),
                                                                    TMIN_mean=mean(TMIN,na.rm=TRUE))




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

#Use the weatherstations list to filter the full NOAA data
#Implement the function
path <- "raw_data/NOAA_WEATHER/"
years <- 1997:2023
stations <- weatherstations

weather_bcr <- process_weather_years(path, years, stations)
#Get annual avg PRCP, TMAX, TMIN by year for the BCRs. This will be joined with the annual abundance data by species for the BCRs.
weather_bcr <- weather_bcr %>% mutate(DATE_YYYYMMDD=ymd(DATE_YYYYMMDD),
                                          Year=year(DATE_YYYYMMDD))

points <- as.data.frame(points)
stations_bcr <- points %>% select(BCR,NearestStation) %>% distinct() %>% rename(ID=NearestStation)
weather_bcr <- merge(weather_bcr,stations_bcr,by="ID")

annual_weather_bcrs <- weather_bcr%>% group_by(BCR,Year) %>% summarize(PRCP_mean=mean(PRCP,na.rm=TRUE), 
                                                                    TMAX_mean=mean(TMAX, na.rm=TRUE),
                                                                    TMIN_mean=mean(TMIN,na.rm=TRUE))

#PRCP: mm
#TMIN and TMAX: degrees celsius

#Calculate annual summary stats
#Filter out hybrids and unknown species
bbs_tx_sub <- bbs_tx_sub %>% filter(                                            
  !grepl(" / ", scientific_name),
  !grepl(" sp\\.", scientific_name),
  !grepl(" x ", scientific_name, ignore.case = TRUE)
)

#Abundance
#State
state_abund <- bbs_tx_sub %>% filter(StateNum==83) %>%
  group_by(Year,scientific_name) %>%
  summarise(Abundance_mean=mean(Abundance)) 
state_abund$region_id <- "TX"
state_abund$spatial_scale <- "State"


#BCR
BCR_abund <- bbs_tx_sub %>% 
  group_by(Year,scientific_name,BCR) %>%
  summarise(Abundance_mean=mean(Abundance)) 

BCR_abund <- BCR_abund %>%
  rename(region_id=BCR)

BCR_abund$spatial_scale <- "BCR"

#Cars observed (get average for the state and for the BCRs)
#State
env_state <- bbs_tx_sub %>% filter(StateNum==83) %>% select(Year,TotalCarObs)%>%
  group_by(Year)%>%
  summarize(CarsObs_mean=mean(TotalCarObs))
  
#BCR
env_bcr<- bbs_tx_sub %>% select(Year,BCR,TotalCarObs)%>%
  group_by(Year,BCR)%>%
  summarise(CarsObs_mean=mean(TotalCarObs)) %>%
  rename(region_id=BCR)

#state trends
state_trends <- merge(state_abund,annual_weather_tx,by="Year")
state_trends <- merge(state_trends,env_state,by="Year")

#bcr trends
annual_weather_bcrs <- annual_weather_bcrs %>% rename(region_id=BCR)
bcr_trends <- merge(BCR_abund,annual_weather_bcrs,by=c("region_id","Year"))
bcr_trends <- merge(bcr_trends,env_bcr,by=c("region_id","Year"))

training_texas <- rbind(state_trends,bcr_trends)

write.csv(training_texas, "training/training_texas.csv",row.names = FALSE)
