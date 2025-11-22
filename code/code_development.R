# Extra code/code used prior to developing functions for automation

#DATA PROCESSING

#subset for 1997-2023 and species of interest
training <- bbs_tx_bcr %>%
  filter(scientific_name=="Chaetura pelagica"|scientific_name=="Chordeiles minor"|scientific_name=="Tachycineta bicolor"|
           scientific_name== "Petrochelidon pyrrhonota"|scientific_name=="Hirundo rustica"|
           scientific_name=="Dendrocygna autumnalis"|scientific_name=="Coragyps atratus"|
           scientific_name=="Ictinia mississippiensis"|scientific_name=="Athene cunicularia") %>%
  filter(Year>=1997)

write.csv(training, "raw_data/training_unfiltered.csv", row.names = FALSE)
#Black vulture
#State Trend - Black vulture
train_state <- training %>% filter(StateNum==83 & scientific_name== "Coragyps atratus") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
vul_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract slope and adjusted R²
slope <- coef(vul_lm)["Year"]
adj_r2 <- summary(vul_lm)$adj.r.squared

# Create formatted text for title
title_text <- paste0(
  "Linear Regression of Adjusted Species Abundance Over Time\n",
  "Slope = ", round(slope, 4),
  "   |   Adjusted R² = ", round(adj_r2, 4)
)

# Plot with stats in the title
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundancen (#Individuals/#Routes)",
     main = title_text)

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Coragyps atratus"
State_Trend <- slope

State_vul <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Black vulture
train_bcr <- training %>% filter(scientific_name== "Coragyps atratus") %>%
  group_by(Year,BCR) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

BCR_vul<- train_bcr %>%
  group_by(BCR) %>%
  group_modify(~ {
    m <- lm(adj_species_abundance ~ Year, data = .x)
    tibble(
      slope = coef(m)["Year"],
      intercept = coef(m)["(Intercept)"],
      summary(m)$adj.r.squared,
    )
  })


write.csv(BCR_vul, "BCR_Trends/black_vulture_BCR.csv")

ggplot(train_bcr, aes(x = Year, y = adj_species_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ BCR, scales = "free_y") +
  labs(
    x = "Year",
    y = "Adjusted Species Abundance",
    title = "Linear Regression of Adjusted Species Abundance Over Time by BCR"
  )

#Common raven
#State Trend - Common raven
train_state <- training %>% filter(StateNum==83 & scientific_name== "Corvus corax") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
bank_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(bank_lm)["Year"] #0.009326113  (Increasing for the state)
summary(bank_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Corvus corax"
State_Trend <- "Increasing"

State_raven <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Common raven
train_bcr <- training %>% filter(scientific_name== "Corvus corax") %>%
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

BCR_trend <- c("Increasing", "Decreasing", "Increasing", "Increasing","Stable",NA) #One with no trend
BCR <- c(18,19,20,21,35,36)
scientific_name <- c("Corvus corax","Corvus corax","Corvus corax","Corvus corax",
                     "Corvus corax","Corvus corax")

BCR_raven <- data.frame(BCR, BCR_trend,scientific_name)


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

#Tree Swallow
#State Trend - Tree Swallow
train_state <- training %>% filter(StateNum==83 & scientific_name== "Tachycineta bicolor") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
bank_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(bank_lm)["Year"] #0.007909698  (Stable)
summary(bank_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Tachycineta bicolor"
State_Trend <- "Stable"

State_tree <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Tree Swallow
train_bcr <- training %>% filter(scientific_name== "Tachycineta bicolor") %>%
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

BCR_trend <- c("Decreasing", "Increasing", NA, "Decreasing","Increasing",NA) #Two with no trend
BCR <- c(18,19,20,21,25,37)
scientific_name <- c("Tachycineta bicolor","Tachycineta bicolor","Tachycineta bicolor",
                     "Tachycineta bicolor","Tachycineta bicolor","Tachycineta bicolor")

BCR_tree <- data.frame(BCR, BCR_trend,scientific_name)

#Cliff Swallow
#State Trend - Cliff Swallow
train_state <- training %>% filter(StateNum==83 & scientific_name== "Petrochelidon pyrrhonota") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
bank_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(bank_lm)["Year"] #-0.5290313   (Decreasing)
summary(bank_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Petrochelidon pyrrhonotar"
State_Trend <- "Decreasing"

State_cliff <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Cliff Swallow
train_bcr <- training %>% filter(scientific_name== "Petrochelidon pyrrhonota") %>%
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

BCR_trend <- c("Decreasing", "Decreasing", "Decreasing", "Decreasing","Increasing","Decreasing","Increasing","Decreasing") 
BCR <- c(18,19,20,21,25,35,36,37)
scientific_name <- c("Petrochelidon pyrrhonotar","Petrochelidon pyrrhonota","Petrochelidon pyrrhonota",
                     "Petrochelidon pyrrhonota","Petrochelidon pyrrhonota","Petrochelidon pyrrhonota","Petrochelidon pyrrhonota",
                     "Petrochelidon pyrrhonota")

BCR_cliff <- data.frame(BCR, BCR_trend,scientific_name)

#Barn Swallow
#State Trend - Barn Swallow
train_state <- training %>% filter(StateNum==83 & scientific_name== "Hirundo rustica") %>%
  group_by(Year) %>%
  summarise(adj_species_abundance=round(sum(SpeciesTotal)/n_distinct(Route)),.groups= "drop")

# Fit linear regression model
bank_lm <- lm(adj_species_abundance ~ Year, data = train_state)

# Extract the slope (coefficient for Year)
coef(bank_lm)["Year"] #-0.3456679    (Decreasing)
summary(bank_lm)
# Plot data and regression line
plot(train_state$Year, train_state$adj_species_abundance,
     xlab = "Year",
     ylab = "Adjusted Species Abundance",
     main = "Linear Regression of Adjusted Species Abundance Over Time")

abline(bank_lm, lwd = 2)

StateNum <- 83
scientific_name <- "Hirundo rustica"
State_Trend <- "Decreasing"

State_barn <- data.frame(StateNum,State_Trend,scientific_name)

#BCR Trend for TX - Barn Swallow
train_bcr <- training %>% filter(scientific_name== "Hirundo rustica") %>%
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

BCR_trend <- c("Increasing", "Decreasing", "Stable", "Decreasing","Decreasing","Decreasing","Decreasing","Decreasing") 
BCR <- c(18,19,20,21,25,35,36,37)
scientific_name <- c("Hirundo rustica","Hirundo rustica", "Hirundo rustica", "Hirundo rustica", "Hirundo rustica", 
                     "Hirundo rustica", "Hirundo rustica", "Hirundo rustica")

BCR_barn <- data.frame(BCR, BCR_trend,scientific_name)

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
State <- merge(State,State_tree, all=TRUE)
State <- merge(State,State_barn, all=TRUE)
State <- merge(State,State_cliff, all=TRUE)
State <- merge(State,State_raven,all=TRUE)
State <- merge(State,State_vul,all=TRUE)

BCR <- merge(BCR_bank,BCR_chimney, all=TRUE)
BCR <- merge(BCR,BCR_night,all=TRUE)
BCR <- merge(BCR,BCR_tree,all=TRUE)
BCR <- merge(BCR,BCR_barn,all=TRUE)
BCR <- merge(BCR,BCR_cliff,all=TRUE)
BCR <- merge(BCR,BCR_raven,all=TRUE)
BCR <- merge(BCR,BCR_vul,all=TRUE)

training$StateNum <- as.numeric(training$StateNum)

training <- training %>% left_join(State, by=c("StateNum","scientific_name"))
training <- training %>% left_join(BCR,by=c("BCR","scientific_name"))

#create csv of the bbs data
write.csv(training,"raw_data/bbs_texas_training.csv", row.names=FALSE)

training <- read.csv("raw_data/bbs_texas_training.csv")


#cars detected is not informative, removed
#env_route_year <- env_route_year %>% select(-cars_detected)

#Create NoiseDetected field (1=Y,0=N)
dist$NoiseDetected <- ifelse(
  rowSums(dist[paste0("Noise", 1:50)] == 1, na.rm = TRUE) > 0,
  1,
  0
)

# Annual summaries per route × year
env_route_year <- env %>%
  group_by(Route, StateNum, BCR, Year) %>%
  summarise(
    Start_temp_mean=mean(StartTemp,na.rm=TRUE),
    End_temp_mean=mean(EndTemp,na.rm=TRUE),
    mean_cars_obs=mean(TotalCarObs,na.rm=TRUE),
    PRCP_sum = sum(PRCP, na.rm = TRUE),
    TMAX_mean = mean(TMAX, na.rm = TRUE),
    TMIN_mean = mean(TMIN, na.rm = TRUE),
    noise_detected = max(NoiseDetected, na.rm = TRUE),  # 1 if any detection
    .groups = "drop"
  )

# Aggregate per-route predictors to State-level
state_env <- env_route_year %>%
  group_by(StateNum, Year) %>%
  summarise(
    Start_temp_stateavg=mean(Start_temp_mean,na.rm=TRUE),
    End_temp_stateavg=mean(End_temp_mean,na.rm=TRUE),
    mean_cars_obs=mean(mean_cars_obs,na.rm=TRUE),
    PRCP_state = mean(PRCP_sum, na.rm = TRUE),
    TMAX_state = mean(TMAX_mean, na.rm = TRUE),
    TMIN_state = mean(TMIN_mean, na.rm = TRUE),
    noise_detected = max(noise_detected, na.rm = TRUE),  # 1 if any detection
    .groups = "drop"
  )

# Aggregate per-route predictors to BCR-level
bcr_env <- env_route_year %>%
  group_by(BCR, Year) %>%
  summarise(
    Start_temp_bcravg=mean(Start_temp_mean,na.rm=TRUE),
    End_temp_bcreavg=mean(End_temp_mean,na.rm=TRUE),
    mean_cars_obs=mean(mean_cars_obs,na.rm=TRUE),
    PRCP_bcr = mean(PRCP_sum, na.rm = TRUE),
    TMAX_bcr = mean(TMAX_mean, na.rm = TRUE),
    TMIN_bcr = mean(TMIN_mean, na.rm = TRUE),
    noise_detected = max(noise_detected, na.rm = TRUE),  # 1 if any detection
    .groups = "drop"
  )

training_state_long <- left_join(env, state_env, by=c("StateNum","Year"))
training_bcr_long <- left_join(env, bcr_env,by=c("BCR","Year"))
abundance <- training_state_long %>% filter(StateNum==83) %>%
  group_by(Year,scientific_name)%>%
  summarise(adj_avg_abundance = mean(sum(SpeciesTotal) /
                                       n_distinct(Route)))
training_state_long <- left_join(training_state_long,abundance,by=c("Year", "scientific_name"))

training_state_long <- training_state_long %>% filter(StateNum==83)


#remove unneeded fields and convert to wide format
training_state_long <- training_state_long %>% select(c(-RouteDataID,-RouteTypeID,-RouteTypeDetailID,-Active,-AOU,-RPID,
                                                        -Month,-Day,-StartTime,-EndTime,-StartTemp,-EndTemp,-Count10,-Count20,-Count30,
                                                        -Count40,-Count50,-StopTotal,-SpeciesTotal,-RecordedCar,-TotalCarObs,-NoiseDetected,
                                                        -Seq, -French_Common_Name, -Genus,-Species,-NearestStation,-DATE_YYYYMMDD,-PRCP,-TMAX,-TMIN,-geometry))

#Need to figure out how to join and get into the wide format!!!!

test<- training_state_long %>% filter(StateNum==83) %>%
  pivot_longer(cols = c(Start_temp_stateavg,End_temp_stateavg, mean_cars_obs,PRCP_state, TMAX_state, TMIN_state,
                        noise_detected,adj_avg_abundance),
               names_to = "var", values_to = "value") %>%
  mutate(feature = paste0(var, "_", Year)) %>%
  select(Route,RouteName,CountryNum,StateNum, Stratum,BCR,Year,English_Common_Name,Order,Family,scientific_name,
         feature, value) %>%
  pivot_wider(names_from = feature, values_from = value)





test <- training_texas%>% 
  pivot_longer(cols=c(adj_avg_abundance,Start_temp_stateavg,End_temp_stateavg,mean_cars_obs,PRCP_state,
                      TMAX_state,TMIN_state,noise_detected),names_to="var",values_to="value") %>%
  mutate(feature=paste0(var,"_",Year)) %>%
  select(Route,RouteName,CountryNum,StateNum,Stratum,BCR,Year,English_Common_Name,Order,Family, Genus, scientific_name,
         feature,value) %>%
  pivot_wider(names_from = feature,values_from=value, values_fill=0)

#remove unnecessary columns
#Data processing and splitting of state and BCR training dataset

#Remove unneeded fields. This will open memory and increase speed because large datasets with many unused columns take more
#memory and longer processing times. Fewer columns also make preprocessing and feature selection steps easier to follow. 
#This also reduces the risk of accidental leakage. Columns not intended as predictors may inadvertently leak information
#about the target (if they contain IDs, dates, or future information - basically, better safe than sorry!) It is best
#practice to keep only the columns you intend to use as predictors or target when feeding data into preprocessing/recipe 
# steps. 

#state training dataset
training_TX <- training %>% filter(StateNum==83) %>% 
  select(RouteName, Active, Stratum, BCR, Month, Day, Year, Avg_Temp, SpeciesTotal, 
         RecordedCar, TotalCarObs, NoiseDetected, English_Common_Name, Order, Family, 
         scientific_name, State_Trend, BCR_trend, PRCP, TMAX, TMIN) %>% rename(BCR_Trend=BCR_trend)


#BCR training dataset
training_BCR <- training %>% 
  select(RouteName, Active, Stratum, BCR, Month, Day, Year, Avg_Temp, SpeciesTotal, 
         RecordedCar, TotalCarObs, NoiseDetected, English_Common_Name, Order, Family, 
         scientific_name, State_Trend, BCR_trend, PRCP, TMAX, TMIN) %>% rename(BCR_Trend=BCR_trend)

#Remove the one BCR that is NA
training_BCR <- training_BCR  %>% filter(!is.na(BCR_Trend))

#write training datasets to training folder
write.csv(training_TX,"training/training_TX.csv",row.names=FALSE)
write.csv(training_BCR,"training/training_BCR.csv",row.names=FALSE)
