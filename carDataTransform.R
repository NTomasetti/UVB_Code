library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

# Data is available in three parts from https://data.transportation.gov/Automobiles/Next-Generation-Simulation-NGSIM-Vehicle-Trajector/8ect-6jqj


cars1 = read.table('vehicle-trajectory-data/trajectories1.txt')

cars2 = read.table('vehicle-trajectory-data/trajectories2.txt')

# IDs restart at one for each data segment, so add 10000/20000 to second / third set so that the IDs are unique for the combined data
cars2[,1] = cars2[,1] + 10000 
cars2[,15] = cars2[,15] + 10000

cars3 = read.table('vehicle-trajectory-data/trajectories3.txt')
cars3[,1] = cars3[,1] + 20000
cars3[,15] = cars3[,15] + 20000

cars = rbind(cars1, cars2, cars3)
rm(cars1, cars2, cars3)

# We are using x = distance from far left of road to far right (horizontal deviation)
# and y = distance travelled from bottom of road to top (ie how far along the car has moved)
# The UVB paper flips this in the plot so x = distance, y = deviatoin

colnames(cars) = c('ID', 'frame', 'totalFrames', 'time', 'x', 'y', 
                   'globalX', 'globalY', 'length', 'width', 'class',
                   'veloc', 'accel', 'lane', 'proceeding', 'following', 
                   'spacing', 'headway')

# Check for lane changes, entering, exiting (lane 6&7 are enter/exit), and get the starting lanes
cars %>%
  group_by(ID) %>%
  summarise(medLane = median(lane),
            changed = any(lane != medLane),
            enterExit = any(lane > 5),
            startLane = head(lane, 1)) %>%
  ungroup() %>%
  right_join(cars, by = "ID") -> cars

# Calcualte the distance travelled 

cars %>%
  filter(enterExit == FALSE) %>%
  group_by(ID) %>%
  mutate(time = time - min(time),
         n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         delta = atan2(y - ylag, x - xlag),
         dist = cumsum(v)) -> cars

# Define spline degree and extract 100 vehicles per lane that did not change lane to estimate spline models
degree <- 50

cars %>%
  filter(changed == FALSE) %>%
  group_by(lane) %>%
  filter(ID %in% head(unique(ID), 100)) %>%
  .$ID %>%
  unique() -> splinesID

# Input pos = (x, y) and centre, a dataframe of distance / spline xhat / spline yhat to transform coordinates
# The function will find the element of the xhat / yhat combo of centre which minimises the Euc. distance to the input (x, y).
transformCoord <- function(pos, centre){
  x <- pos[1]
  y <- pos[2]
  centre <- centre[abs(centre$yhat - y) < 15,] 
  centre$dist <- sqrt((x - centre$xhat)^2 + (y - centre$yhat)^2)
  closest <- centre[which.min(centre$dist),]
  relX <- sign(x - closest$xhat) * closest$dist
  relY <- closest$d
  return(c(relX, relY))
}

# Loop through each lane to:
# 1) Fit spline models
# 2) Evaluate spline on a fine grid of distances
# 3) Loop through each vehicle in the lane to calculate relative coordinates
# Result: Relcoord dataframe with relX / relY / time / ID
# relX is the lateral lane deviation
# relY is the distance travelled along the lane (note that lateral lane position is given by y in the UVB paper)
relCoord <- NULL
skip <- NULL
for(i in 1:5){
  cars %>%
    filter(lane == i  & ID %in% splinesID) -> carSubset
  modX <- smooth.spline(carSubset$dist, carSubset$x, df = degree)
  modY <- smooth.spline(carSubset$dist, carSubset$y, df = degree)
  
  cars %>%
    filter(startLane == i & !ID %in%splinesID) %>% 
    ungroup() -> toPredict
  
  centrePath <- data.frame(d = seq(min(toPredict$dist), max(toPredict$dist), length.out = 1000000))
  centrePath$xhat <- predict(modX, centrePath$d)$y
  centrePath$yhat <- predict(modY, centrePath$d)$y
  
  for(j in seq_along(unique(toPredict$ID))){
    carJ <- filter(toPredict, ID == unique(toPredict$ID)[j])
    rel <- select(carJ, x, y) %>%
      apply(1, function(row) transformCoord(c(row[1], row[2]), centrePath)) %>% 
      t() %>% as.data.frame()
    if(nrow(rel) == 1){
      skip <- c(skip, carJ$ID[1])
      next
    }
    
    colnames(rel) <- c('relX', 'relY')
    rel$ID <- carJ$ID
    rel$time <- carJ$time
    
    relCoord <- rbind(relCoord, rel)
    print(paste(i, j))
  }
}

# Combine datasets to add missing info back into data

cars %>% 
  filter(!ID %in% c(skip, splinesID)) %>%
  left_join(relCoord) %>%
  group_by(ID) %>%
  mutate(n = seq_along(time)) %>%
  rename(laneDev = relx) %>%
  filter(n > 2) %>%
  ungroup() %>% 
  select(ID, laneDev, changed, lane, startLane, x, y, time) -> cars

# Ensure we only select cars with more than 250  obs

cars %>%
  filter(changed == FALSE) %>%
  group_by(ID) %>%
  summarise(n = n()) %>%
  filter(n >= 250) %>%
  .$ID -> highObs

# Select 500 of these cars
subID <- sample(highObs, 500)

# Transform to matrix for VB 
cars %>%
  filter(ID %in% subID) %>%
  group_by(ID) %>%
  mutate(n = seq_along(ID)) %>%
  filter(n <= 250) %>%
  select(ID, n, laneDev) %>%
  spread(ID, laneDev) %>% 
  select(-n) %>%
  as.matrix() -> subMat

saveRDS(subMat, 'subCarsUVB.RDS')

# Pick the first car travelling in each lane

cars %>%
  filter(changed == FALSE) %>%
  group_by(lane) %>%
  filter(time == min(time)) %>%
  .$ID -> firstCars

# Plot the movement of these cars using real coordinates
# as laneDev = distance from vehicle to midline, we can get the midline as x - laneDev
# Note that we flip x / y from the system provided in the data for the plot
cars %>%
  filter(ID %in% firstCars) %>%
  mutate(mid = x - laneDev) %>%
  ggplot() + geom_path(aes(y, x, group = ID)) + 
  geom_path(aes(y, mid, group = ID), colour = 'red', linetype = 'dashed') + 
  theme_bw() + 
  labs(y = 'Lateral Lane Position (feet)', x = 'Distance Travelled (feet)')
