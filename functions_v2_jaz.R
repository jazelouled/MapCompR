# Read data ---------------------------------------------------------------------------
library(raster)
library(data.table)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
aa

# Here this is still manual
MEDITS_data <- read.csv("~/Dropbox/2023_MapComparisonPackage_FishMIP/00input/MEDITS_data/Abundancia_biomasa_GSA1_2_6.csv")
MEDITS_data_ <- MEDITS_data[MEDITS_data$ESPECIE %like% "Engraulis encrasicolus",]

MEDITS_data_$longitud <- rowMeans(MEDITS_data_[, c("LONGITUD_INI", "LONGITUD_VIR")])
MEDITS_data_$latitud <- rowMeans(MEDITS_data_[, c("X.LATITUD_INI", "LATITUD_VIR")])

# MEDITS_data_ <- MEDITS_data[MEDITS_data$ESPECIE %like% "Sardina pilchardus",]
Ecospace_data <- dir("~/Dropbox/2023_MapComparisonPackage_FishMIP/00input/asc/", pattern = "Adult anchovy", full.names = T)
# Ecospace_data <- dir("~/Dropbox/2023_MapComparisonPackage_FishMIP/00input/asc/", pattern = "Adult sardine", full.names = T)
Ecospace_data <- Ecospace_data[Ecospace_data %like% "Biomass"]

#############
# obsInPred #
#############

# Keep only observed data within the extent of the predicted value by the  MEM --------
# Rationale: this function subsets the observed locations to keep only the ones that are within the predicted
# extent by the MEM. It uses a spatial mask of the MEM, and deletes the points that are out of it.

# The function has three arguments: 
# observed: object containing the observed data (MEDITS in this case)
# observedLon: name of the column with longitude data
# observedLat: name of the column with latitude data
# predicted: raster with the extent we want to work in. It is used to create a mask and delete
#            the observed points out of the raster
observed <- MEDITS_data_

predicted <- Ecospace_data

# Keep observations within predicted extent

# Version 1
obsInPred <- function(observed, observedLon, observedLat, predicted){
  predicted_data_mask <- raster(predicted[1])
  predicted_mask <- predicted_data_mask/predicted_data_mask
  polygonDeletePoints <- rasterToPolygons(predicted_mask, dissolve = T)
  polygonDeletePoints_ <- st_as_sf(polygonDeletePoints)
  observedLon <- observedLon
  observedLat <- observedLat
  coordinatesAndAttributes <- sf::st_as_sf(observed, coords = c(observedLon, observedLat))
  xysf <- st_as_sf(as.data.frame(observed), coords = c(observedLon, observedLat), crs = 4326)
  xy_intersect <- st_intersection(polygonDeletePoints_, xysf)
  xy_intersect_ <- tidyr::extract(xy_intersect, geometry, into = c('Lon', 'Lat'), '\\((.*),(.*)\\)', conv = T)
  
  return(xy_intersect_)
}

system.time(observed <- obsInPred(MEDITS_data_, observedLon = "longitud", observedLat = "latitud", Ecospace_data))


# Version 2
obsInPred <- function(observed, observedLon, observedLat, predicted){
  # Polygons 
  mask <- raster(predicted[1]) 
  polygon <- rasterToPolygons(mask, dissolve = F)
  sf_polygon <- st_as_sf(polygon)
  boundary <- st_union(sf_polygon)
  # Points 
  points <- SpatialPoints(observed[,c(observedLon,observedLat)])
  sf_points <- st_as_sf(points)
  
  # Intersect 
  st_crs(sf_points) <- st_crs(boundary)
  xy_intersect <- st_intersection(boundary, sf_points)
  coordinates <- lapply(xy_intersect, function(point) st_coordinates(point))
  xy_intersect_ <- do.call(rbind, coordinates)
  colnames(xy_intersect_) <- c(observedLon, observedLat)
  xy_intersect_ <- as.data.frame(xy_intersect_)
  
  # Joint the intersected 
  observed_ <- semi_join(observed, xy_intersect_, by = c(observedLon, observedLat))
  return(observed_)
}

observed <- obsInPred(observed = MEDITS_data_,
                      observedLon = "longitud",
                      observedLat = "latitud",
                      predicted = Ecospace_data)



####################
# obsPresBiomassDF #
####################

# Put data frame in the format we need for comparison ---------------------
# Rationale: this function will take the dataset obtained with the obsInPred function and extract the biomass
# values in the observed latitude and longitude from the predicted values by the MEM. It also adds the depth
# of the location by using a bathymetry layer.
# *how to not have the bathymetry in local so it can be downloaded. Maybe github?

# The function has three arguments: 
# observed: data frame obtained with the obsInPred function --> need to change this argument's name
# observedBiomass: name of the biomass column
# predicted: vector of paths to the MEM prediction we want to work with. 
# Afegir argument raster


bathymetry <-"~/Dropbox/2022_SouthHemisphere_SDM/bathymetry/"
bathy <- raster(paste0(bathymetry, "/GEBCO_2014_2D.nc"))  # bathymetry: This we have to check how to put it online so it is downloaded with the package



obsPresBiomassDF <- function(observed, observedLon, observedLat, observedBiomass, predicted, rasterVariable, bufferExtract){
  listExtract <- list()
  listExtract_ <- list()
  
  for (j in 1:length(predicted)){
    if (j==24) next
    predictedRaster <- predicted[j]
    fileName <- sub("\\..*", "", predictedRaster)
    fileName_ <- substr(fileName, nchar(fileName) - 4, nchar(fileName))
    predictedRaster_ <- raster(predictedRaster)
    numberMonth <- as.numeric(fileName_)
    originDate <- lubridate::ymd("2000-01-01")
    date <- lubridate::ymd(originDate+months(numberMonth-1))
    year <- year(date)
    observed_ <- observed[observed$Year %like% year,]

    
    
        for (i in 1:length(observed_$Year)){
          print(i)
          columnInfo <- observed_[[observedBiomass]]
          observedLon_ <- observed_[[observedLon]]
          observedLat_ <- observed_[[observedLat]]

          bufferExtract <- bufferExtract*1000
          
          extracted_env <- raster::extract(predictedRaster_, cbind(observedLon_[i], observedLat_[i]), buffer=0, fun=mean, na.rm=TRUE)  # extract data
          extracted_rasterVariable <- raster::extract(rasterVariable, cbind(observedLon_[i], observedLat_[i]), bufferExtract=bufferExtract, fun=mean, na.rm=TRUE)  # extract biomass data
          df_predObs <- data.frame("longitude" = numeric(), "latitude" = numeric(), "biomass_observed" = numeric(), "biomass_predicted" = numeric(), "year" = numeric(), "variable" = numeric())
          df_predObs[1,1] <- observedLon_[i]
          df_predObs[1,2] <- observedLat_[i]
          df_predObs[1,3] <- columnInfo[i]
          df_predObs[1,4] <- extracted_env
          df_predObs[1,5] <- year
          df_predObs[1,6] <- extracted_rasterVariable
          listExtract[[i]] <- df_predObs
    }
    
      aa <- do.call(rbind, listExtract)
      listExtract_[[j]] <- aa
      
  }
  
  DF_standardized <- do.call(rbind, listExtract_)
  
  DF_standardized_ <- DF_standardized[complete.cases(DF_standardized),]
  
  
  return(DF_standardized_)
  
}

hh <- obsPresBiomassDF(observed, observedLon = "longitud", observedLat = "latitud", observedBiomass = "BIOMASA", predicted = Ecospace_data, rasterVariable = bathy, bufferExtract = 20)



#####################
# visual_inspection #
#####################

# Visualize general trends in the data ------------------------------------
# Rationale: this function will take the dataset obtained by obsPresBiomassDF and provide some plots that show general trends. By
# now we have (poorly) implemented temporal trends and predicted vs. observed. 

# The function has five arguments: 
# df: data frame obtained with the obsPresBiomassDF function
# observed: column of observed data
# predicted: column of predicted data
# year: column of year data
# depth: column of depth data





# (plots have to be re-thought and improved)
visual_inspection <- function(df, observed, predicted, year, depth){
  
  observed <- pull(df, {{ observed }})
  predicted <- pull(df, {{ predicted }})
  year <- pull(df, {{ year }})
  depth <- pull(df, {{ depth }})
  
  
  plot_temporal <-  ggplot() +
                            geom_point(data=df, aes(year, predicted)) +
                            geom_smooth(data=df, aes(year, predicted))+
                            geom_point(data=df, aes(year, observed)) +
                            geom_smooth(data=df, aes(year, observed)
                            )
  
  
  plot_depths <-  ggplot() +
                           geom_point(data=df, aes(depth, predicted)) +
                           geom_point(data=df, aes(depth, observed)) +
                           geom_smooth(data=df, aes(depth, observed)) +
                           geom_smooth(data=df, aes(depth, predicted)
                           )
  
  
 plot_vs <- ggplot() +
                  geom_point(data=df, aes(observed, predicted*100)) +
                  geom_smooth(data=df, aes(observed, predicted*100))
  
  
  
  grid_arrange <- grid.arrange(plot_temporal, plot_depths, plot_vs, nrow = 3) # Arrange plots in a grid
  plots <- grid_arrange
  
  return(plots)
  
}

visual_inspection(hh, biomass_observed, biomass_predicted, year, depth) 




#####################
# calculate_metrics #
#####################

# Calculate metrics from Hipsey et al., 2020 ------------------------------
# This function will calculate metrics that allow us to compare the predicted and observed values. They 
# are extracted from Hipsey et al., 2020. Plots are provided too (bias & residuals for the moment).

# The function has three arguments: 
# df: data frame obtained with the obsPresBiomassDF function
# observed: column of observed data
# predicted: column of predicted data




calculate_metrics <- function(df, observed, predicted) {
  
  observed <- pull(df, {{ observed }})
  predicted <- pull(df, {{ predicted }})
  
  # Calculate Bias
  bias <- mean(predicted - observed)
  
  # Calculate RMSE
  rmse <- sqrt(mean((predicted - observed)^2))
  
  # Calculate MEF
  mef <- 1 - (rmse^2 / var(observed))
  
  # Calculate NMAE --> MAE/mean(observed) --> is this the good function?
  nmae <- mean(abs(predicted - observed)) / (max(observed) - min(observed))
  
  # Calculate MAE
  mae <- mean(abs(predicted - observed))
  
  # Create plots
  plot_bias <- ggplot(data.frame(x = observed, y = predicted), aes(x, y)) +
                      geom_point() +
                      geom_abline(intercept = bias, slope = 1, color = "red") +
                      labs(title = "Bias Plot", x = "Observed", y = "Predicted")
  
  plot_residuals <- ggplot(data.frame(residuals = predicted - observed), aes(residuals)) +
                           geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
                           labs(title = "Residuals Plot", x = "Residuals", y = "Frequency")
  
  # Arrange plots in a grid
  grid_arrange <- grid.arrange(plot_bias, plot_residuals, nrow = 1)
  
  # Create a dataframe with the metrics
  result <- data.frame(bias = bias, rmse = rmse, mef = mef, nmae = nmae, mae = mae)
  return(result)
  
  plots <- grid_arrange
  return(plots)
  
}


calculate_metrics(hh, biomass_observed, biomass_predicted)








