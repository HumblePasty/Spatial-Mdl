# this script means to attach new column to the main data
# the new column is the population density of the land values
# data source: https://placeslab.org/fmv-usa/

# libs
library(raster)
# library(rgdal)
library(sp)
library(dplyr)

# load the data
landvalue = raster("./data/landvalue/places_fmv_all.tif")

# load the main data
data = read.csv("./data/GaN2023_full.csv")

# get the locations in the main data
locs = data %>% select(Longitude, Latitude)

# project the locations to the land value data
# set the projection system of the locations to WGS84
spatial_loc = SpatialPoints(locs, proj4string = CRS("+proj=longlat +datum=WGS84"))
spatial_loc = spTransform(spatial_loc, CRS(proj4string(landvalue)))

# sample the land value
landvaluelist = extract(landvalue, spatial_loc)

# merge the column to the main data
data$landvalue_dollar_ha = landvaluelist

# export the merged data
write.csv(data, "./data/GaN2023_full.csv")
