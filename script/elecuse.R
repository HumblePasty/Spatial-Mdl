# this script means to attach new column to the main data
# the new column is the population density of the electricity consumption
# data source: https://figshare.com/articles/dataset/Global_1_km_1_km_gridded_revised_real_gross_domestic_product_and_electricity_consumption_during_1992-2019_based_on_calibrated_nighttime_light_data/17004523/1

# libs
library(raster)
# library(rgdal)
library(sp)
library(dplyr)

# load the data
elecuse = raster("./data/elecuse/EC2019.tif")

# load the main data
data = read.csv("./data/GaN2023_full.csv")

# get the locations in the main data
locs = data %>% select(Longitude, Latitude)

# project the locations to the land value data
# set the projection system of the locations to WGS84
spatial_loc = SpatialPoints(locs, proj4string = CRS("+proj=longlat +datum=WGS84"))
spatial_loc = spTransform(spatial_loc, CRS(proj4string(elecuse)))

# sample the land value
elecuselist = extract(elecuse, spatial_loc)

# merge the column to the main data
data$elec_use = elecuselist

# export the merged data
write.csv(data, "./data/GaN2023_full.csv")
