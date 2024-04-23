# this script means to attach new column to the main data
# the new column is the population density of the locations
# data source: https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11/data-download

# libs
library(raster)
# library(rgdal)
library(sp)
library(dplyr)

# load the data
popden  = raster("./data/popden/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")

# load the main data
data = read.csv("./data/GaN2023.csv")

# get the locations in the main data
locs = data %>% select(Longitude, Latitude)

# sample the population density
popdenlist = extract(popden, locs)

# merge the column to the main data
data$popden_km2 = popdenlist

# export the merged data
write.csv(data, "./data/GaN2023_full.csv")

# show a simple plot over the United States
# adjust the extent to the region of interest
# plot(popden, ext = extent(-125, -65, 25, 50))
# points(locs$Longitude, locs$Latitude, col = "red", pch = 20)
