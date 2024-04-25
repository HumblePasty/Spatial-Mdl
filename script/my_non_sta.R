library(BMNUS)
library(RepeatedMeasurements)
library(PairedMeasurements)
library(MappingUtilities)
library(BasicCodaFunctions)
library(sp)
library(ggplot2)
library(Matrix)
library(rstan)
library(dplyr)
library(tidyr)
library(tibble)

# load our own dataset
# work dir
setwd("./data")

# load the data
data = read.csv('GaN2023_full.csv')

# filter the data with "SQMReading" between 16 and 22 and within north America
f_data = data %>%
  select(ID, Latitude, Longitude, SQMReading) %>%
  filter(SQMReading >= 16 & SQMReading <= 22) %>%
  filter(Longitude >= -125 & Longitude <= -65 & Latitude >= 15) %>%
  drop_na()

rm('data')

measurements1 = f_data$SQMReading

# estimate the standard deviation
oEvaluation = RM_Evaluation(measurements1, 'SQM Reading')

plot(oEvaluation)

summary(oEvaluation)

rm("measurements1", "oEvaluation")

# sd = 1.37114
Add_Geography <- function(boundary_color = "gray70") {
  # add boundaries to the map
  map_df <- ggplot2::map_data("state")
  
  return (ggplot2::geom_path(ggplot2::aes(x = long, y = lat, group = group),
                            data = map_df,
                            colour = boundary_color))
}

# organize my dataset
my_dataset = select(f_data, ID, Longitude, Latitude, SQMReading)
colnames(my_dataset) = c("labno", "long", "lat", "conc")
# convert to tbl_df class
my_dataset = as_tibble(my_dataset)

ggplot() +
  Add_Geography(boundary_color = "gray50") +
  # Add_Path(Example_domain, "blue") +
  # Add_Path(Example_roi, "red") +
  Add_Points(my_dataset) +
  Refine_Map("lambert", 25, 40, latLimits = c(20, 50))


me_sd = 1.37
Value <- scaledLogit(my_dataset$conc, constSumValue = 1e6)
spatialData <-
  sp::SpatialPointsDataFrame(data.frame(long = my_dataset$long,
                                        lat = my_dataset$lat,
                                        row.names = my_dataset$labno),
                             data.frame(Value = Value,
                                        IndValue = "no",
                                        me_sd = rep.int(me_sd, length(Value)),
                                        row.names = my_dataset$labno),
                             proj4string=sp::CRS(Example_CRS_arg_longlat),
                             match.ID=TRUE)

# construct the basis function
spacings <- c("020", "025", "030", "035", "040",
              "045", "050", "055", "060",
              "070", "080", "090", "100")

spacings = "100"
# The domain: two columns, the first is the longitude and the second is the latitude.
my_domain = cbind(c(-125, -65, -65, -125), c(15, 15, 50, 50))
colnames(my_domain) = c("long", "lat")
# convert to spec_tbl_df
my_domain = as.data.frame(my_domain)

my_roi = cbind(c(-125, -65, -65, -125), c(15, 15, 50, 50))
colnames(my_roi) = c("long", "lat")
my_roi = as.data.frame(my_roi)

oBasisFunctions <- BasisFunctions(spacings,
                                  my_domain,
                                  spatialData,
                                  Example_CRS_arg_utm,
                                  seed = 777)
summary(oBasisFunctions)

# plot locations of the basis function centers
plot(oBasisFunctions, 100) +
  Add_Geography(boundary_color = "gray50") +
  # add a red rectangle to the map
  # geom_polygon(data = as.data.frame(my_domain), aes(x = long, y = lat), fill = "red") +
  Add_Path(my_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(20, 50))

# cross validation
filename <- paste0(normalizePath(path.package("BMNUS")), "\\stan\\BMNUS.stan")

# Translate Stan program BMNUS into C++ functions and then compile
# these functions
tr <- stanc(file = filename, model_name = "BMNUS")
BMNUS_sm <- stan_model(stanc_ret = tr, verbose=FALSE)

nFolds <- 10
nCpuCores <- 3

# Construct the class for the cross validation test.
# These computations require a long time.
oCrossValidation1 <- CrossValidation(spatialData,
                                     BMNUS_sm,
                                     oBasisFunctions,
                                     nFolds, nCpuCores)

# Save the results because they compuation require so long.
save(oCrossValidation1, file = paste0("CrossValidation_10-Folds\\Object.dat"))

# Plot the results
plot(oCrossValidation1)

# Print a summary of the cross validation test
summary(oCrossValidation1)

# Repeat for 20 folds
nFolds <- 20
oCrossValidation2 <- CrossValidation(spatialData,
                                     BMNUS_sm,
                                     oBasisFunctions,
                                     nFolds, nCpuCores,
                                     selected_spacings = c("020", "025", "030", "035", "040"))

save(oCrossValidation2, file = paste0("CrossValidation_20-Folds\\Object.dat"))
plot(oCrossValidation2)
summary(oCrossValidation2)

rm("filename", "tr", "nFolds", "nCpuCores")

load("BasisFunctions\\Spacing_030.dat")

csrBasisFuncs <- as(basis_functions$basisFuncs, "RsparseMatrix")

# Calculate quantities that are needed for the sparse CAR model;
# these quantities are derived from the adjacency matrix
carQuantities <- SparseCarQuantities(basis_functions)

stanData <- list( N = nrow(spatialData@data),
                  M = ncol(csrBasisFuncs),
                  UU = csrBasisFuncs@p + 1,
                  VV = csrBasisFuncs@j + 1,
                  WW = csrBasisFuncs@x,
                  nNeighborPairs = carQuantities$nNeighborPairs,
                  W_sparse = carQuantities$W_sparse,
                  D_sparse = carQuantities$D_sparse,
                  lambda = carQuantities$lambda,
                  X = spatialData@data$Value,
                  me_sd = spatialData@data$me_sd,
                  areLeftCensored = as.integer(spatialData@data$IndValue == "left"),
                  betaPar = c(2.5, 1.2),
                  gammaPar = c(2, 0.3),
                  cauchyPar = 3)