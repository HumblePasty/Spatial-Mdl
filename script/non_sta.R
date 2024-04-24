##
## This file lists the R scripts that are used in the User's Guide.
##

##
## Section "Preparatory Steps"
##

library(BMNUS)
library(RepeatedMeasurements)
library(PairedMeasurements)
library(MappingUtilities)
library(BasicCodaFunctions)
library(sp)
library(ggplot2)
library(Matrix)
library(rstan)

##
## Section "Estimate the Standard Deviation of the Measurement Error"
##

# Transform the titanium concentrations (which are repeated measurements
# of a reference sample) using the isometric log-ratio transform
measurements1 <- scaledLogit(Example_repeat_meas1$conc, constSumValue = 1.0e6)

# Construct the class to evaluate the repeated measurements
oEvaluation1 <- RM_Evaluation(measurements1, "Ilr-transformed concentrations")

# Plot the distribution of the repeated measurements and plot a
# quantile-quantile plot for which the theoretical distribution is
# normal
plot(oEvaluation1)

# Print summary statistics of the repeated measurements
summary(oEvaluation1)

measurements2 <- scaledLogit(Example_repeat_meas2$conc, constSumValue = 1.0e6)
oEvaluation2 <- RM_Evaluation(measurements2, "Ilr-transformed concentrations")
plot(oEvaluation2)
summary(oEvaluation2)

# Selected value for the standard deviation of the measurement error
me_sd <- 0.032

rm("measurements1", "oEvaluation1", "measurements2", "oEvaluation2")

##
## Section "Define the Region of Interest and the Domain"
##

# Add_Geography
#
# Add pertinent geography to the map for the statistical modeling.
#
# This function has been developed specifically for the data in
# package BMNUS, for which the processing is descripted in the User's Guide.
# This function must be modified to make it suitable for other areas.
#
# The function value is a list with two ggplot objects. One object pertains to
# the state boundaries, and the other to geographic names.
#
Add_Geography <- function(boundary_color = "gray70") {
  
  map_df <- ggplot2::map_data("state", region = c("Virginia", "North Carolina",
                                                  "South Carolina", "Georgia",
                                                  "Alabama", "Florida"))
  
  text_df <- data.frame(x = c(-86.7, -84, -82, -81, -80, -81.7, -78, - 86),
                        y = c(34.2, 34.2, 34.8, 35.9, 37, 28.5, 31, 28.5),
                        label = c("Alabama", "Georgia",
                                  "South\nCarolina", "North\nCarolina",
                                  "Virginia", "Florida", "Atlantic Ocean",
                                  "Gulf of Mexico"))
  return(
    list(ggplot2::geom_path(ggplot2::aes(x = long, y = lat, group = group),
                            data = map_df,
                            colour = boundary_color),
         ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = label),
                            size = 2.5,
                            data = text_df)))
}

# Plot the domain and the region of interest
ggplot() +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_domain, "blue") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

##
## Section "Organize the Spatial Data"
##
Example_dataset = Example_dataset
head(Example_dataset)

# Plot the locations of the measurements, as well as the domain and
# the region of interest
ggplot() +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_domain, "blue") +
  Add_Path(Example_roi, "red") +
  Add_Points(Example_dataset) +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

# Transform the titanium concentrations using the isometric log-ratio
# transform
Value <- scaledLogit(Example_dataset$conc, constSumValue = 1e6)

# Store the transformed concentrations, the associated censor indicators,
# and the associated standard deviations of the measurement errors in a
# spatial points dataframe.
spatialData <-
  sp::SpatialPointsDataFrame(data.frame(long = Example_dataset$long,
                                        lat = Example_dataset$lat,
                                        row.names = Example_dataset$labno),
                             data.frame(Value = Value,
                                        IndValue = "no",
                                        me_sd = rep.int(me_sd, length(Value)),
                                        row.names = Example_dataset$labno),
                             proj4string=sp::CRS(Example_CRS_arg_longlat),
                             match.ID=TRUE)

rm("me_sd", "Value")

##
## Section "Generate the Basis Functions"
##

# Maximum distance (in kilometers) between the transect and a field sample
maxOffset <- 10

# Location of end points (specified with latitude and longitude)
# of the first transect
end_points1 <- data.frame(long = c(-83.48355, -81.21104),
                          lat = c(33.30348, 29.69067))

# Construct a class to view the transect
oViewTransect1 <- ViewTransect(spatialData, maxOffset, end_points1,
                               Example_CRS_arg_utm)

# Plot the transect and the field samples associated with it
plotSampleLocations(oViewTransect1) +
  Add_Geography() +
  Add_Path(Example_domain, color = "blue") +
  Refine_Map("lambert", 25, 40,
             latLimits = c(29, 34),
             longLimits = c(-85.5, -80))

# Plot the values of the field samples against distance along the transect
plot(oViewTransect1, "Ilr-transformed\nconcentration", span = 0.60, se = FALSE)

# Second transect
end_points2 <- data.frame(long = c(-86, -87),
                          lat = c(33, 31.3))

oViewTransect2 <- ViewTransect(spatialData, maxOffset, end_points2,
                               Example_CRS_arg_utm)

plotSampleLocations(oViewTransect2) +
  Add_Geography() +
  Add_Path(Example_domain, color = "blue") +
  Refine_Map("lambert", 25, 40,
             latLimits = c(30,34),
             longLimits = c(-89,-84))

plot(oViewTransect2, "Ilr-transformed\nconcentration", span = 0.60, se = FALSE)

# Third transect
end_points3 <- data.frame(long = c(-82.3, -80),
                          lat = c(34, 32.5))

oViewTransect3 <- ViewTransect(spatialData, maxOffset, end_points3,
                               Example_CRS_arg_utm)

plotSampleLocations(oViewTransect3) +
  Add_Geography() +
  Add_Path(Example_domain, color = "blue") +
  Refine_Map("lambert", 25, 40,
             latLimits = c(31, 35),
             longLimits = c(-83, -78))

plot(oViewTransect3, "Ilr-transformed\nconcentration", span = 0.60, se = FALSE)

# Fourth transect
end_points4 <- data.frame(long = c(-79.5, -77.4),
                          lat = c(36, 33.8))

oViewTransect4 <- ViewTransect(spatialData, maxOffset, end_points4,
                               Example_CRS_arg_utm)

plotSampleLocations(oViewTransect4) +
  Add_Geography() +
  Add_Path(Example_domain, color = "blue") +
  Refine_Map("lambert", 25, 40,
             latLimits = c(33.5, 36),
             longLimits = c(-81,-77))

plot(oViewTransect4, "Ilr-transformed\nconcentration", span = 0.60, se = FALSE)

rm("maxOffset", "end_points1", "end_points2", "end_points3", "end_points4",
   "oViewTransect1", "oViewTransect2", "oViewTransect3", "oViewTransect4")


# Separation (in kilometers) between the centers of the basis functions
spacings <- c("020", "025", "030", "035", "040",
              "045", "050", "055", "060",
              "070", "080", "090", "100")

spacings <- "030"

# Construct a class that generates the basis functions
oBasisFunctions <- BasisFunctions(spacings,
                                  Example_domain,
                                  spatialData,
                                  Example_CRS_arg_utm,
                                  seed = 777)

# Print summary statistics of the basis functions for all spacings
summary(oBasisFunctions)

# Plot the locations of the basis function centers (with a spacing of
# 20 kilometers), as well as the domain and the region of interest
plot(oBasisFunctions, 20) +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

plot(oBasisFunctions, 30) +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

plot(oBasisFunctions, 50) +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

plot(oBasisFunctions, 70) +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

plot(oBasisFunctions, 100) +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

rm("spacings")

##
## Conduct the cross validation test
##

# Determine the file name (including the complete path) for
# Stan program BMNUS
filename <- paste0(normalizePath(path.package("BMNUS")), "\\stan\\BMNUS.stan")

# Translate Stan program BMNUS into C++ functions and then compile
# these functions
tr <- stanc(file = filename, model_name = "BMNUS")
BMNUS_sm <- stan_model(stanc_ret = tr, verbose=FALSE)

nFolds <- 10
nCpuCores <- 5

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

##
## Section "Sample the Posterior Probability Density Function"
##

# Prepare the data
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

# Function gen_int initializes some variables for the sampling
gen_init <- function(){
  phi1 <- rnorm(stanData$M, mean = 0.0, sd = 1e-6)
  phi2 <- rnorm(stanData$M, mean = 0.0, sd = 1e-6)
  return(list(phi1 = phi1, phi2 = phi2))
}

nChains <- 3
if(nChains < parallel::detectCores()){
  oldPar <- options(mc.cores = nChains)
} else{
  oldPar <- options(mc.cores = parallel::detectCores()-1)
}


rawSamples <- sampling(BMNUS_sm, data = stanData, chains = nChains,
                       iter = 2000, warmup = 500,
                       refresh = 100,
                       seed = 7,
                       init = gen_init,
                       pars=c("phi1", "tau1", "alpha1",
                              "phi2", "tau2", "alpha2", "rho",
                              "Y_mean", "Y_sd"))

options(oldPar)

# The following warning message may appear.
#
# Warning message:
#   In throw_sampler_warnings(nfit) :
#   Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess
#
# This warning message is caused by low
# effective sample size for one or more of the following model parameters:
# tau2, alpha2, tau1, alpha1, or row. Usually this message may be ignored.


# Archive the samples for subsequent analysis
if(!dir.exists("SamplingResults")){
  dir.create("SamplingResults")
}

save(stanData, rawSamples, file = "SamplingResults\\RawSamples.dat")

rm("csrBasisFuncs", "carQuantities", "gen_init", "nChains", "oldPar", "BMNUS_sm")

##
## Section "Check Convergence"
##

# print a summary of the sampling to a text file
filename <- paste("SamplingResults\\Summary.txt", sep = "")
oldPars <- options(width = 140, max.print = 100000)
sink(filename)

print(rawSamples, digits = 3,
      pars = c("tau1", "alpha1", "tau2", "alpha2", "rho"))

sink()
options(width = oldPars$width, max.print = oldPars$max.print)

rm("filename", "oldPars")


# additional way to assess convergence
library(shinystan)
launch_shinystan(rawSamples)

##
## Section "Check the Statistical Model"
##

# Check #1---analyze the fit between the model and the data
# along one or more transects

# Location of end points (specified with latitude and longitude)
# of the transect
end_points <- data.frame(long = c(-83.48355, -81.21104),
                         lat = c(33.30348, 29.69067))

# Maximum distance (in kilometers) between the transect and a field sample
maxOffset <- 10

# Construct the class for the transect
oCheckTransect <- CheckTransect(spatialData, maxOffset, end_points,
                                Example_CRS_arg_utm)

# Plot the transect and the field samples that are associated with it
plotSampleLocations(oCheckTransect) +
  Add_Geography() +
  Add_Path(Example_domain, color = "blue") +
  Refine_Map("lambert", 25, 40,
             latLimits = c(29, 34),
             longLimits = c(-85.5, -80))

# Plot four graphs to analyze the fit
plot(oCheckTransect, rawSamples, "Ilr-transformed\nconcentration")

# Check #2---compare, along a transect, the measured values to
# values that are simulated with the model.

# Plot eight graphs: 1 graph of the measured values along the transect
# and 7 graphs of simulated values along the transect
plotSimulatedTransects(oCheckTransect, rawSamples,
                       "Ilr-transformed\nconcentration")

# Check #3---analyze the statistical properties of the standardized
# residuals

# Construct the class for the standarized residuals
oCheckResiduals <- CheckResiduals(spatialData, rawSamples)

# Plot three graphs to analyze the standardized residuals
plot(oCheckResiduals, spatialData, Example_CRS_arg_utm,
     variogram_breaks = seq(from = 0, to = 150, length.out = 16))

# Print a table summarizing the standardized residuals
summary(oCheckResiduals)

# check #4---map the numerical signa of the standardized residuals
plotResidualMap(oCheckResiduals, spatialData, shape = 16, size = 0.9) +
  Add_Geography(boundary_color = "gray50") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5)) +
  Refine_Legend(c(0.02,0.99), c(0,1)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(legend.key = element_rect(fill = "gray90"))

rm("end_points", "maxOffset")

##
## Section "Map the Predicted Quantities"
##

# set the threshold for the map of exceedance probability
threshold <- scaledLogit(8000, constSumValue = 1e6)

# predict the quantities on a grid of uniformly-spaced points
# that cover the domain.
oPredictions <- Prediction(Example_domain,
                           spatialData,
                           Example_CRS_arg_utm,
                           basis_functions,
                           rawSamples,
                           threshold)

# plot one-third of the points
plotPredictionLocations(oPredictions, fraction = 1/4, fill = "black", shape = ".") +
  Add_Geography(boundary_color = "gray50") +
  Add_Path(Example_roi, "red") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5))

# Create a color palette for the maps
plotPalette <- colorRampPalette(c("blue", "green",
                                  "yellow", "orange", "red", "black"))(12)

# plot of the map of the process mean
plot(oPredictions, plotPalette, quantity = "Y_mean", isPaletteScaled = TRUE) +
  Add_Path(Example_roi, "black") +
  Add_Geography(boundary_color = "gray50") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5)) +
  Refine_Legend(c(0.02,0.99), c(0,1), legend.direction = "horizontal") +
  theme(plot.title=element_text(hjust = 0, face = "bold.italic", size = 10))

# plot of the map of the process standard deviation
plot(oPredictions, plotPalette, quantity = "Y_sd", isPaletteScaled = FALSE) +
  Add_Path(Example_roi, "black") +
  Add_Geography(boundary_color = "gray50") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5)) +
  Refine_Legend(c(0.02,0.99), c(0,1), legend.direction = "horizontal") +
  theme(plot.title=element_text(hjust = 0, face = "bold.italic", size = 10))

# plot of the map of the exceedance probability
plot(oPredictions, plotPalette, quantity = "prob", isPaletteScaled = FALSE) +
  Add_Path(Example_roi, "black") +
  Add_Geography(boundary_color = "gray50") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5)) +
  Refine_Legend(c(0.02,0.99), c(0,1), legend.direction = "horizontal") +
  theme(plot.title=element_text(hjust = 0, face = "bold.italic", size = 10))

# plot of the map of the quantile
plot(oPredictions, plotPalette, quantity = "quantile", isPaletteScaled = TRUE) +
  Add_Path(Example_roi, "black") +
  Add_Geography(boundary_color = "gray50") +
  Refine_Map("lambert", 25, 40, latLimits = c(28, 39.5)) +
  Refine_Legend(c(0.02,0.99), c(0,1), legend.direction = "horizontal") +
  theme(plot.title=element_text(hjust = 0, face = "bold.italic", size = 10))

# Re-express the axis labels for the process mean as titanium concentration
# with units of milligrams per kilogram
invScaledLogit(c(-4.5, -4.0, -3.5), constSumValue = 1e6)

##
## Scripts for Appendix 1
##

library(tibble)
library(BMNUS)
library(PairedMeasurements)

measurements <- tibble(first = scaledLogit(PairedMeasurements::ExampleData$first, constSumValue = 1.0e6),
                       second = scaledLogit(PairedMeasurements::ExampleData$second, constSumValue = 1.0e6))

oEvaluation <- PM_Evaluation(measurements, "ilr-transformed concentration")
plot(oEvaluation)
summary(oEvaluation)




