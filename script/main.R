library(ggplot2)
library(gridExtra)
library(reshape2)

# work dir
setwd("D:/UMich/Win24/BIOSTAT 696/FinalProject")

# load the data
data = read.csv('D:/UMich/Win24/BIOSTAT 696/FinalProject/data/GaN2023.csv')

# filter the data with "SQMReading" between 16 and 22
# drop NA values
data = data[!is.na(data$SQMReading),]
f_data = data[data$SQMReading >= 16 & data$SQMReading <= 22,]

# plot the locations of the data
p1 = ggplot(f_data, aes(x = Longitude, y = Latitude)) +
  geom_point() +
  labs(title = "Locations of the data")

# crop the data within USA
f_data = f_data[f_data$Longitude >= -125 & f_data$Longitude <= -65 & f_data$Latitude >= 15,]
p1

# plot a semivarogram
Ext_cp <- f_data %>%
  select(y_coord=Latitude, x_coord=Longitude, y=SQMReading) %>%
  drop_na()

sv2 <- Ext_cp %>% with(geoR::variog(data = y, coords = cbind(x_coord, y_coord), messages=FALSE))

sv2_df <- data.frame(dists = sv2$u, variogram = sv2$v, npairs = sv2$n, sd = sv2$sd)
sv2_plot <- ggplot(sv2_df, aes(x=dists, y=variogram)) +
  geom_point(size=2, shape=8) +
  ggtitle("Semivariogram of the Point Reference Data") +
  theme_minimal()
sv2_plot
