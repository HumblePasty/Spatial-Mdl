library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(tidyr)
library(GGally)

# work dir
setwd("./data")
# setwd("./BIOSTAT 696/Final")
# load the data
data = read.csv('GaN2023_full.csv')

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


# Bayesian Modelling
lm_model <- lm(SQMReading ~ Elevation + popden_km2 + elec_use -1 , data = f_data)
summary(lm_model)
var(lm_model$residuals)

#Semi-variogram of the residuals
Ext_cp <- f_data %>%
  select(y_coord=Latitude, x_coord=Longitude, y=SQMReading, x1 = Elevation, x2= popden_km2, x3 = elec_use) %>%
  drop_na()

sv_residual <- Ext_cp %>% with(geoR::variog(data = lm_model$residuals, coords = cbind(x_coord, y_coord), messages=FALSE))

sv_res_df <- data.frame(dists = sv_residual$u, variogram = sv_residual$v, npairs = sv_residual$n, sd = sv_residual$sd)

sv_res_plot <- ggplot(sv_res_df, aes(x=dists, y=variogram)) +
  geom_point(size=2, shape=8) +
  ggtitle("Semivariogram of the Residuals") +
  theme_minimal()

sv_res_plot

ini.vals <- expand.grid(seq(0.1,1,l=5), seq(0.1,1,l=5))
variofit <- geoR::variofit(sv_residual, ini=ini.vals, fix.nugget=TRUE, cov.model="exponential")

plot(sv_residual)
lines(variofit)

print(c(variofit$cov.pars[1],variofit$cov.pars[2]))

# Bayes Model
library(spBayes)




X <- Ext_cp %>%  select(c(x1,x2,x3)) %>% as.matrix()
X[,3] <- X[,3]/10000

Y <- Ext_cp$y %>% as.matrix()

sample <- sample.int(length(Y), size=floor(.8*length(Y)),replace=FALSE)

Y_train <- Y[sample,] %>% as.numeric()
X_train <- X[sample,] %>% as.matrix()
Y_test <- Y[-sample,] %>% as.numeric()
X_test <- X[-sample,] %>% as.matrix()



coords <- cbind(Ext_cp$y_coord,Ext_cp$x_coord)
dup <- duplicated(coords)
coords[dup] <- coords[dup] + runif(sum(dup), 0, 1e-2)

n.samples <- 20000

coords_train <- coords[sample,]
coords_test <- coords[-sample,]

starting <- list("phi"=1/0.5, "sigma.sq"=70, "tau.sq"=1)
tuning <- list("phi"=0.05, "sigma.sq"=0.05, "tau.sq"=0.05)

priors <- list("beta.Norm"=list(rep(0,ncol(X_train)), diag(1000,ncol(X_train))),
               "phi.Unif"=c(1/1, 1/0.1), "sigma.sq.IG"=c(2, 2),
               "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"


system.time({
  model_fit <- spLM(Y_train~X_train-1, coords=coords_train, starting=starting,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=n.samples, verbose=FALSE, n.report=500)
})


model_fit$acceptance

png(filename = "Theta_Trace.png", res=800, width=5000, height=3000)

par(mfrow=c(2,2))

ts.plot(model_fit$p.theta.samples[,1],main="sigmasq",ylab="", xlim=c(15000,nrow(model_fit$p.theta.samples)),ylim=c(85,135))

ts.plot(model_fit$p.theta.samples[,2],main="tausq",ylab="", xlim=c(15000,nrow(model_fit$p.theta.samples)),ylim=c(0.10,0.3))

ts.plot(model_fit$p.theta.samples[,3],main="phi",ylab="", xlim=c(15000,nrow(model_fit$p.theta.samples)),ylim=c(1,1.05))

dev.off()

library(gridExtra)
burn.in<-0.75*n.samples
model_recover<-spRecover(model_fit,start=burn.in)
w_post <- model_recover$p.w.recover.samples
b_post <- model_recover$p.beta.recover.samples
XB_post <- X_train %*% t(b_post)
XBW_mean <- (XB_post + w_post) %>% apply(1, mean)
model_residual <- Y_train - XBW_mean

train <- cbind(Y_train,X_train,coords_train) %>% as.data.frame()

png(filename = "Residual_Bayes.png" ,res=800, width=6000, height=3000)
df <- train %>% mutate(residual = model_residual, Longitude = V6, Latitude = V5)
residual_map <- ggplot(df, aes(Latitude, Longitude, color=residual)) +
  geom_point() + 
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Residual")

sv <- df %>% with(geoR::variog(data = residual, coords = coords_train, messages=FALSE))

sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=2, shape=8) +
  theme_minimal()

m <-grid.arrange(residual_map, sv_plot, nrow=1)

dev.off()




#Predictions
pred <- spPredict(model_recover, pred.covars = X_test, pred.coords = coords_test, start = 0.5* n.samples)

y.hat <- apply(pred$p.y.predictive.samples, 1, mean)

quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}

y.hat <- apply(pred$p.y.predictive.samples, 1, quant)

png(filename = "PredictvsTest_Bayes.png" ,res=800, width=6000, height=3000)
# Convert your predictions and observed data to a data frame for ggplot
data_for_plot <- data.frame(
  observed = Y_test,
  predicted = y.hat[2,],
  lower = y.hat[1,],
  upper = y.hat[3,]
)

ggplot(data_for_plot, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.6) +  # Add points with slight transparency
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, alpha = 0.6) +  # Add error bars for intervals
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +  # Add a 1:1 reference line
  theme_minimal() +  # Use a minimal theme
  labs(x = "Observed Y", y = "Predicted Y", title = "Predicted vs Observed Y") + 
  theme(
    text = element_text(size = 12),  # Change text size
    plot.title = element_text(face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold")  # Bold axis titles
  )

dev.off()

pred_points <- expand.grid(list(Longitude=seq(min(coords[,2]), max(coords[,2]), length.out=50),
                                Latitude=seq(min(coords[,1]), max(coords[,1]), length.out=50))) 
outcoords <- as.matrix(pred_points)


#MCMC diagnosis
library(coda)
mcmc <- as.mcmc(model_fit$p.theta.samples)

model

summary(mcmc)

# Define different starting values for each chain
starting_values <- list(
  list("phi" = 1/0.5, "sigma.sq" = 70, "tau.sq" = 1),   # Original starting values
  list("phi" = 1/0.3, "sigma.sq" = 100, "tau.sq" = 2),  # Modified starting values for second chain
  list("phi" = 1/0.7, "sigma.sq" = 50, "tau.sq" = 0.5)  # Modified starting values for third chain
)

# Define different starting values for each chain
starting_values <- list(
  list("phi" = 1/0.5, "sigma.sq" = 70, "tau.sq" = 1),   # Original starting values
  list("phi" = 1/0.3, "sigma.sq" = 100, "tau.sq" = 2),  # Modified starting values for second chain
  list("phi" = 1/0.7, "sigma.sq" = 50, "tau.sq" = 0.5)  # Modified starting values for third chain
)

# Run multiple chains
model_fits <- lapply(starting_values, function(s) {
  spLM(Y_train ~ X_train - 1, coords = coords_train, starting = s,
       tuning = tuning, priors = priors, cov.model = cov.model,
       n.samples = 20000, n.burnin = 2000, n.thin = 10)
})

# Combine the MCMC samples from each chain into a single mcmc.list object
mcmc_list <- mcmc.list(lapply(model_fits, function(fit) as.mcmc(fit$p.theta.samples)))


gelman.plot(mcmc)

#Low rank GP

system.time({
  lrmodel_fit <- spLM(Y~X-1, coords=coords, starting=starting,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=n.samples, verbose=FALSE, n.report=500,knots= c(6,6,0.1))
})

# Nearest Neighbor GP
library(spNNGP)
