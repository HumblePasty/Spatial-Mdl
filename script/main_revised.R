library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(tidyr)
library(GGally)
library(coda)
library(scico)
library(MBA)
library(viridis)
library(fields)
# work dir
# setwd("./data")
setwd("./BIOSTAT 696/Final")
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
X <- scale(X)

Y <- Ext_cp$y %>% as.matrix()
Y <- scale(Y)

sample <- sample.int(length(Y), size=floor(.8*length(Y)),replace=FALSE)

Y_train <- Y[sample,] %>% as.numeric()
X_train <- X[sample,] %>% as.matrix()
Y_test <- Y[-sample,] %>% as.numeric()
X_test <- X[-sample,] %>% as.matrix()



coords <- cbind(Ext_cp$y_coord,Ext_cp$x_coord)
dup <- duplicated(coords)
coords[dup] <- coords[dup] + runif(sum(dup), 0, 1e-2)

n.samples <- 50000

coords_train <- coords[sample,]
coords_test <- coords[-sample,]

starting <- list("phi"=1/0.5, "sigma.sq"=2, "tau.sq"=1)
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

png(filename = "Theta_Trace_sigma_Bayes.png", res=400, width=2000, height=1500)


ts.plot(model_fit$p.theta.samples[,1],main="sigmasq",ylab="", xlim=c(40000,nrow(model_fit$p.theta.samples)),ylim=c(0,1))

dev.off()

png(filename = "Theta_Trace_tau_Bayes.png", res=400, width=2000, height=1500)

ts.plot(model_fit$p.theta.samples[,2],main="tausq",ylab="", xlim=c(40000,nrow(model_fit$p.theta.samples)),ylim=c(0.1,0.25))

dev.off()

png(filename = "Theta_Trace_phi_Bayes.png", res=400, width=2000, height=1500)
ts.plot(model_fit$p.theta.samples[,3],main="phi",ylab="", xlim=c(40000,nrow(model_fit$p.theta.samples)),ylim=c(8,10.5))

dev.off()

mc_test <- as.mcmc(model_fit$p.theta.samples)
traceplot(mc_test,xlim=c(15000,nrow(model_fit$p.theta.samples)),col=1:2)

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
df <- train %>% mutate(residual = model_residual, Longitude = V5, Latitude = V6)
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

png(filename = "Compare_LM_Bayes.png" ,res=800, width=6000, height=3000)
par(mfrow=c(1,2))
lm.resids <- resid(lm_model)
resid.surf <- mba.surf(cbind(coords, lm.resids), no.X=200, no.Y=200, extend=TRUE)$xyz.est
image.plot(resid.surf, col=viridis(100), main="LM residuals\n(interpolated surface)")

quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

w.summary <- apply(model_recover$p.w.recover.samples, 1, quants)
w.mu.surf <- mba.surf(cbind(coords_train, w.summary[1,]), no.X=200, no.Y=200, extend=TRUE)$xyz.est
image.plot(w.mu.surf, col=viridis(100), main="spBayes()\nspatial random effects, w")
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

mcmc <- as.mcmc(model_fit$p.theta.samples)

model

summary(mcmc)

# Define different starting values for each chain
starting_values <- list(
  list("phi" = 1/0.5, "sigma.sq" = 2, "tau.sq" = 1),   # Original starting values
  list("phi" = 1/0.3, "sigma.sq" = 4, "tau.sq" = 2),  # Modified starting values for second chain
  list("phi" = 1/0.7, "sigma.sq" = 6, "tau.sq" = 0.5)  # Modified starting values for third chain
)

# Run multiple chains
model_fits <- lapply(starting_values, function(s) {
  spLM(Y_train ~ X_train - 1, coords = coords_train, starting = s,
       tuning = tuning, priors = priors, cov.model = cov.model,
       n.samples = 20000)
})

# Combine the MCMC samples from each chain into a single mcmc.list object
mcmc_list <- mcmc.list(lapply(model_fits, function(fit) as.mcmc(fit$p.theta.samples)))

png(filename = "GelmanDiag_Bayes.png" ,res=800, width=6000, height=5000)
gelman.plot(mcmc_list)
dev.off()

gelman.diag(mcmc_list)


# Predictions on new locations
pred.coords <- expand.grid(seq(36,38,0.05),seq(-81,-78, 0.1))
m <- nrow(pred.coords)
pred.covars <- mkMvX(list(matrix(1,m,1), matrix(1,m,1), matrix(1,m,1)))

out <- spPredict(model_recover,start = 200, pred.coords = pred.coords, pred.covars = pred.covars,n.omp.threads = 8, thin = 20, joint = TRUE)

colnames(pred.coords) <- c("Latitude", "Longitude")

predicted_values <- apply(out$p.y.predictive.samples, 1, mean)

prediction_data <- cbind(pred.coords, Predicted_Values = predicted_values)

prediction_data <- as.data.frame(prediction_data)
colnames(prediction_data) <- c("Latitude", "Longitude", "Predicted_Values")

ggplot(Ext_cp, aes(x = x_coord, y = y_coord)) +
  geom_raster(data = prediction_data, aes(x = Longitude, y = Latitude, fill=Predicted_Values)) +  # Use geom_tile() if data represents a grid
  scale_fill_scico(palette="batlowK") +  # Apply a color scale
  labs(title = "Spatial Prediction Results", x = "Longitude", y = "Latitude") +
  theme_minimal()

y_pred <- rowMeans(out$p.y.predictive.samples)

par(mfrow=c(1,2))
Ext_cp.crop <- Ext_cp[Ext_cp$y_coord>36 & Ext_cp$y_coord <38 & Ext_cp$x_coord > -81 & Ext_cp$x_coord < -78,]
coords_crop <- cbind(Ext_cp.crop$x_coord, Ext_cp.crop$y_coord)
Y_crop <- Ext_cp.crop$y
surf <- mba.surf(cbind(coords_crop,Y_crop),no.X=100, no.Y=100,extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, add=TRUE)
surf <- mba.surf(cbind(pred.coords , y_pred),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Predicted counts")
contour(surf, add=TRUE)
points(pred.coords,pch=19,cex=0.5,col="blue")



#Low rank GP

system.time({
  lrmodel_fit <- spLM(Y_train~X_train-1, coords=coords_train, starting=starting,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=20000, verbose=FALSE, n.report=500,knots= c(6,6,0.1))
})


burn.in<-0.75*20000
lrmodel_recover<-spRecover(lrmodel_fit,start=burn.in)
w_post <- lrmodel_recover$p.w.recover.samples
b_post <- lrmodel_recover$p.beta.recover.samples
XB_post <- X_train %*% t(b_post)
XBW_mean <- (XB_post + w_post) %>% apply(1, mean)
lrmodel_residual <- Y_train - XBW_mean

rmse.lr <- sqrt(mean(lrmodel_residual^2))
rmse.lr

lrtrain <- cbind(Y_train,X_train,coords_train) %>% as.data.frame()

png(filename = "Low Rank Residual_Bayes.png" ,res=800, width=6000, height=3000)
df <- lrtrain %>% mutate(residual = lrmodel_residual, Longitude = V6, Latitude = V5)
residual_map <- ggplot(df, aes(Latitude, Longitude, color=residual)) +
  geom_point() + 
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Low-Rank Residual")

sv <- df %>% with(geoR::variog(data = residual, coords = coords_train, messages=FALSE))

sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=2, shape=8) +
  theme_minimal()

m <-grid.arrange(residual_map, sv_plot, nrow=1)

dev.off()

# Define different starting values for each chain
starting_values <- list(
  list("phi" = 1/0.5, "sigma.sq" = 2, "tau.sq" = 1),   # Original starting values
  list("phi" = 1/0.3, "sigma.sq" = 4, "tau.sq" = 2),  # Modified starting values for second chain
  list("phi" = 1/0.7, "sigma.sq" = 6, "tau.sq" = 0.5)  # Modified starting values for third chain
)

# Run multiple chains
lrmodel_fits <- lapply(starting_values, function(s) {
  spLM(Y_train ~ X_train - 1, coords = coords_train, starting = s,
       tuning = tuning, priors = priors, cov.model = cov.model,
       n.samples = 20000,knots= c(6,6,0.1))
})

# Combine the MCMC samples from each chain into a single mcmc.list object
lrmcmc_list <- mcmc.list(lapply(lrmodel_fits, function(fit) as.mcmc(fit$p.theta.samples)))

png(filename = "GelmanDiag_LRBayes.png" ,res=800, width=6000, height=5000)
gelman.plot(lrmcmc_list)
dev.off()

gelman.diag(lrmcmc_list)


#Predictions
lrpred <- spPredict(lrmodel_recover, pred.covars = X_test, pred.coords = coords_test, start = 0.5* n.samples)

y.hat.lr <- apply(lrpred$p.y.predictive.samples, 1, mean)

quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}

y.hat.lr <- apply(lrpred$p.y.predictive.samples, 1, quant)

png(filename = "PredictvsTest_LRBayes.png" ,res=800, width=6000, height=3000)
# Convert your predictions and observed data to a data frame for ggplot
data_for_plot <- data.frame(
  observed = Y_test,
  predicted = y.hat.lr[2,],
  lower = y.hat.lr[1,],
  upper = y.hat.lr[3,]
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


# Nearest Neighbor GP
library(spNNGP)

tuning.NN <- list("phi"=2, "sigma.sq"=0.05, "tau.sq"=0.05)
system.time({
nnmodel_fit.latent <- spNNGP(Y_train~X_train-1, coords=coords_train, starting=starting, method="response", n.neighbors=10,
                       tuning=tuning, priors=priors, cov.model=cov.model,
                       n.samples=n.samples, n.omp.threads=4,fit.rep=TRUE,sub.sample=list(start=4000,thin=10))
})

png(filename = "Trace_NNBayes(Response).png" ,res=800, width=6000, height=6000)
plot(nnmodel_fit.latent$p.theta.samples)
dev.off()


#MCMC diagnosis


# Define different starting values for each chain
starting_values <- list(
  list("phi" = 1/0.5, "sigma.sq" = 0.2, "tau.sq" = 0.1),   # Original starting values
  list("phi" = 1/0.3, "sigma.sq" = 0.4, "tau.sq" = 0.2),  # Modified starting values for second chain
  list("phi" = 1/0.7, "sigma.sq" = 0.6, "tau.sq" = 0.5)  # Modified starting values for third chain
)

# Run multiple chains
model_fits <- lapply(starting_values, function(s) {
  spNNGP(Y_train ~ X_train - 1, coords = coords_train, starting = s,
        tuning=tuning.NN, priors = priors, cov.model = cov.model,
       n.neighbors=10, method='latent',
       n.samples=20000, n.omp.threads=4,fit.rep=TRUE,sub.sample=list(start=4000,thin=10))
})

# Combine the MCMC samples from each chain into a single mcmc.list object
mcmc_list <- mcmc.list(lapply(model_fits, function(fit) as.mcmc(fit$p.theta.samples)))

png(filename = "GelmanDiag_NNBayes(Response).png" ,res=800, width=6000, height=5000)
gelman.plot(mcmc_list)
dev.off()

gelman.diag(mcmc_list)

#Residuals
XBW_mean <- nnmodel_fit.latent$y.hat.samples %>% apply(1, mean)
nnmodel_residual.latent <- Y_train - XBW_mean

rmse.nn <- sqrt(mean(nnmodel_residual.latent^2))
rmse.nn

nntrain <- cbind(Y_train,X_train,coords_train) %>% as.data.frame()

png(filename = "NN_Residual_Bayes(Response).png" ,res=800, width=6000, height=3000)
df <- nntrain %>% mutate(residual = nnmodel_residual.latent, Longitude = V5, Latitude = V6)
residual_map <- ggplot(df, aes(Latitude, Longitude, color=residual)) +
  geom_point() + 
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Nearest Neighbor Residual(Response)")

sv <- df %>% with(geoR::variog(data = residual, coords = coords_train, messages=FALSE))

sv_df <- data.frame(dists = sv$u, variogram = sv$v, npairs = sv$n, sd = sv$sd)
sv_plot <- ggplot(sv_df, aes(x=dists, y=variogram)) + geom_point(size=2, shape=8) +
  theme_minimal()

m <-grid.arrange(residual_map, sv_plot, nrow=1)

dev.off()



#Predictions
nnpred <- predict(nnmodel_fit.latent, X.0 = X_test, coords.0 = coords_test,n.omp.threads = 4, n.report=1000)

y.hat.lr <- apply(lrpred$p.y.predictive.samples, 1, mean)

quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}

y.hat.lr <- apply(lrpred$p.y.predictive.samples, 1, quant)

png(filename = "PredictvsTest_LRBayes.png" ,res=800, width=6000, height=3000)
# Convert your predictions and observed data to a data frame for ggplot
data_for_plot <- data.frame(
  observed = Y_test,
  predicted = y.hat.lr[2,],
  lower = y.hat.lr[1,],
  upper = y.hat.lr[3,]
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
