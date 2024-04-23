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
data = read.csv('GaN2023.csv')

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
lm_model <- lm(SQMReading ~ Elevation - 1 , data = f_data)
summary(lm_model)
var(lm_model$residuals)

#Semi-variogram of the residuals
Ext_cp <- f_data %>%
  select(y_coord=Latitude, x_coord=Longitude, y=SQMReading) %>%
  drop_na()

sv_residual <- Ext_cp %>% with(geoR::variog(data = lm_model$residuals, coords = cbind(x_coord, y_coord), messages=FALSE))

sv_res_df <- data.frame(dists = sv_residual$u, variogram = sv_residual$v, npairs = sv_residual$n, sd = sv_residual$sd)

sv_res_plot <- ggplot(sv_res_df, aes(x=dists, y=variogram)) +
  geom_point(size=2, shape=8) +
  ggtitle("Semivariogram of the Residuals") +
  theme_minimal()

sv_res_plot

# Bayes Model
library(spBayes)




X <- f_data %>%  select(Elevation) %>% as.matrix()

Y <- f_data$SQMReading %>% as.matrix()

sample <- sample.int(length(Y), size=floor(.8*length(Y)),replace=FALSE)

Y_train <- Y[sample,] %>% as.numeric()
X_train <- X[sample,] 
Y_test <- Y[-sample,] %>% as.numeric()
X_test <- X[-sample,] %>% as.matrix()



coords <- cbind(f_data$Latitude, f_data$Longitude)
dup <- duplicated(coords)
coords[dup] <- coords[dup] + runif(sum(dup), 0, 1e-2)

n.samples <- 20000

coords_train <- coords[sample,]
coords_test <- coords[-sample,]

starting <- list("phi"=15, "sigma.sq"=145, "tau.sq"=1)
tuning <- list("phi"=0.03, "sigma.sq"=0.03, "tau.sq"=0.03)

priors <- list("beta.Norm"=list(rep(0,ncol(X)), diag(1000,ncol(X))),
               "phi.Unif"=c(10/1, 10/0.1), "sigma.sq.IG"=c(2, 200),
               "tau.sq.IG"=c(3, 300))

cov.model <- "exponential"

system.time({
  model_fit <- spLM(Y_train~X_train-1, coords=coords_train, starting=starting,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=n.samples, verbose=FALSE, n.report=500)
})


model_fit$acceptance

png(filename = "Theta_Trace.png", res=800, width=5000, height=3000)

par(mfrow=c(2,2))

ts.plot(model_fit$p.theta.samples[,1],main="sigmasq",ylab="", xlim=c(40000,nrow(model_fit$p.theta.samples)),ylim=c(110,200))

ts.plot(model_fit$p.theta.samples[,2],main="tausq",ylab="", xlim=c(40000,nrow(model_fit$p.theta.samples)),ylim=c(1.5,2.8))

ts.plot(model_fit$p.theta.samples[,3],main="phi",ylab="", xlim=c(40000,nrow(model_fit$p.theta.samples)),ylim=c(9.9,10.4))

dev.off()

library(gridExtra)
burn.in<-0.75*n.samples
model_recover<-spRecover(model_fit,start=burn.in)
w_post <- model_recover$p.w.recover.samples
b_post <- model_recover$p.beta.recover.samples
XB_post <- X_train %*% t(b_post)
XBW_mean <- (XB_post + w_post) %>% apply(1, mean)
model_residual <- Y_train - XBW_mean



png(filename = "Residual_Bayes.png" ,res=800, width=6000, height=3000)
df <- f_data %>% mutate(residual = model_residual)
residual_map <- ggplot(df, aes(Latitude, Longitude, color=residual)) +
  geom_point() + 
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Residual")

sv <- df %>% with(geoR::variog(data = residual, coords = coords, messages=FALSE))

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

plot(Y_test, y.hat[2,], pch=19, cex=0.5, xlab="observed y", ylab="predicted y")
arrows(Y_test[-mod], y.hat[2,-mod], Y_test, y.hat[1,-mod], angle=90, length=0.05)
arrows(Y_test[-mod], y.hat[2,-mod], Y_test, y.hat[3,-mod], angle=90, length=0.05)

#Low rank GP

system.time({
  lrmodel_fit <- spLM(Y~X-1, coords=coords, starting=starting,
                    tuning=tuning, priors=priors, cov.model=cov.model,
                    n.samples=n.samples, verbose=FALSE, n.report=500,knots= c(6,6,0.1))
})

# Nearest Neighbor GP
library(spNNGP)
