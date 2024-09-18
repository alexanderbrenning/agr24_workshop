####################################################
# Landslide Susceptibility Modeling with R
# Workshop at the conference of the
# Romanian Association of Geomorphologists
# Iasi, 19 September 2024
# (c) Alexander Brenning - University of Jena
####################################################
# Statistical learning for classification
# Case study: landslides in Ecuador
####################################################

##############################################
# Prepare the training and test samples:
# Transform some of the variables
##############################################

# This is where my raster data from Ecuador is located,
# relative to the RStudio project's folder:
setwd("data")

library("pacman")
p_load("ROCR")
p_load("mgcv")
p_load("rpart")
p_load("randomForest")
p_load("RSAGA")

# Load the saved training/test data:
d <- readRDS("landslides.rds")
d$slides89 <- factor(d$slides89)

# A function for plotting ROC curves and
# calculating the area under the ROC curve (AUROC)
# using the ROCR package:
auroc <- function(pred, obs, plot=FALSE) {
  stopifnot(is.logical(obs) | is.factor(obs))
  stopifnot(is.numeric(pred))
  stopifnot(length(pred) == length(obs))
  if (is.factor(obs)) stopifnot(nlevels(obs)==2)
  require(ROCR)
  predobj <- prediction( pred, obs )
  if (plot) plot(performance(predobj, "tpr", "fpr"))
  auroc <- performance( predobj, measure="auc" )@y.values[[1]]
  return(auroc)
}

# Plot a nice landslide susceptibility map:
nice_map <- function(file, title = "Landslide susceptibility", 
                     hillshade = "hillshade.asc",
                     pts = NULL) {
  require("terra")
  # Make nice plot of results
  plot(terra::rast(hillshade), 
       col = gray.colors(10, start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL),
       legend = FALSE, main = title)
  plot(rast(file), alpha = 0.60, add=TRUE)
  if (!is.null(pts)) points(pts, cex = 0.3)
  invisible(NULL)
}

# Apply some transformations to the data:
my.trafo <- function(x) {
  # Same for >300m from deforestation:
  x$distdeforest[ x$distdeforest > 300 ] <- 300
  # ...and >300m from road:
  x$distroad[ x$distroad > 300 ] <- 300
  # Convert radians to degrees:
  x$slope <- x$slope * 180 / pi
  x$cslope <- x$cslope * 180 / pi
  # Log-transform size of catchment area - it is extremely skewed:
  x$log.carea <- log10(x$carea)
  # Plan and profile curvature have outliers -> trim both variables:
  x$plancurv[ x$plancurv < -0.2 ] <- -0.2
  x$plancurv[ x$plancurv >  0.2 ] <-  0.2
  x$profcurv[ x$profcurv < -0.04 ] <- -0.04
  x$profcurv[ x$profcurv >  0.04 ] <-  0.04
  return(x)
}

d <- my.trafo(d)


##############################################
# Exploratory Data Analysis
##############################################

# Spinograms are a very fine tool for displaying
# conditional probabilities:
par(mfrow=c(2,4),mex=0.7,mar=c(5,4,2,3))
with(d, {
  #spineplot(slides89 ~ distslidespast, ylevels = c("TRUE", "FALSE"), xlab = "Dist. from past landslides", ylab = "", yaxlabels = c("Landslide", "No lsl."),
  #        breaks = c(0,10,40,80,100))
  spineplot(slides89 ~ slope, ylevels = c("TRUE", "FALSE"), xlab = "Slope angle", ylab = "", yaxlabels = c("Landslide", "No landsl."))
  spineplot(slides89 ~ plancurv, ylevels = c("TRUE", "FALSE"), xlab = "Plan curvature", ylab = "", yaxlabels = c("Landslide", "No landsl."))
  spineplot(slides89 ~ profcurv, ylevels = c("TRUE", "FALSE"), xlab = "Profile curvature", ylab = "", yaxlabels = c("Landslide", "No landsl."))
  spineplot(slides89 ~ log.carea, ylevels = c("TRUE", "FALSE"), xlab = "Log. contributing area", ylab = "", yaxlabels = c("Lsl.", "No landslide"))
  spineplot(slides89 ~ cslope, ylevels = c("TRUE", "FALSE"), xlab = "Catchment slope", ylab = "", yaxlabels = c("Landslide", "No landsl."))
  spineplot(slides89 ~ distdeforest, ylevels = c("TRUE", "FALSE"), xlab = "Dist. to deforestation", ylab = "", yaxlabels = c("Landslide", "No lsl."),
            breaks = c(0,25,50,100,250,300))
  spineplot(slides89 ~ distroad, ylevels = c("TRUE", "FALSE"), xlab = "Distance to road", ylab = "", yaxlabels = c("Landslide", "No lsl."),
            breaks = c(0,25,50,100,250,300))
} )



##############################################
# Generalized Additive Model
##############################################

# Fit a generalized additive model (GAM) with a logistic
# link function by stepwise forward variable selection.
# GAMs are semi-parametric models that combine linear and
# nonlinear terms; they are more flexible than the generalized
# linear model, but still interpretable because of their
# additive structure.

p_load("mgcv")
# Model building with a shrinkage penalty that is capable of penalizing predictors
# out of the model:
fit <- gam(slides89 ~ s(slope, k = 7) + s(plancurv) + s(profcurv, k = 7) +
             s(log.carea) + s(cslope) + 
             I(distdeforest == 0) + s(distdeforest) + s(distroad),
           data = d, family = binomial, select = TRUE)


summary(fit)

par(mfrow = c(2,2))
gam.check(fit)
# Low k index for distdeforest and distroad, 
# but edf are not really close to k' = 9.

par(mfrow = c(3,3))
plot(fit)

# Create ROC plots and calculate AUROC on the training set:
par(mfrow=c(1,1))
pred <- as.vector(predict(fit, newdata = d, type = "response"))
auroc(pred, d$slides89=="TRUE", plot=TRUE )

# Apply fitted GAM to the entire study area:
library("RSAGA")
multi.local.function(
  in.grids = c("slope", "carea", "cslope", "plancurv", "profcurv", "distroad", "distdeforest"),
  out.varnames = "gampred",
  fun = grid.predict, 
  control.predict = list(type = "response"),
  fit = fit, trafo = my.trafo, quiet = FALSE )

# Plot a simple prediction map:
# p_load("raster")
# gam.raster <- raster("gampred.asc")
# plot(gam.raster)

# nicer verion of map:
nice_map("gampred.asc", title = "GAM prediction")


##############################################
# Classification Tree
##############################################

p_load("rpart")
fit <- rpart(slides89 ~ slope + plancurv + profcurv +
              log.carea + cslope + distdeforest + distroad, data = d,
            control = rpart.control(cp = 0.005))
par(xpd = TRUE, mfrow=c(1,1))
plot(fit)
text(fit, use.n = TRUE)
par(xpd = FALSE)

# Training set AUROC:
pred <- predict(fit, newdata = d, type = "prob")[ , "TRUE" ]
auroc(pred, d$slides89=="TRUE", plot=TRUE )

multi.local.function(
  in.grids = c("slope", "plancurv", "profcurv", "carea", "cslope", "distroad", "distdeforest"),
  out.varnames = "ctpred",
  fun = grid.predict, control.predict = list(type = "prob"),
  fit = fit, trafo = my.trafo, predict.column = "TRUE",
  quiet = FALSE )

nice_map("ctpred.asc", title = "Classification tree prediction")


##############################################
# Random Forest
##############################################

p_load("randomForest") # alternatives: 'ranger', 'cforest'
fit <- randomForest(slides89 ~ slope + plancurv + profcurv + log.carea + cslope + 
                     distdeforest + distroad, data = d, importance = TRUE)
fit

# Variable importance plot - will talk about this later:
varImpPlot(fit,type=1)

# Training set AUROC:
pred <- predict(fit, newdata = d, type = "prob")[ , "TRUE" ]
auroc(pred, d$slides89 == "TRUE", plot=TRUE )
# AUROC: 1.00!!! (i.e. nicely overfitting...)

hist(pred, br = 100)

d$rf_pred <- pred
d$rf_class <- pred > 0.5
# d[c(1700:1705), c("slides89", "rf_class")]
mean((d$slides89 == "TRUE") != d$rf_class)

multi.local.function(
  in.grids = c("slope", "plancurv", "profcurv", "carea", "cslope", "distroad", "distdeforest"),
  out.varnames = "rfpred",
  fun = grid.predict, control.predict = list(type = "prob"),
  fit = fit, trafo = my.trafo, predict.column = "TRUE",
  quiet = FALSE )

nice_map("rfpred.asc", title = "Random forest prediction")




##############################################
# Support Vector machine
##############################################

p_load("e1071") # see also 'kernlab'
fit <- svm(slides89 ~ slope + plancurv + profcurv + log.carea + cslope + 
                      distdeforest + distroad, data = d, 
           cost = 3,    #  default: 1
           gamma = 1/2, # default: 1/nfeatures, i.e. 1/7
           probability = TRUE)
fit

# pred <- predict(fit, newdata = d, probability = TRUE)
# predprob <- attr(pred, "probabilities")[,"TRUE"]
# hist(predprob, br = 100)

mypredictsvm <- function(object, newdata) {
  pred <- predict(object, newdata, probability = TRUE)
  unname(attr(pred, "probabilities")[,"TRUE"])
}

# Training set AUROC:
pred <- mypredictsvm(fit, newdata = d)
auroc(pred, d$slides89=="TRUE", plot=TRUE )


multi.local.function(
    in.grids = c("slope", "plancurv", "profcurv", "carea", "cslope", "distroad", "distdeforest"),
    out.varnames = "svmpred",
    fun = grid.predict, control.predict = list(),
    predfun = mypredictsvm,
    fit = fit, trafo = my.trafo,
    quiet = FALSE )

nice_map("svmpred.asc", title = "SVM prediction")


##############################################
# Putting it all together:
##############################################

par(mfrow = c(2,2))
nice_map("gampred.asc", title = "GAM prediction")
nice_map("svmpred.asc", title = "SVM prediction")
nice_map("ctpred.asc", title = "Classification tree prediction")
nice_map("rfpred.asc", title = "Random forest prediction")
par(mfrow = c(1,1)) # switch back to default



#####################################################
# Apply prediction alternative with terra package
#####################################################
p_load("terra")
# see ?terra-package for an overview of the functions,
# e.g. apply arithmetic operations to rasters, project them, etc.

# Load input rasters - can be other file formats, e.g. geotiff
plancurv <- rast('plancurv.asc')
profcurv <- rast('profcurv.asc')

plancurv[ plancurv < -0.2 ] <- -0.2
plancurv[ plancurv >  0.2 ] <-  0.2
profcurv[ profcurv < -0.04 ] <- -0.04
profcurv[ profcurv >  0.04 ] <-  0.04


# Don't forget the transformations we applied to the data
log.carea <- log10(rast('carea.asc'))
 # or
 # carea <- raster('carea.asc')
 # log.carea <- calc(carea, fun= function(x){log10(x)})
names(log.carea) <- "log.carea"

slope <- rast('slope.asc') * 180 / pi
cslope <- rast('cslope.asc') * 180 / pi

distroad <- rast('distroad.asc')
distroad[distroad > 300] <- 300
distdeforest <- rast('distdeforest.asc')
distdeforest[distdeforest > 300] <- 300


# Stack grids
layers <- c(slope, log.carea, plancurv, profcurv, cslope, distroad, distdeforest)
nlyr(layers)
layers
plot(layers)

# Random forest example
predfun.rf <- function(object, newdata) {
	predict(object, newdata = newdata, type = "prob")[,2]
}

rf <- terra::predict(layers, fit, fun = predfun.rf)
terra::writeRaster(rf, "pred_rf.tif", overwrite = TRUE)
# Note that you cannot save SpatRaster objects from terra
# using R's saveRDS function!

plot(rf)

# nicer version of map:
nice_map("pred_rf.tif", pts = d[d$slides89 == "TRUE",], 
         title = "Random forest prediction")

# This the end of an overview of the use of statistical learning techniques
# for spatial modeling of a binary response variable!
