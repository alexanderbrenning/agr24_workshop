####################################################
# Landslide Susceptibility Modeling with R
# Workshop at the conference of the
# Romanian Association of Geomorphologists
# Iasi, 19 September 2024
# (c) Alexander Brenning - University of Jena
####################################################
# Hyperparameter tuning
# Case study: landslides in Ecuador
#####################################################

library("pacman")
p_load("mgcv")
p_load("rpart")
p_load("randomForest")
p_load("sperrorest")

# Fewer decimal places - works better for instructor:
options(digits=4)

setwd("data")

# Load the data set:
d <- readRDS("landslides_transformed.rds")

# Make sure we use the exact same formula in all models:
fo <- slides89 ~ slope + plancurv + profcurv + log.carea + 
  cslope + distroad + distdeforest



#####################################################
# Start by exploring sensitivity of classification
# tree to maxdepth hyperparameter;
# turn built-in tuning off: cp = 0
#####################################################

p_load("rpart")

# Classification trees with different maxdepth, no pruning:
fit <- rpart(fo, data = d, control = rpart.control(cp=0, maxdepth=3))
par(xpd = TRUE, mfrow = c(1,1))
plot(fit)
text(fit, use.n = TRUE)

fit <- rpart(fo, data = d, control = rpart.control(cp = 0, maxdepth = 10))
plot(fit)
text(fit, use.n = TRUE)

# What strategy works best?
# a) "Guess" the "right" maxdepth (or cp) value?
# b) Optimize it as part of the model building process?

# Until now we have only learned how to perform cross-validation
# on models with fixed hyperparameters:

# Control parameters for rpart:
ctrl <- rpart.control(cp = 0, maxdepth = 3)

# A wrapper function that extracts the predicted probabilities from rpart predictions:
mypred_rpart <- function(object, newdata) {
  predict(object, newdata, type = "prob")[,2]
}

# Perform 5-repeated 10-fold spatial cross-validation:
res <- sperrorest(formula = fo, data = d, 
                  coords = c("x","y"),
                  model_fun = rpart, 
                  model_args = list(control = ctrl),
                  pred_fun = mypred_rpart,
                  smp_fun = partition_kmeans, 
                  smp_args = list(repetition = 1:5, nfold = 10))
# Cross-validation estimate of AUROC:
summary(res$error_rep)["test_auroc",1]

# We can use this to explore how hyperparameter values influence
# cross-validation estimates of the accuracy Let's do that, as a first step.

# Let's put this inside a wrapper function:

estimate_auroc <- function(maxdepth) {
  # Control parameters for rpart:
  ctrl <- rpart.control(cp = 0, maxdepth = maxdepth)
  # Perform 1-repeated 5-fold spatial cross-validation:
  res <- sperrorest(formula = fo, data = d, coords = c("x","y"),
                    model_fun = rpart, 
                    model_args = list(control = ctrl), 
                    pred_fun = mypred_rpart,
                    smp_fun = partition_kmeans, 
                    # simplified settings:
                    smp_args = list(repetition = 1:5, nfold = 10),
                    progress = FALSE)
  # Cross-validation estimate of AUROC:
  auroc <- summary(res$error_rep)["test_auroc",1]
  return(auroc)
}
# Note that this is a bit 'dirty' because this function accesses global
# variables such as 'd' and 'fo'.

# Let's give it a try:
estimate_auroc(7)

maxdepths <- c(1,2,3,5,7,10,15,20)
aurocs <- sapply(maxdepths, estimate_auroc)
plot(aurocs ~ maxdepths, type="l")
points(aurocs ~ maxdepths, pch="+")
cbind(maxdepths, aurocs)

# What have we learned?



#####################################################
# *Properly* tuning hyperparameters for
# predictive modelling
#####################################################

selftuning_rpart <- function(formula, data, 
                             maxdepths = c(1,2,3,5,7,10,15,20))
{
  aurocs <- rep(NA, length(maxdepths))
  for (i in 1:length(maxdepths)) {
    ctrl <- rpart.control(cp = 0, maxdepth = maxdepths[i])
    # Note that this sperrorest call MUST use the 'formula' and 'data'
    # objects that were passed to this function - not the 'fo' and 'd'
    # objects!
    res <- sperrorest(formula = formula, data = data, coords = c("x","y"),
                      model_fun = rpart, 
                      model_args = list(control = ctrl), 
                      pred_fun = mypred_rpart, 
                      smp_fun = partition_kmeans,
                      # using simplified settings for speed-up:
                      smp_args = list(repetition = 1:5, nfold = 10),
                      progress = FALSE)
    aurocs[i] <- summary(res$error_rep)["test_auroc",1]
  }
  
  opt.maxdepth <- maxdepths[ which.max(aurocs) ]
  cat("Optimal tree depth:", opt.maxdepth, "\n")
  # Possible improvement: Use smallest maxdepth value
  # whose performance is nearly as good as the optimal
  # one. -> Occam's razor!

  # Note that this rpart call MUST use the 'formula' and 'data'
  # objects that were passed to this function - not the 
  # global 'fo' and 'd' objects!
  fit <- rpart(formula, data = data, 
               control = rpart.control(cp = 0, maxdepth = opt.maxdepth))
  return(fit)
}

fit <- selftuning_rpart(fo, d)
plot(fit)
text(fit, use.n = TRUE)

# How accurate is this particular fitted model, fit?
# Can we even tell?



# Perform spatial cross-validation on selftuning_rpart:
# future::plan(future.callr::callr, workers = 10)
res <- sperrorest(formula = fo, data = d, coords = c("x","y"),
                  model_fun = selftuning_rpart, 
                  pred_fun = mypred_rpart, 
                  # mode_fold = "future", mode_rep = "sequential",
                  smp_fun = partition_kmeans, 
                  # simplified settings, not patient enough:
                  smp_args = list(repetition = 1:1, nfold = 10))
# Note that the estimated optimal tree depths are not all the same:
# there is some random variability.
# Cross-validation estimate of overall accuracy:
summary(res$error_rep)["test_auroc",1]


# Next steps:
# - Adapt this to random forest
# - Apply to more than one hyperparameter using a matrix instead of 
#   an array of candidate hyperparameter values (two-dim. grid search).




#####################################################
# Now re-do this for the random forest model,
# focusing on the nodesize hyperparameter, which
# controls the complexity of the trees.
#####################################################

p_load("randomForest")

# Random forest models with different nodesize hyperparameter settings,
# all other parameters at default values:
fit <- randomForest(fo, data = d, nodesize = 1) # grow trees until nodes are pure
fit

fit <- randomForest(fo, data = d, nodesize = 5) # default for classification
fit

fit <- randomForest(fo, data = d, nodesize = 50) # might underfit the data?
fit

# A wrapper function that extracts the predicted 
# probabilities from randomForest predictions:
mypred_rf <- function(object, newdata) {
  predict(object, newdata, type = "prob")[,2]
}

# Perform spatial cross-validation:
res <- sperrorest(formula = fo, data = d, 
                  coords = c("x","y"),
                  model_fun = randomForest, 
                  model_args = list(nodesize = 5),
                  pred_fun = mypred_rf,
                  smp_fun = partition_kmeans, 
                  smp_args = list(repetition = 1:5, nfold = 10))
# Cross-validation estimate of AUROC:
summary(res$error_rep)["test_auroc",1]

# Let's put this inside a wrapper function:
estimate_auroc <- function(nodesize) {
  # Perform 1-repeated 5-fold spatial cross-validation:
  res <- sperrorest(formula = fo, data = d, coords = c("x","y"),
                    model_fun = randomForest, 
                    model_args = list(nodesize = nodesize), 
                    pred_fun = mypred_rf,
                    smp_fun = partition_kmeans, 
                    # simplified settings:
                    smp_args = list(repetition = 1:5, nfold = 10),
                    progress = FALSE)
  # Cross-validation estimate of accuracy:
  auroc <- summary(res$error_rep)["test_auroc",1]
  return(auroc)
}
# Note that this is a bit 'dirty' because this function accesses global
# variables such as 'd' and 'fo'.

# Let's give it a try:
estimate_auroc(7)

nodesizes <- c(1,5,10,20,50,100,200)
aurocs <- sapply(nodesizes, estimate_auroc)
plot(aurocs ~ nodesizes, type="l")
points(aurocs ~ nodesizes, pch="+")
cbind(nodesizes, aurocs)

# What have we learned?



#####################################################
# *Properly* tuning hyperparameters for
# predictive modelling
#####################################################

selftuning_randomforest <- function(formula, data, 
                             nodesizes = c(5,10,20,50,100,200))
{
  aurocs <- rep(NA, length(nodesizes))
  for (i in 1:length(nodesizes)) {
    # Note that this sperrorest call MUST use the 'formula' and 'data'
    # objects that were passed to this function - not the 'fo' and 'd'
    # objects!
    res <- sperrorest(formula = formula, data = data, coords = c("x","y"),
                      model_fun = randomForest, 
                      model_args = list(nodesize = nodesizes[i]), 
                      pred_fun = mypred_rf, 
                      smp_fun = partition_kmeans,
                      # using simplified settings for speed-up:
                      smp_args = list(repetition = 1:5, nfold = 10),
                      progress = FALSE)
    aurocs[i] <- summary(res$error_rep)["test_auroc",1]
  }
  
  opt.nodesize <- nodesizes[ which.max(aurocs) ]
  cat("Optimal tree depth:", opt.nodesize, "\n")
  # Possible improvement: Use largest nodesize value
  # whose performance is nearly as good as the optimal
  # one. -> Occam's razor!
  
  # Note that this randomForest call MUST use the 'formula' and 'data'
  # objects that were passed to this function - not the 
  # global 'fo' and 'd' objects!
  fit <- randomForest(formula, data = data, 
               nodesize = opt.nodesize)
  return(fit)
}

fit <- selftuning_randomforest(fo, d)
plot(fit)
text(fit, use.n = TRUE)


p_load("RSAGA")
multi.local.function(
    in.grids = c("slope", "plancurv", "profcurv", "carea", "cslope", "distroad", "distdeforest"),
    out.varnames = "rfoptpred",
    fun = grid.predict, control.predict = list(type = "prob"),
    fit = fit, trafo = my.trafo, predict.column = "TRUE",
    quiet = FALSE )

# Using the nice_map function from one of the previous scripts:
nice_map("rfoptpred.asc", title = "Optimized random forest prediction")

# How accurate is this particular fitted model, fit?
# Can we even tell?


### The following step is VEEEERY slow!

# Perform spatial cross-validation on selftuning_randomforest:
# future::plan(future.callr::callr, workers = 10)
res <- sperrorest(formula = fo, data = d, coords = c("x","y"),
                  model_fun = selftuning_randomforest,
                  pred_fun = mypred_rf, 
                  # mode_fold = "future", mode_rep = "sequential",
                  smp_fun = partition_kmeans, 
                  # simplified settings, not patient enough:
                  smp_args = list(repetition = 1:1, nfold = 10))
# Cross-validation estimate of AUROC:
summary(res$error_rep)["test_auroc",1]


# Next steps:
# - Apply to more than one hyperparameter using a matrix instead of 
#   an array of candidate hyperparameter values (two-dim. grid search).
