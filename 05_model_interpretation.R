####################################################
# Landslide Susceptibility Modeling with R
# Workshop at the conference of the
# Romanian Association of Geomorphologists
# Iasi, 18 September 2024
# (c) Alexander Brenning - University of Jena
#####################################################
# Model interpretation
# Case study: landslides in Ecuador
#####################################################

setwd("data")

# Fewer decimal places - works better for instructor:
options(digits=4)

library("pacman")
p_load("sperrorest")    # PVI in spatial CV
p_load("randomForest")  # Random forest
p_load("rpart")         # CART
p_load("ggplot2")       # Plotting
p_load("patchwork")     # Plotting
p_load("iml")           # ML model interpretation

# Load the data set:
d <- readRDS("landslides_transformed.rds")

fo <- slides89 ~ slope + plancurv + profcurv + log.carea + 
  cslope + distroad + distdeforest


############################################################
### Random Forest:
### start with tools from the randomForest package
############################################################

fit <- randomForest(fo, data = d, importance = TRUE)

# Random forest package comes with two built-in variable importance criteria,
# one of which (mean decrease in Gini purity) is model-specific:
par(mfrow = c(1,2))
randomForest::varImpPlot(fit, type = 1, pch = 19)
randomForest::varImpPlot(fit, type = 2, pch = 19)

# Examples of a PDP using function from randomForest package:
par(mfrow = c(2,2))
randomForest::partialPlot(fit, x.var = "distroad", pred.data = d)
randomForest::partialPlot(fit, x.var = "distdeforest", pred.data = d)
randomForest::partialPlot(fit, x.var = "cslope", pred.data = d)
randomForest::partialPlot(fit, x.var = "log.carea", pred.data = d)

rf_imp <- randomForest::importance(fit)
top3 <- rownames(rf_imp[order(rf_imp[,"MeanDecreaseAccuracy"], decreasing = TRUE),])[1:3]
randomForest::partialPlot(fit, x.var = top3[1], pred.data = d, xlab = top3[1])
randomForest::partialPlot(fit, x.var = top3[2], pred.data = d, xlab = top3[2])
randomForest::partialPlot(fit, x.var = top3[3], pred.data = d, xlab = top3[3])



############################################################
### Random Forest:
### Use sperrorest for (spatial) variable importance
############################################################

# Calculate importance for these variables:
imp_vars <- all.vars(fo)[-1]
imp_vars

# Perform spatial cross-validation; using simplified settings to test the code:
res <- sperrorest(formula = fo, data = d, coords = c("x","y"),
                  model_fun = randomForest, 
                  model_args = list(ntree = 100, nodesize = 100), 
                  pred_fun = mypred_rf,
                  smp_fun = partition_kmeans, 
                  smp_args = list(repetition = 1:2, nfold = 10),
                  importance = TRUE, imp_permutations = 10,
                  imp_variables = imp_vars)
# Cross-validation estimate of error rate & accuracy:
imp <- summary(res$importance)
# ... a data.frame with results...

# E.g. mean decrease in auroc when permuting ndvi01:
imp["slope", "mean.auroc"]
# Its standard deviation over the repetitions:
imp["slope", "sd.auroc"]

imp <- imp[order(imp$mean.auroc, decreasing = TRUE),]
imp[1:5, c("mean.auroc", "sd.auroc")]

# Create a barplot - looks better with greater importances at the top:
imp <- imp[order(imp$mean.auroc, decreasing = FALSE),]
imp[1:5, c("mean.auroc", "sd.auroc")]
par(mar = c(5,7,1,1)) # bottom, left, top, right margins
barplot(imp$mean.auroc, names.arg = rownames(imp), horiz = TRUE, las = 1, 
        xlab = "Mean decrease in AUROC")



############################################################
### Random Forest:
### use iml package for partial dependence and ALE plot
############################################################

p_load("iml")        # ML model interpretation
p_load("patchwork")  # Plotting

fit <- randomForest(fo, data = d, importance = TRUE)

predvars <- all.vars(fo)[-1]
predictors <- d[, predvars]
response <- d$slides89

predictor <- iml::Predictor$new(fit, data = predictors, y = response)

# Permutation variable importance on the training set (!):
imp <- iml::FeatureImp$new(predictor, loss = "ce", n.repetitions = 100, compare = "difference")
plot(imp)

# PVI, randomForest style (for comparison)
par(mfrow = c(1,1))
randomForest::varImpPlot(fit, type = 1)


# Calculate SHAP feature importance from Shapley values:
# This is very SLOW and not parallelized!
extract_shapleys <- function(x) 
  Shapley$new(predictor, 
              x.interest = predictors[x,],
              sample.size = 10 # use default of 100! -> slower
             )$results$phi
shaps <- sapply(1:nrow(d), extract_shapleys)
shap <- apply(abs(shaps), 1, mean)
# average over all four classes:
shap <- tapply(shap, rep(predvars, 2), mean)
barplot(shap[order(shap)], horiz = TRUE, las = 1, xlab = "SHAP feature importance")



# Interaction strength using iml package:
interac <- iml::Interaction$new(predictor)
plot(interac)

# log.carea and cslope and ndvi04 showed the strongest interactions with other
# predictors, therefore take a closer look at how they interact:

interac <- Interaction$new(predictor, feature = "log.carea")
plot(interac)

interac <- Interaction$new(predictor, feature = "cslope")
plot(interac)


# PDP / ALE plot using iml package (use method argument!):
effs <- iml::FeatureEffects$new(predictor, method = "pdp", features = top3)
plot(effs)

effs <- iml::FeatureEffects$new(predictor, method = "ale", features = top3, grid.size = 40)
plot(effs)

