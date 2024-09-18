####################################################
# Landslide Susceptibility Modeling with R
# Workshop at the conference of the
# Romanian Association of Geomorphologists
# Iasi, 19 September 2024
# (c) Alexander Brenning - University of Jena
####################################################
# Install required R packages
####################################################

# install.packages('pacman')
library("pacman")
p_load("terra")
p_load("sf")
p_load("RSAGA")
p_load("ROCR")
p_load("mgcv")
p_load("rpart")
p_load("randomForest")
p_load("sperrorest")
p_load("ggplot2")
p_load("patchwork")
p_load("iml")
