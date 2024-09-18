####################################################
# Landslide Susceptibility Modeling with R
# Workshop at the conference of the
# Romanian Association of Geomorphologists
# Iasi, 19 September 2024
# (c) Alexander Brenning - University of Jena
####################################################
# Using R as a GIS to prepare the landslide dataset
####################################################

# This is where my raster data from Ecuador is located,
# relative to the RStudio project's folder:
setwd("data")

# install.packages('pacman')
library("pacman")
p_load("terra") # see https://r.geocompx.org/spatial-class
p_load("sf")    # see https://r.geocompx.org/spatial-class
p_load("RSAGA") # see https://cran.r-project.org/web/packages/RSAGA/vignettes/RSAGA.html


##############################################
# First steps with 'terra'
##############################################

# Digital elevation model from our study area in southern Ecuador:
dem <- terra::rast('dem.asc')
class(dem)
dem
# try res(dem), ncell(dem), ext(dem), and crs(dem)
plot(dem)
# Use packages 'tmap' and 'rasterVis' for more advanced displays.

# A mask defining the study area can be useful:
mask <- !is.na(dem)
names(mask) <- "mask"
mask
summary(mask) # NA = outside study area
plot(mask)
# If they had a different extent or resolution, we'd
# want to rearrange them using terra::resample() or
# terra::crop() functions.

# We also have landslide presence/absence data and
# a raster masking the inventoried area as raster data
# with same extent (how convenient!):
slides89 <- rast("slides89.asc")
slides89
plot(slides89) # NA = background; values 1-6 not relevant for us
slides89 <- !is.na(slides89)
slides89[!mask] <- NA

# and also 'distroad' and 'distdeforest', if you want to
# check them out...
distroad <- rast("distroad.asc")
distdeforest <- rast("distdeforest.asc")
slope <- rast("slope.asc")
cslope <- rast("cslope.asc")
carea <- rast("carea.asc")
plancurv <- rast("plancurv.asc")
profcurv <- rast("profcurv.asc")

# We can stack them all:
ecua <- c(dem, slides89, mask, 
          distroad, distdeforest, 
          slope, cslope, carea, plancurv, profcurv)
ecua
# ...but we don't always have to ;-)

# Either way, we can do arithmetic operations on them
# (local function):
dem_feet <- dem * 3.28084
plot(dem_feet)

# Focal functions (on moving windows) can also be applied
# using terra::focal() ...

# Or we can subset rasters for further analyses...
dist_slides <- as.vector(distroad[slides89])
dist_no_slides <- as.vector(distroad[!slides89])
boxplot(list("No slide" = dist_no_slides, "Slide" = dist_slides),
        ylab = "Distance to road [m]")



##############################################
# Prepare the training and test samples:
# Create random point sample and pick values
# of predictor variables from grids
##############################################

# Generate an initial random sample of grid cells:
d <- terra::spatSample(ecua, size = 10000, method = "random", 
                       na.rm = TRUE, xy = TRUE)

# Take an equal number of landslide and non-landslide grid cells:
N.sli <- sum(d$slides89)
sel <- c( sample(which(d$slides89), size = N.sli), 
          sample(which(!d$slides89), size = N.sli))
d <- d[sel,]

# R supports simple features through the 'sf' package:
# Points, lines, polygons, linestrings, multi-...
# see https://r.geocompx.org/spatial-class

# Create point geometries from x/y coordinates:
pts <- sf::st_geometry(st_as_sf(d, coords = c("x","y")))
plot(pts)

# Create point feature class with attribute data:
st_geometry(d) <- pts
d

# Save the result for subsequent analyses:
# saveRDS(d, file = "landslides.rds")


### A quick look at some maps:

# Plot map of slides89 attribute:
plot(d["slides89"], pch = 19, cex = 0.7)

# Vary point size with slope angle:
plot(d["slides89"], cex = d$slope)
# I guess a boxplot would be more suitable in this case...

# Combine raster and vector map:
plot(distroad)
cols <- c("red", "black")[(d$slides89 == "TRUE") + 1]
plot(d["slides89"], add = TRUE, col = cols, pch = 19, cex = 0.7)

# Check out the following chapter introducing 'tmap'
# and 'leaflet' for interactive map displays:
# https://r.geocompx.org/adv-map




##############################################
# Calculate terrain attributes using (R)SAGA
##############################################

# But where did we get the terrain attributes from?
# Let's couple R with SAGA GIS for advanced GIS analyses...

# Let RSAGA find SAGA GIS:
env <- RSAGA::rsaga.env()
# If it can't find it, use the path argument to specify its location.
# On my computer(s) e.g.:
# env <- rsaga.env(path = "C:/Progra~1/SAGA")

# See what you've got:
env

# Explore available libraries and modules (slow!):
libs <- rsaga.get.libraries(path = env$modules)
libs
modules <- rsaga.get.modules(libs = "ta_morphometry", env = env)
modules
# We are also going to use tools from libraries "io_grid"
# and "ta_hydrology".

# Convert the DEM ASCII grid to SAGA grid format (.sgrd):
rsaga.esri.to.sgrd("dem.asc", env = env)

# To be able to see what SAGA_CMD call was actually run, we
# can enable the display.command argument:
rsaga.esri.to.sgrd("dem.asc", env = env, display.command = TRUE)

# Or, we could do the same thing 'the long way', using
# the low-level rsaga.geoprocessor:
rsaga.get.usage("io_grid", "Import ESRI Arc/Info Grid", env = env)
rsaga.geoprocessor("io_grid", "Import ESRI Arc/Info Grid",
    param = list(GRID = "dem.sgrd", FILE = "dem.asc"), 
    display.command = TRUE, env = env)
# ...but we'd rather prefer the more convenient higher-level front-end
# function in RSAGA...


# Now let's calculate a bunch of terrain attributes
# that could be useful for modeling landslide 
# susceptibility:

# Slope, aspect, plan and profile curvature:
rsaga.slope.asp.curv("dem", out.slope = "slope", out.cplan = "plancurv",
                     out.cprof = "profcurv", method = "poly2zevenbergen", env = env)
# Create a hydrologically more meaningful DEM by filling any local sinks:
rsaga.sink.removal("dem", out.dem = "sdem", method = "fill", env=env)
# Use the multiple flow direction algorithm to calculate catchment area and slope:
rsaga.wetness.index("sdem", out.carea = "carea",
            out.cslope = "cslope", 
            area.type = "absolute", slope.type = "catchment", env=env)
# Convert terrain attribute grids into ASCII format:
rsaga.sgrd.to.esri(c("slope", "plancurv", "profcurv", "carea", "cslope"), env=env)

# Delete SAGA grids, we don't need them any more, unless you want to view them in SAGA GIS:
# shell("del *.sgrd"); shell("del *.mgrd"); shell("del *.sdat")



# Now you have a general overview of some of the GIS capabilities of R!
# The book 'Geocomputation with R' is your go-to place
# if you want to learn more about using R as a GIS!
# https://r.geocompx.org/
