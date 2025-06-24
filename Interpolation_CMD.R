# -------------------------------------------------------------------------------
# Content: Automatic interpolation of sensor point data in the command line (Windows)   
# Date: 2025-06-24
# Author: Ingmar Schroeter
# E-mail: ingmar.schroeter@hnee.de
# Version 1.0
# -------------------------------------------------------------------------------

# Load packages ----
library(sf)
library(gstat)
library(terra)
library(dplyr)
library(concaveman)
library(cowplot)
library(ggplot2)
library(optparse)

# Helper functions ----

## creatRast ----
createRast <- function(boundary, res = c(10, 10), buffer = 10, value = 0) {
  if (class(boundary)[1] == "sf") {
    boundary <- vect(boundary)
  }
  ext <- as.numeric(terra::ext(boundary)[1:4])
  ext_xmin <- round(ext[1] - buffer)
  ext_ymin <- round(ext[3] - buffer)
  ext_xmax <- round(ext[2] + buffer)
  ext_ymax <- round(ext[4] + buffer)
  xcol <- (ext_xmax - ext_xmin)/res[1]
  yrow <- (ext_ymax - ext_ymin)/res[2]
  ext_xmax <- ext_xmin + ceiling(xcol) * res[1]
  ext_ymax <- ext_ymin + (ceiling(yrow) * res[2])
  
  ext_r <- terra::ext(ext_xmin, ext_xmax, ext_ymin, ext_ymax)
  
  base_rast <- rast(ext_r, ceiling(yrow), ceiling(xcol))
  terra::values(base_rast) <- NA
  terra::crs(base_rast) <- terra::crs(boundary)
  boundary_buf <- terra::buffer(boundary, width = buffer)
  terra::values(base_rast) <- value
  base_rast <- terra::mask(x = base_rast, mask = boundary_buf)
  return(base_rast)
}

## autofitVariogram3 ----
autofitVariogram3 <- function(formula, 
                              input_data, 
                              model = c("Sph", "Exp", "Gau", "Ste"), 
                              kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), 
                              fix.values = c(NA, NA, NA), 
                              verbose = FALSE, 
                              GLS.model = NA, 
                              start_vals = c(NA, NA, NA), 
                              miscFitOptions = list(),
                              fit.method = 7, ...)
{
  # Check for anisotropy parameters
  if ("alpha" %in% names(list(...))) 
    warning("Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.")
  
  # Take the misc fit options and overwrite the defaults by the user specified ones
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, 
                                min.np.bin = 5,
                                orig.behavior = TRUE,
                                num.bins = NA, 
                                equal.width.bins = FALSE,
                                equal.np.bins = FALSE)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
  
  if (is(input_data, "Spatial")) {
    longlat = !is.projected(input_data)
    if (is.na(longlat)) 
      longlat = FALSE
    diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1, 2]
  }
  else {
    longlat = st_is_longlat(input_data)
    if (is.na(longlat)) 
      longlat = FALSE
    bb = sf::st_bbox(input_data)
    diagonal = sqrt(((bb$xmax - bb$xmin)^2) + ((bb$ymax - bb$ymin)^2))
  }
  
  ### BEGIN MODIFICATIONS ###
  
  if(miscFitOptions[["orig.behavior"]]){
    if(verbose) cat ("Boundaries as defined by original autofitVariogram...\n\n")
    # compute boundaries the old way
    boundaries = c(2,4,6,9,12,15,25,35,50,65,80,100) * diagonal/100         # Boundaries for the bins in km
    # If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
    # request by Jon Skoien
    if(miscFitOptions[["merge.small.bins"]]) {
      # bin the old way
      if(verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
      while(TRUE) {
        if(length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
        boundaries = boundaries[2:length(boundaries)]
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
        } else {
          experimental_variogram = variogram(g, boundaries = boundaries, ...)
        }
      }
    }
    
    ### equal-width bins (approximately, variogram does its own binning too) ###
  } else if(miscFitOptions[["equal.width.bins"]]){
    if(verbose) cat("Using equal-width bins...\n")
    if('width' %in% names(list(...))) stop('Cannot pass width when equal.width.bins = TRUE. Supply "init.width" in
                                             miscFitOptions instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # user must supply either width or num.bins
    if(!'init.width' %in% names(miscFitOptions)){
      if(is.na(miscFitOptions[['num.bins']])) stop('when equal.width.bins = TRUE, user must supply either init.width or num.bins as well.')
      width <- diagonal/miscFitOptions[['num.bins']]
      if(verbose) cat("initial width not provided. Calculating using num.bins.\n")
    } else {
      width <- miscFitOptions[['init.width']]
      if(verbose) cat("initial width provided.\n")
    }
    # get the experimental variogram
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, width = width, ...)
    } else {
      experimental_variogram = variogram(g, width = width, ...)
    }
    # merge small bins if requested
    if(miscFitOptions[['merge.small.bins']]){
      if(verbose) cat("Checking if any bins have less than ", miscFitOptions[["min.np.bin"]], " points, merging bins when necessary...\n")
      iter <- 0
      maxiter <- 1000
      while(TRUE){
        if(!any(experimental_variogram$np < miscFitOptions[["min.np.bin"]])) break
        # increase width by 10% and try again
        width <- width*1.1
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, width = width, ...)
        } else {
          experimental_variogram = variogram(g, width = width, ...)
        }
        iter <- iter + 1
        if(iter > maxiter){
          cat('maximum number of interations reached. Try decreasing min.np.bin or init.width.\n\n')
          break
        }
      }
    }
    
    ### equal observation count bins ###
  } else if(miscFitOptions[["equal.np.bins"]]){
    if(verbose) cat("Using bins of equal observation counts...\n")
    if('boundaries' %in% names(list(...))) stop('Cannot pass boundaries when equal.np.bins is TRUE. Pass num.bins or min.np.bin instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # get a sorted list of distances
    dists <- sort(spDists(input_data))
    
    # apply the cutoff
    dists <- dists[dists < diagonal & dists > 0]
    
    # split the data into bins based on number of observations
    if(is.na(miscFitOptions[['num.bins']])){
      
      # compute number of bins based on the minimum number of observations per bin
      miscFitOptions[['num.bins']] <- floor(0.5*length(dists)/miscFitOptions[['min.np.bin']])
      if(verbose) cat("num.bins not supplied. Setting num.bins =", miscFitOptions[['num.bins']], 'based on min.np.bin.\n')
    }
    
    cat("checking bins, decreasing num.bins if necessary... \n")
    while(TRUE){
      
      # compute interval based on the number of bins
      interval <- length(dists)/miscFitOptions[['num.bins']]
      # define boundaries
      boundaries <- rep(NA, miscFitOptions[['num.bins']])
      for(i in 1:miscFitOptions[['num.bins']]){
        boundaries[i] <- dists[round(i*interval)]
      }
      if(length(boundaries == length(unique(boundaries)))) break
      # reduce number of bins
      miscFitOptions[['num.bins']] <- miscFitOptions[['num.bins']] - 1
    }
    
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    } 
  } else {
    # default behavior of variogram
    cat("No binning action specified in miscFitOptions.\n\n")
    
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, ...)
    }
  }
  
  ### END MODIFICATIONS ###  
  
  # set initial values
  if (is.na(start_vals[1])) { 
    initial_nugget = min(experimental_variogram$gamma)
  } else {
    initial_nugget = start_vals[1]
  }
  if (is.na(start_vals[2])) { # Range
    initial_range = 0.1 * diagonal  # 0.10 times the length of the central axis through the area
  } else {
    initial_range = start_vals[2]
  }
  if (is.na(start_vals[3])) { # Sill
    initial_sill = mean(c(max(experimental_variogram$gamma), 
                          median(experimental_variogram$gamma)))
  } else {
    initial_sill = start_vals[3]
  }
  
  # Determine what should be automatically fitted and what should be fixed
  # Nugget
  if (!is.na(fix.values[1])) {
    fit_nugget = FALSE
    initial_nugget = fix.values[1]
  } else 
    fit_nugget = TRUE
  
  # Range
  if (!is.na(fix.values[2])) {
    fit_range = FALSE
    initial_range = fix.values[2]
  } else 
    fit_range = TRUE
  
  # Partial sill
  if (!is.na(fix.values[3])) {
    fit_sill = FALSE
    initial_sill = fix.values[3]
  } else 
    fit_sill = TRUE
  getModel = function(psill, 
                      model, 
                      range, 
                      kappa, 
                      nugget, 
                      fit_range, 
                      fit_sill, 
                      fit_nugget, 
                      fit.method, 
                      verbose) {
    if (verbose) 
      debug.level = 1
    else debug.level = 0
    if (model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if (is.na(start_vals[1])) 
        nugget = 0
      if (is.na(start_vals[2])) 
        range = 1   # If a power mode, range == 1 is a better start value
      if (is.na(start_vals[3])) 
        sill = 1
    }
    if(fit.method==5){
      #fit.variogram.reml(formula, locations, data, model, debug.level = 1, set, degree = 0)
      #	obj = try(fit.variogram.reml(experimental_variogram,
      #					        model=vgm(psill=psill, model=model, range=range,
      #							nugget=nugget,kappa = kappa),
      #			  fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
      #			  debug.level = 0, fit.method=fit.method), TRUE)
    }
    else if(fit.method==8){
      #fit.variogram.gls(formula, data, model, maxiter = 30,
      #eps = .01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf,
      #plot = FALSE
      #obj = try()
    } else {
      obj = try(fit.variogram(experimental_variogram,
                              model=vgm(psill=psill, 
                                        model=model, 
                                        range=range,
                                        nugget=nugget,
                                        kappa = kappa),
                              fit.ranges = c(fit_range), 
                              fit.sills = c(fit_nugget, fit_sill),
                              debug.level = 0, 
                              fit.method = fit.method), TRUE)
    }
    if("try-error" %in% class(obj)) {
      #print(traceback())
      warning("An error has occured during variogram fitting. Used:\n",
              "\tnugget:\t", nugget,
              "\n\tmodel:\t", model,
              "\n\tpsill:\t", psill,
              "\n\trange:\t", range,
              "\n\tkappa:\t",ifelse(kappa == 0, NA, kappa),
              "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj)
      return(NULL)
    } else return(obj)
  }
  
  
  # Automatically testing different models, the one with the smallest sums-of-squares is chosen
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1
  for (m in test_models) {
    if (m != "Mat" && m != "Ste") {  # If not Matern and not Stein
      model_fit = getModel(initial_sill - initial_nugget, 
                           m, initial_range, 
                           kappa = 0, 
                           initial_nugget, 
                           fit_range, 
                           fit_sill, 
                           fit_nugget, 
                           verbose = verbose, 
                           fit.method=fit.method)
      if(!is.null(model_fit)) {       # skip models that failed
        vgm_list[[counter]] = model_fit
        SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
      counter = counter + 1
    } else {                 # Else loop also over kappa values
      for(k in kappa) {
        model_fit = getModel(initial_sill - initial_nugget, 
                             m, 
                             initial_range, 
                             k, 
                             initial_nugget, 
                             fit_range, 
                             fit_sill, 
                             fit_nugget, 
                             verbose = verbose, 
                             fit.method = fit.method)
        if (!is.null(model_fit)) {
          vgm_list[[counter]] = model_fit
          SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))
        }
        counter = counter + 1
      }
    }
  }
  
  # Check for negative values in sill or range coming from fit.variogram
  # and NULL values in vgm_list, and remove those with a warning
  strange_entries = sapply(vgm_list, function(v) any(c(v$psill, 
                                                       v$range) < 0) | is.null(v))
  if (any(strange_entries)) {
    if (verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }
  if (verbose) {
    cat("Selected:\n")
    print(vgm_list[[which.min(SSerr_list)]])
    cat("\nTested models, best first:\n")
    tested = data.frame(`Tested models` = sapply(vgm_list, 
                                                 function(x) as.character(x[2, 1])), kappa = sapply(vgm_list,
                                                                                                    function(x) as.character(x[2, 4])), SSerror = SSerr_list)
    tested = tested[order(tested$SSerror), ]
    print(tested)
  }
  result = list(exp_var = experimental_variogram, 
                var_model = vgm_list[[which.min(SSerr_list)]], 
                sserr = min(SSerr_list))
  class(result) = c("autofitVariogram", "list")
  return(result)
  
}


## Interpolation ----
#' @title Automatically perform ordinary kriging
#'
#' @description This function performs automatic kriging on the given dataset. The variogram is generated automatically using a modified version of the \code{autofitVariogram} function of the \code{automap} package.
#' @param x An sf object with geometry type points or a SpatialPointsDataFrame class containing all the data points to be interpolated.
#' @param boundary object of class \code{sf} or \code{sp} that defines the boundary of a field. The projection (CRS) should be in UTM coordinates.
#' @param res vector of two positive integers. To, optionally, set the output resolution of the cells in each direction (horizontally and vertically) in meter (unit). The default setting is c(2, 2).
#' @param block vector of two positive integers. Use this parameter to pass on a specification for the block size. Default is c(10,10).
#' @param mod List of models that will be tested during automatic variogram fitting.
#' @param maxdist for local kriging: only observations within a distance of maxdist from the prediction location are used for prediction or simulation; if combined with nmax, both criteria apply
#' @param nmax for local kriging: the number of nearest observations that should be used for a kriging prediction or
#' simulation, where nearest is defined in terms of the space of the spatial locations. By default, \code{nmax = 40}
#' @param cutoff spatial separation distance up to which point pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three.
#' @param width the width of subsequent distance intervals into which data point pairs are grouped for semivariance estimates
#'
#' @return An object which is a list with components:
#' @return \item{autoKrigMaps}{SpatRaster of interpolated data}
#' @return \item{models}{Results of fitted variogram models}
#' @return \item{formula}{Formula used for ordinary and simple kriging e.g. z~1}
#' @return \item{variogramPlots}{Plot containing all fitted variogram models of the given dataset.}
#'
#' @examples
#' autoKrigMap(x = points_sf, boundary = field_sf)
#'
#' @import concaveman 
#' @import cowplot 
#' @import terra
#' @import sp
#' @import sf
#' @import dplyr
#' @import gstat
#'
#' @export

autoKrigMap <- function(x, 
                        boundary = NULL, 
                        res = c(2,2) , 
                        block = c(10,10), 
                        maxdist = Inf, 
                        nmax = 40, 
                        cutoff = NULL, 
                        mod = c("Sph", "Exp", "Gau"), 
                        width = 10,
                        epsg_code = 25833){
  
  if(class(x)[1] == "SpatialPointsDataFrame"){
    x <- st_as_sf(x)
  }
  
  if (is.na(sf::st_crs(x))) {
    message("Error: CRS is missing. Please assign a coordinate reference system.")
    stop("Missing CRS detected.")
  }
  
  if (is.na(sf::st_crs(boundary))) {
    message("Error: CRS for boundary is missing. Please assign a coordinate reference system.")
    stop("Missing CRS detected for boundary.")
  }
  
  if (sf::st_crs(x)$epsg == 4326){
    message(paste0("Spatial data were transformed from Geodetic coordinate system (EPSG:4326) into UTM coordinates ",
                   "(EPSG:", epsg_code,")."))
    
    x <- st_transform(x = x, crs = epsg_code)
  } 
  
  # Remove columns with non-numeric values
  x <- dplyr::select(x, where(is.numeric))
  
  if (length(boundary) > 0){
    if (is.na(sf::st_crs(boundary))) {
      message("Error: CRS is missing. Please assign a coordinate reference system.")
      stop("Missing CRS detected for boundary.")
    }
  }
  
  # remove duplicate locations ----
  data_sf <- x[!duplicated(st_geometry(x)), ]
  
  
  # check if boundary was set 
  if (length(boundary) == 0){
    
    # more complex hull curve 
    boundary <- concaveman::concaveman(data_sf, concavity = 2, length_threshold = 100)
  }
  
  # Create kriging surface ----
  basis.r <- createRast(boundary = boundary, res = c(res[1], res[2]), buffer = 50, value = -9999)
  base_df <- as.data.frame(basis.r, xy=T)
  basis.grid <- st_as_sf(x = base_df, coords = c("x", "y"), crs = epsg_code)
  
  # Create list object for results ----
  formula <- list()
  model <- list()
  ok_st <- list()
  
  for (i in seq_len(ncol(data_sf) - 1)) {
    message("Processing variable: ", names(data_sf)[i])
    #i <- 2
    time <- Sys.time()
    formula[[i]] <- paste0(names(data_sf)[i], " ~ ", 1)
    outcome <- names(data_sf)[i]
    
    # remove NA values 
    clean_sf <- na.omit(data_sf[,i])     
    
    if (length(cutoff) == 1){
      model[[i]] <- autofitVariogram3(formula = as.formula(paste(outcome, "~ 1")),
                                      input_data = clean_sf,
                                      model = mod,
                                      cutoff = cutoff,
                                      miscFitOptions = list(equal.width.bins = T, orig.behavior = F , init.width = width))
    } else {
      model[[i]] <- autofitVariogram3(formula = as.formula(paste(outcome, "~ 1")), 
                                      input_data = clean_sf, 
                                      model = mod)
    }
    
    ok.results <- gstat::krige(formula = as.formula(paste(outcome, "~ 1")),
                               locations = clean_sf,
                               newdata = basis.grid,          # basis.grid 
                               model = model[[i]]$var_model,
                               maxdist = maxdist,
                               nmax= nmax,
                               block = block)
    
    
    ok_df <- data.frame(st_coordinates(ok.results), st_drop_geometry(ok.results))
    ok_r <- rast(ok_df, type="xyz")
    crs(ok_r) <- paste0("epsg:", st_crs(boundary)$epsg)
    ok_st[[i]] <- ok_r$var1.pred
    names(ok_st[[i]]) <- names(x)[i]
    print(Sys.time()-time)
  }
  
  
  r_st <- rast(ok_st)
  
  # Check if CRS are the same  
  if (crs(r_st, describe = TRUE)$code != st_crs(boundary)$epsg) 
    stop(paste("Projections of r_st and boundary do not match:\n"))
  #crs(r_st)
  
  # clip raster stack to mask
  r.mask <- terra::mask(x = r_st, mask = boundary)
  
  # create plot of all variograms
  var.list <- list()
  
  # Create variogram plots with ggplot 
  for (i in seq_along(model)) {
    # Extract data
    exp_var <- model[[i]]$exp_var
    fit_model <- model[[i]]$var_model
    
    # Extract parameters
    model_name <- fit_model[2, "model"]
    nugget     <- fit_model[1, "psill"]
    psill      <- fit_model[2, "psill"]
    range_val  <- fit_model[2, "range"]
    
    # Build model line: sequence of distances and corresponding gamma values
    dist_seq <- seq(0, max(exp_var$dist), length.out = 200)
    gamma_vals <- gstat::variogramLine(fit_model, dist_vector = dist_seq)$gamma
    model_line <- data.frame(dist = dist_seq, gamma = gamma_vals)
    preds <- variogramLine(fit_model, maxdist = max(exp_var$dist))
    
    
    # Create ggplot
    p <- 
      ggplot() +
      geom_point(data = exp_var, aes(x = dist, y = gamma), 
                 color = "dodgerblue", 
                 size = 2.5) +
      ylim(0, max(exp_var$gamma)) +
      geom_line(data = preds, aes(x = dist, y = gamma), 
                linewidth = 1, 
                alpha = 0.8)+
      labs(
        title = names(r.mask)[i],
        x = "Distance",
        y = "Semi-variance"
      ) +
      annotate("text", 
               x = Inf, 
               y = -Inf, 
               hjust = 1.1, 
               vjust = -0.5, 
               size = 3.5, 
               fontface = "bold",
               label = paste0("Model: ", model_name, "\n",
                              "Nugget: ", round(nugget, 2), "\n",
                              "PSill: ", round(psill, 2), "\n",
                              "Range: ", round(range_val, 2))) +
      theme_bw()
    
    var.list[[i]] <- p
  }
  
  variogramPlots <- cowplot::plot_grid(plotlist = var.list)
  
  # Export results ----
  crs(r.mask)  <- paste0("epsg:", epsg_code)
  
  results <- list()
  results$maps <- r.mask
  results$models <- model
  results$formulas <- formula
  results$varioPlots <- variogramPlots
  return(results)
}


# Parameter transfer in the CMD ----
args <- commandArgs(trailingOnly=TRUE)
farm_name  <- args[1] # string for farm name (e.g. "LPP") 
field_name <- args[2] # string for field name (e.g. "1392")
boundary   <- args[3] # path to field boundary Shapefile in WGS84
bd         <- args[4] # Integer of [1] boundary exist (default) or [0] boundary is auto created    
points_sf  <- args[5] # path to sensor point data in WGS84 (EPSG: 4326) or UTM (EPSG: 25833) 
date       <- as.character(args[4]) # Datum der Sensorkartierung

# Examples 
#farm_name  <- "LPP"
#field_name <- "1392"
#boundary   <- "data/1392/PP_1392_Schlaggrenze_Neu.shp" 
#bd         <-  1   
#points_sf  <- "data/1392/LPP_1392_EPSG_4326_Geophilus_2017-09-27.shp"
#date       <- "2017-09-27"   


#  Boundary festlegen ----

if (bd == 1){
  # Einladen der Schlaggrenze ----
  sg_sf <- sf::st_read(boundary)
} else {
  sg_sf <- NULL
}

# Sensordaten einladen ----
sdat_sf <- sf::st_read(points_sf)

# Interpolation starten ----
results_akm <- autoKrigMap(x = sdat_sf[2], 
                       boundary = sg_sf, 
                       nmax = 40)


# Save as GeoTiff ----

# Set output path
dir.create("output")
path.out <- "output"
r_stack <- results_akm$maps 

for (i in seq_len(nlyr(r_stack))) {
  terra::writeRaster(r_stack[[i]], paste0(path.out,"/", farm_name, "_", field_name, "_OK_",names(r_stack)[i],"_",date,".tif"), overwrite=TRUE)
}


sessionInfo()

# ---
# End ----
# ---

