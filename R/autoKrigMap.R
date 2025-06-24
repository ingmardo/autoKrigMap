# # Load packages ----
# library(sf)
# library(raster)
# library(gstat)
# #library(automap)
# library(concaveman)
# library(lattice)
# library(grid)
# library(cowplot)

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
#' autoKrigMap(x = msp3.sp, boundary = field.sf)
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
  
  


