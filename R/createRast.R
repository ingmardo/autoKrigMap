#' @title Create a RasterStack object of specified resolution and extent of multiple regionalised digital soil mapping data.
#'
#' @description Create a RasterStack object of specified resolution based on a boundary and selected RasterLayer objects.
#'
#' @param l path Character. Vector of full path names; the default corresponds to the working directory, getwd().
#' @param y  Extent object, or any object (e.g. sp, sf objects) from which an Extent object can be extracted.
#' @param res vector of two positive integers. To, optionally, set the output resolution of cells in the RasterStack object in meter (unit). The default setting is c(2, 2).
#' @param width numeric; buffer distance to adapt the coverage of the RasterStack according to the extend object (y).  The default value is \code{0}.
#'
#' @return SpatRaster object
#'
#' @examples
#' createRasterStack(path, y = field, res = c(10, 10), width = 0)
#'
#' @import terra
#' @rawNamespace
#' @export


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
