# Functions for acquiring environmental data programmatically from common sources.

# get_bathymetry ----------------------------------------------------------

#' Function to download bathymetry data
#'
#' @param out_dir File path to the directory where raster will be saved. String.
#' @param region Bounding box for region as xmin, ymin, xmax, ymax.
#' @param overwrite Should existing files be overwritten. Logical.
#' @param stride Rate at which to sample underlying raster resolution. Integer.
#' @param mask_land Value of elevation that will be used to mask out land. Integer. Default is NULL.
#' @param add_slope Should a layer for slope be calculated. Logical.
#' @param neighbors Number of neighbours to use in slope calculation. Either 4 or 8.
#' @param slope_unit Unit for the slope calculation. One of "degrees" or "radians".
#' @param add_coastdist Should a layer for distance to coastline be calculated. Logical.
#' @param dist_unit Unit for the distances One of "m" or "km".
#'
#' @details
#'
#' This function uses rerddap to download the GEBCO_2020 Grid, which is a
#' continuous terrain model for oceans and land at 15 arc-second intervals (0.004 deg).
#'
#' Stride is passed to rerddap if you want to downsample the resolution of the
#' raster data before downloading. For example, a stride of 2 will only acquire
#' every second grid cell and a stride of 10 will acquire every 10th grid cell.
#'
#' If mask_land is not null, all elevation values greater than mask_land will be converted to NA.
#'
#' @return A netcdf file named 'bathymetry.nc' saved to out_dir. If add_slope or
#' add_coast are set to TRUE, then additional files named 'slope.nc' and 'coastdist.nc'
#' will be saved. Test git
#' @export
#' @seealso [get_slope()], [get_coastdist()]
#'
#' @examples
#'# Just get bathymetry
#'get_bathymetry(out_dir = 'tmp',
#'               region = c(-68, 44, -63, 46),
#'               overwrite = FALSE,
#'               stride = 1,
#'               mask_land = 0,
#'               add_slope = FALSE,
#'               neighbors = 8,
#'               slope_unit = c('degrees','radians'),
#'               add_coastdist = FALSE,
#'               dist_unit = c('m','km'))
#'r <- terra::rast('tmp/bathymetry.nc')
#'terra::plot(r)
#'
#'# Get bathymetry and calculate slope and distance from coast
#' get_bathymetry(out_dir = 'tmp',
#'               region = c(-68, 44, -63, 46),
#'               overwrite = FALSE,
#'               stride = 1,
#'               mask_land = 0,
#'               add_slope = TRUE,
#'               neighbors = 8,
#'               slope_unit = 'degrees',
#'               add_coastdist = TRUE,
#'               dist_unit = 'km')
#' r <- terra::rast('tmp/slope.nc')
#' terra::plot(r, main = terra::longnames(r))
#' r <- terra::rast('tmp/coastdist.nc')
#' terra::plot(r, main = terra::longnames(r))
#' unlink('tmp', recursive = TRUE)
get_bathymetry <- function(out_dir = NULL,
                           region,
                           overwrite = FALSE,
                           stride = 1,
                           mask_land = NULL,
                           add_slope = FALSE,
                           neighbors = 8,
                           slope_unit = c('degrees','radians'),
                           add_coastdist = FALSE,
                           dist_unit = c('m','km')
) {

  bath_file <- paste0(out_dir,'/bathymetry.nc')

  if (file.exists(bath_file) == FALSE | overwrite == TRUE) {
    b <- rerddap::griddap(rerddap::info('GEBCO_2020'),
                          latitude = region[c(2,4)],
                          longitude = region[c(1,3)],
                          stride = stride,
                          store = rerddap::disk(path = out_dir,
                                                overwrite = TRUE)
    )
    file.rename(b$summary$filename, paste0(dirname(b$summary$filename), '/bathymetry.nc'))
    message(paste('Saving:', bath_file))
  } else message(paste('Skipping:', bath_file, 'exists in destination'))

  if (is.null(mask_land) == FALSE) {
    if (is.numeric(mask_land) == FALSE | length(mask_land) > 1) stop('mask_land must be a single numeric value')
    r <- terra::rast(bath_file)
    r <- terra::ifel(r > mask_land, NA, r)
    terra::writeCDF(r, bath_file,
                    longname = terra::longnames(r),
                    unit = terra::units(r),
                    varname = terra::varnames(r),
                    overwrite = TRUE)
    print(message('Masking elevation values greater than:', mask_land, terra::units(r)))
  }

  if (add_slope) {
    r <- terra::rast(bath_file)
    get_slope(bath_file, out_dir, neighbors = neighbors, slope_unit = slope_unit, overwrite = overwrite)
  }

  if (add_coastdist) {
    r <- terra::rast(bath_file)
    get_coastdist(bath_file, out_dir, mask_land = mask_land, dist_unit = dist_unit, overwrite = overwrite)
  }

}

# get_slope ---------------------------------------------------------------

#' Function to calculate slope from a bathymetry file
#'
#' @param bath_file File path to saved bathymetry raster. String.
#' @param out_dir File path to the directory where slope raster will be saved. String.
#' @param neighbors Number of neighbours to use in slope calculation. Either 4 or 8.
#' @param slope_unit Unit for the slope calculation. One of "degrees" or "radians".
#' @param overwrite Should existing files be overwritten. Logical.
#'
#' @details
#' This function is a wrapper for terra::terrain(), which will calculate slope
#' and save a standardized file to ensure consistency with other outputs from PelagicMagic.
#' The function can be used directly with a bathymetry file as input or it can be galled with get_bathymetry()
#' when doing the initial download of bathymetry data.
#'
#' See terra::terrain() for detailed description of how slope is calcualted.
#'
#' @return A netcdf file named 'slope.nc' saved to out_dir
#' @export
#' @seealso [get_bathymetry()], [terra::terrain()]
#'
#' @examples
#'get_bathymetry(out_dir = 'tmp',
#'               region = c(-68, 44, -63, 46),
#'               overwrite = FALSE,
#'               stride = 1,
#'               mask_land = 0)
#'get_slope(bath_file = 'tmp/bathymetry.nc',
#'          out_dir = 'tmp',
#'          neighbors = 8,
#'          slope_unit = 'degrees',
#'          overwrite = FALSE)
#'r <- terra::rast('tmp/slope.nc')
#'terra::plot(r)
#'unlink('tmp', recursive = TRUE)
#'
get_slope <- function(bath_file,
                      out_dir,
                      neighbors = 8,
                      slope_unit = c('degrees','radians'),
                      overwrite = FALSE
                      ) {

  if (file.exists(bath_file) == FALSE) stop(paste(bath_file), 'not found. Use get_bathymetry() to download bathymetry data.')
  if (neighbors %in% c(4,8) == FALSE) stop('neighbours must be one of 4 or 8')

  b <- terra::rast(bath_file)
  slope_file <- paste0(out_dir,'/slope.nc')

  if (file.exists(slope_file) == FALSE | overwrite == TRUE) {
    slope <- terra::terrain(b, neighbors = neighbors, v = 'slope', unit = slope_unit[1])
    terra::writeCDF(slope, slope_file,
                    unit = slope_unit[1],
                    longname = paste0('slope (',slope_unit[1],')'),
                    varname = 'slope',
                    overwrite = TRUE)
    message(paste('Saving slope to', slope_file))
  } else message(paste(slope_file, 'file already exists'))
}

# get_coastdist -----------------------------------------------------------

#' Function to calculate slope from a bathymetry file
#'
#' @param bath_file File path to saved bathymetry raster. String.
#' @param out_dir File path to the directory where slope raster will be saved. String.
#' @param mask_land Value of elevation that will be used to mask out land. Integer. Default is 0
#' @param dist_unit Unit for the distances One of "m" or "km".
#' @param overwrite Should existing files be overwritten. Logical.
#'
#' @details
#' This function is a wrapper for terra::costDist(). This function requires a bathymetry
#' file as input, where elevation values considered to be on land have been converted
#' to NA.
#'
#' If the input raster includes on land elevations, then the mask_land argument,
#' can be used to mask values greater than mask_land to NA.
#'
#' @return A netcdf file named 'coastdist.nc' saved to out_dir
#' @export
#' @seealso [get_bathymetry()], [terra::costDist()]

#'
#' @examples
#'get_bathymetry(out_dir = 'tmp',
#'               region = c(-68, 44, -63, 46),
#'               overwrite = FALSE,
#'               stride = 1,
#'               mask_land = 0)
#'get_coastdist(bath_file = 'tmp/bathymetry.nc',
#'              out_dir = 'tmp',
#'              mask_land = 0,
#'              dist_unit = 'km',
#'              overwrite = TRUE)
#'r <- terra::rast('tmp/coastdist.nc')
#'terra::plot(r)
get_coastdist <- function(bath_file,
                          out_dir,
                          mask_land = 0,
                          dist_unit = c('m','km'),
                          overwrite = FALSE
                          ) {
  if (file.exists(bath_file) == FALSE) stop(paste(bath_file), 'not found. Use get_bathymetry() to download bathymetry data.')

  b <- terra::rast(bath_file)
  if (is.null(mask_land) == FALSE) {
    if (is.numeric(mask_land) == FALSE | length(mask_land) > 1) stop('mask_land must be a single numeric value')
    b <- terra::ifel(b > mask_land, NA, b)
  }

  coast_file <- paste0(out_dir,'/coastdist.nc')
  if (file.exists(coast_file) == FALSE | overwrite == TRUE) {
    count_na <- sum(is.na(terra::values(b, na.rm = F)))
    if (sum(count_na) == 0) stop('Land values must be NA for get_coastdist(). Use mask_land argument to set land elevation.')
    r <- terra::ifel(is.na(b), 0, 1)
    d <- terra::costDist(r, target = 0)
    if (dist_unit[1] == 'km') d <- d/1000
    d <- terra::ifel(d == 0, NA, d)
    terra::longnames(d) <- 'distance from coast (km)'
    terra::varnames(d) <- 'coastdist'
    terra::writeCDF(d, coast_file,
                    unit = dist_unit[1],
                    varname = 'coastdist',
                    longname = paste0('distance from coast (',dist_unit[1],')'),
                    overwrite = T)
  } else (message(paste(coast_file, 'file already exists, skipping download')))
}
