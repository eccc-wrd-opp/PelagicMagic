# Functions for acquiring environmental data programmatically from common sources.

# get_bathymetry ----------------------------------------------------------

#' Function to download bathymetry data and calculate related layers: slope, coastdist, shelfdist
#'
#' @param out_dir File path to the directory where raster will be saved. String.
#' @param region Bounding box for region as: xmin, xmax, ymin, ymax.
#' @param overwrite Should existing files be overwritten. Logical.
#' @param stride Rate at which to sample underlying raster resolution. Integer.
#' @param vars List of bathymetry related variables to calculate. Default is c('slope', 'coastdist','shelfdist').
#' @param neighbors Number of neighbours to use in slope calculation. One of 4 or 8. Default is 8.
#' @param slope_unit Unit for the slope calculation. One of "degrees" or "radians". Default is degrees.
#' @param land_elev Value of elevation, in m, that will be used to mask out land when calculating coastdist. Integer. Default is 0
#' @param shelf_elev Value of elevation, in m, that will be used to calculate shelf distance. Integer. Default is -200 m
#' @param dist_unit Unit for the distance values One of "m" or "km". Default is m
#' @param maxiter The maximum number of iterations when calculating coastdist and shelfdist.
#' Increase this number if you get the warning that costDistance did not converge.
#'
#' @details
#'
#' This function uses rerddap to download the GEBCO_2020 elevation grid, which is a
#' continuous terrain model for oceans and land at 15 arc-second intervals (0.004 deg)
#' with elevation values in meters. The stride argument is passed to rerddap if
#' you want to down sample the resolution of the raster data before downloading.
#' For example, a stride of 2 will only acquire every second grid cell and a stride
#' of 10 will acquire every 10th grid cell.
#'
#' Vars is a list of common layers that can be derived from a bathymetry raster.
#' Current options include slope, distance from coastline (coastdist), and distance from shelf break (shlefdist).
#'
#' If vars includes slope, then slope is calculated using the terra::terrain() function.
#' The arguments neighbors and slope_unit are arguments passed to terra::terrain()
#' when calculating slope from a bathymetry layer. See ?terra::terrain for details.
#'
#' If vars includes coastdist, then distance to coast is calculated the terra::costDist() function.
#' The arguments land_elev and dist_unit are used when calculating coastdist.
#' All elevation values greater than land_elev will be converted to NA and coastdist
#' is calculated as distance to the closest NA value using terra::costDist(). The
#' distance raster will be returned in dist_unit (m or km).
#'
#' If vars includes shelfdist, then distance to shelf break is calculated using
#' the terra::costDist() function.The arguments shelf_elev and dist_unit are used
#' when calculating shelfdist. All elevation values greater than shelf_elev (i.e.
#' regions inshore of the shelf break) will be assigned negative distance values in dist_unit.
#' All elevation values less than shelf_elev (i.e. regions offshore of the shelf
#' break) will be assigned positive distance values in dist_unit. If shelf_elev is
#' outside of the range of elevation values in the bathymetry layer an error will be returned.
#'
#' @note
#'
#' When calculating shelfdist, distances close to the edge of the raster layer may
#' not reflect the true distance from the shelf break if the bathymetry layer is cropped.
#'
#' @return A netcdf file named 'bathymetry.nc' saved to out_dir. If additional
#' variables are included using vars, then additional files named will be saved.
#'
#' @export
#'
#' @seealso [terra::terrain()], [terra::costDist()], [rerddap::griddap()]
#'
#' @examples
#'# Gets all variables
#'get_bathymetry(out_dir = 'tmp',
#'               region = c(-68, -63, 44, 46),
#'               overwrite = FALSE,
#'               stride = 1,
#'               vars = c('slope','coastdist','shelfdist'),
#'               neighbors = 8,
#'               slope_unit = c('degrees','radians'),
#'               land_elev = 0,
#'               dist_unit = 'km',
#'               shelf_elev = -50)
#'r <- terra::rast('tmp/bathymetry.nc')
#'terra::plot(r, main = terra::longnames(r))
#'r <- terra::rast('tmp/slope.nc')
#'terra::plot(r, main = terra::longnames(r))
#'r <- terra::rast('tmp/coastdist.nc')
#'terra::plot(r, main = terra::longnames(r))
#'r <- terra::rast('tmp/shelfdist.nc')
#'terra::plot(r, main = terra::longnames(r))
#'
#' unlink('tmp', recursive = TRUE)
get_bathymetry <- function(out_dir = NULL,
                           region,
                           overwrite = FALSE,
                           stride = 1,
                           vars = c('slope','coastdist','shelfdist'),
                           neighbors = 8,
                           slope_unit = c('degrees','radians'),
                           land_elev = 0,
                           shelf_elev = -200,
                           dist_unit = c('m','km'),
                           maxiter = 50
) {

  bath_file <- paste0(out_dir,'/bathymetry.nc')

  if (file.exists(bath_file) == TRUE & overwrite == FALSE) {
    b <- terra::rast(bath_file)
    bext <- terra::ext(b)
    if (sum(sapply(1:4, function(x) round(bext[x]) == region[x]))!=4) {
      warning('Extent of existing bathymetry file does not match region. Do you want to overwrite existing files?',call. = FALSE)
    }
  }

  if (file.exists(bath_file) == FALSE | overwrite == TRUE) {
    b <- rerddap::griddap(rerddap::info('gebco2021', url = 'https://erddap.ifremer.fr/erddap'),
                          latitude = region[c(3,4)],
                          longitude = region[c(1,2)],
                          stride = stride,
                          store = rerddap::disk(path = out_dir,
                                                overwrite = TRUE)
    )
    file.rename(b$summary$filename, paste0(dirname(b$summary$filename), '/bathymetry.nc'))
    message(paste('Saving:', bath_file))
  } else message(paste('Skipping:', bath_file, 'exists in destination'))

  if ('slope' %in% vars) {
    get_slope(bath_file, out_dir, neighbors = neighbors, slope_unit = slope_unit, overwrite = overwrite)
  }

  if ('coastdist' %in% vars) {
    get_coastdist(bath_file, out_dir, land_elev = land_elev, dist_unit = dist_unit, overwrite = overwrite)
  }

  if ('shelfdist' %in% vars) {
    get_shelfdist(bath_file, out_dir, shelf_elev = shelf_elev, dist_unit = dist_unit, overwrite = overwrite)
  }


}

# get_slope ---------------------------------------------------------------

#' Internal function to calculate slope from a bathymetry file
#'
#' @param bath_file File path to saved bathymetry raster. String.
#' @param out_dir File path to the directory where slope raster will be saved. String.
#' @param neighbors Number of neighbours to use in slope calculation. Either 4 or 8.
#' @param slope_unit Unit for the slope calculation. One of "degrees" or "radians".
#' @param overwrite Should existing files be overwritten. Logical.
#'
#' @details
#' This function is a wrapper for terra::terrain(), which will calculate slope
#' and save a standardized file to ensure consistency with other outputs from get_bathymetry().
#' This is an internal function called by get_bathymetry().
#'
#' See terra::terrain() for detailed description of how slope is calculated.
#'
#' @return A netcdf file named 'slope.nc' saved to out_dir
#' @seealso [get_bathymetry()], [terra::terrain()]
#' @noRd
#'
get_slope <- function(bath_file,
                      out_dir,
                      neighbors = 8,
                      slope_unit = c('degrees','radians'),
                      overwrite = FALSE
) {

  if (file.exists(bath_file) == FALSE) stop(paste(bath_file), ' not found. Add bathymetry to vars argument to download bathymetry data.')
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

#' Internal function to calculate slope from a bathymetry file
#'
#' @param bath_file File path to saved bathymetry raster. String.
#' @param out_dir File path to the directory where slope raster will be saved. String.
#' @param land_elev Value of elevation that will be used to mask out land. Integer. Default is 0
#' @param dist_unit Unit for the distances One of "m" or "km".
#' @param overwrite Should existing files be overwritten. Logical.
#' @param maxiter The maximum number of iterations when calculating coastdist. Increase this number if you get the warning that costDistance did not converge
#'
#' @details
#' This function is a wrapper for terra::costDist(). This is an internal function
#' called by get_bathymetry().
#'
#' @return A netcdf file named 'coastdist.nc' saved to out_dir
#' @noRd
#' @seealso [get_bathymetry()], [terra::costDist()]
#'
get_coastdist <- function(bath_file,
                          out_dir,
                          land_elev = 0,
                          dist_unit = c('m','km'),
                          overwrite = FALSE,
                          maxiter = 50
) {
  if (file.exists(bath_file) == FALSE) stop(paste(bath_file), ' not found. Add bathymetry to vars argument to download bathymetry data.')

  b <- terra::rast(bath_file)
  if (is.null(land_elev) == FALSE) {
    if (is.numeric(land_elev) == FALSE | length(land_elev) > 1) stop('land_elev must be a single numeric value')
    b <- terra::ifel(b > land_elev, NA, b)
  }

  coast_file <- paste0(out_dir,'/coastdist.nc')
  if (file.exists(coast_file) == FALSE | overwrite == TRUE) {
    count_na <- sum(is.na(terra::values(b, na.rm = F)))
    if (sum(count_na) == 0) stop('Land values must be NA for get_coastdist(). Use land_elev argument to set land elevation.')
    r <- terra::ifel(is.na(b), 0, 1)
    d <- terra::costDist(r, target = 0, overwrite = TRUE, maxiter = maxiter)
    if (dist_unit[1] == 'km') d <- d/1000
    d <- terra::ifel(d == 0, NA, d)
    terra::longnames(d) <- 'distance from coast (km)'
    terra::varnames(d) <- 'coastdist'
    terra::writeCDF(d, coast_file,
                    unit = dist_unit[1],
                    varname = 'coastdist',
                    longname = paste0('distance from coast (',dist_unit[1],')'),
                    overwrite = TRUE)
  } else (message(paste(coast_file, 'file already exists, skipping download')))
}

# get_coastdist -----------------------------------------------------------

#' Internal function to calculate slope from a bathymetry file
#'
#' @param bath_file File path to saved bathymetry raster. String.
#' @param out_dir File path to the directory where slope raster will be saved. String.
#' @param land_elev Value of elevation that will be used to mask out land. Integer. Default is 0
#' @param dist_unit Unit for the distances One of "m" or "km".
#' @param overwrite Should existing files be overwritten. Logical.
#' @param maxiter The maximum number of iterations when calculating shelfdist. Increase this number if you get the warning that costDistance did not converge
#'
#' @details
#' This function is a wrapper for terra::costDist(). This is an internal function
#' called by get_bathymetry().
#'
#' @return A netcdf file named 'coastdist.nc' saved to out_dir
#' @noRd
#' @seealso [get_bathymetry()], [terra::costDist()]
#'
get_shelfdist <- function(bath_file,
                          out_dir,
                          shelf_elev = -200,
                          dist_unit = c('m','km'),
                          overwrite = FALSE,
                          maxiter = 50
) {

  if (file.exists(bath_file) == FALSE) stop(paste(bath_file), ' not found. Add bathymetry to vars argument to download bathymetry data.')

  b <- terra::rast(bath_file)

  b_range <- range(terra::values(b, na.rm = T))
  if (shelf_elev < b_range[1]) stop(paste0('shelf_elev [',shelf_elev,' m] is below minimum bathymetry [',b_range[1],' m]'))
  if (shelf_elev > b_range[2]) stop(paste0('shelf_elev [',shelf_elev,' m] is above maximum bathymetry [',b_range[1],' m]'))

  shelf_file <- paste0(out_dir,'/shelfdist.nc')

  if (file.exists(shelf_file) == FALSE | overwrite == TRUE) {

    rn <- terra::ifel(b < shelf_elev, 0, 1)
    inshore <- terra::costDist(rn, target = 0, overwrite = TRUE, maxiter = maxiter)
    if (dist_unit[1] == 'km') inshore <- inshore/1000
    inshore <- -inshore
    rn <- terra::ifel(b > shelf_elev, 0, 1)
    offshore <- terra::costDist(rn, target = 0, overwrite = TRUE, maxiter = maxiter)
    if (dist_unit[1] == 'km') offshore <- offshore/1000
    d <- inshore + offshore
    names(d) <- 'shelfdist'
    terra::writeCDF(d, shelf_file,
                    unit = dist_unit[1],
                    varname = 'shelfdist',
                    longname = paste0('Distance from shelf break [', shelf_elev,' m] (',dist_unit[1],')'),
                    overwrite = TRUE)
  } else (message(paste(shelf_file, 'file already exists, skipping download')))
}




# check_colony ---------------------------------------------------------

#' Function to check if colony coordinates touch water and, if not, adjust coordinates
#' to nearest water point. Internal function for within get_colonydist
#'
#' @param landmask A SpatRaster object. Land is converted to NA and water to 1 within function.
#' @param coords Data frame with columns labelled 'location', 'x', and 'y'. Where location is the
#' @param search_dist Optional search distance (in m) to look for water locations if a location is located on land.
#' location name used for saving each distance raster, x and y are the location coordinates in the same units as landmask.
#'
#' @returns An adjusted colony coordinate data frame to be used in get_colonydist.
#'
#' @noRd
#' @seealso [get_colonydist()]]
#'
#' @examples
#'
#' # Get bathymetry
#' get_bathymetry(out_dir = 'tmp',
#'                region = c(-68, -63, 44, 46),
#'                overwrite = FALSE,
#'                stride = 1,
#'                vars = NULL)
#'
#' m <- terra::rast('tmp/bathymetry.nc')
#'
#' my_cols <- data.frame(location = c('loc1', 'loc2'),
#'                       x = c(-66.718, -66.387),
#'                       y = c(44.702, 44.240)
#' )
#'
#' get_colonydist(landmask = m,
#'                coords = my_cols,
#'                out_dir = 'tmp/colony_dist',
#'                dist_unit = c('km'),
#'                search_dist = NULL,
#'                plot = TRUE,
#'                overwrite = TRUE)
#'

check_colony <- function(landmask, coords, search_dist = NULL) {

  # # Make sure landmask is in a binary style (land = NA and water = 1)
  # landmask <- terra::ifel(landmask >=0, NA, 1)

  # Check if the extracted value is NA
  if (is.na(terra::extract(landmask, coords[,c('x','y')], ID = FALSE)[,1])) {

    # If NA, find the closest 1 value
    # Create a mask for all cells with value 1
    mask_1 <- landmask == 1

    #convert loc to SpatVector
    ml <- terra::vect(coords, geom=c("x", "y"), crs="epsg:4326")

    if (!is.null(search_dist)) {
      # Crop to within search_dist for speed
      loc_buff <- terra::buffer(ml, search_dist)
      mask_1 <- terra::crop(mask_1, loc_buff)
    }

    # Select closest overwater location
    distances <- terra::distance(mask_1, ml)
    distances <- terra::mask(distances, mask_1)
    names(distances) <- 'distance'

    coords_new <- data.frame(
      location = coords[,1],
      terra::crds(distances, na.rm = T),
      terra::values(distances, na.rm = T)
    ) |>
      dplyr::slice_min(distance) |>
      dplyr::select(names(coords)) |>
      data.frame()

    # Output the new coordinates
    print(paste0('Shifting ', coords[,1], ' to be over water. New location at: ', round(coords_new[,2],3), ' , ',round(coords_new[,3],3)))
    coords_new

  } else {
    # If the extracted value is 1, just return the original coordinates
    coords
  }

}

# get_colonydist ---------------------------------------------------------

#' Function to calculate overwater distance from a colony (single point location or set of point locations)
#'
#' @param landmask A SpatRaster object. Raster cells that are inaccessible (land) should have NA values. All other cells
#' should have values corresponding to the cost of travelling over that location. In most cases this will be values of
#' 1 (water), however if you want to account for species avoiding certain areas (like ice) cells can have weights to
#' increase the cost of travel. See example for converting bathymetry to NA/1 raster before running the function.
#' @param coords Data frame with columns labelled 'location', 'x', and 'y'. Where location is the
#' location name used for saving each distance raster, x and y are the location coordinates in the same units as landmask.
#' @param out_dir File path to the directory where slope raster will be saved. String.
#' @param dist_unit Unit for the distances One of "m" or "km".
#' @param search_dist Optional search distance (in m) to look for water locations if a location is located on land.
#' @param plot Should the distance rasters be plotted. Logical.
#' @param overwrite Should existing files be overwritten. Logical.
#'
#' @return A separate ncdf file with overwater distance calculate for each record in coords. Output saved to out_dir.
#' @export
#'
#' @examples
#'
#' # Get bathymetry
#' get_bathymetry(out_dir = 'tmp',
#'                region = c(-68, -63, 44, 46),
#'                overwrite = FALSE,
#'                stride = 1,
#'                vars = NULL)
#'
#' m <- terra::rast('tmp/bathymetry.nc')
#' m <- terra::ifel(m >=0, NA, 1)
#'
#' my_cols <- data.frame(location = c('loc1', 'loc2'),
#'                       x = c(-66.718, -66.387),
#'                       y = c(44.702, 44.240)
#' )
#'
#' get_colonydist(landmask = m,
#'                coords = my_cols,
#'                out_dir = 'tmp/colony_dist',
#'                dist_unit = c('km'),
#'                search_dist = NULL,
#'                plot = TRUE,
#'                overwrite = TRUE)
#'
#' unlink('tmp', recursive = TRUE)
#'
#'
#'
get_colonydist <- function(landmask, coords, out_dir, dist_unit = c('m','km'),
                           search_dist = NULL, plot = TRUE, overwrite = FALSE) {

  if (class(landmask)[1] != 'SpatRaster') stop('landmask must be a SpatRaster', call. = FALSE)
  if (!('x' %in% names(coords)) | !('y' %in% names(coords)) | !('location' %in% names(coords))) stop('coords must include columns named: location, x, and y', call. = FALSE)
  if (dir.exists(out_dir) == FALSE) dir.create(out_dir, recursive = TRUE)

  #coords$idx <- terra::cellFromXY(landmask, xy = coords[,c('x','y')])

  for (i in 1:nrow(coords)) {

    out_file <- paste0(out_dir,'/',coords$location[i], '.nc')
    if (file.exists(out_file) == FALSE | overwrite == TRUE) {

      # Check and then adjust (if necessary) the location of the colony
      coords_adj <- check_colony(landmask = landmask, coords = coords[i,], search_dist = search_dist)
      coords_adj$idx <- terra::cellFromXY(landmask, xy = coords_adj[,c('x','y')])

      # # Convert land to NA and water to 1
      # if(max(terra::values(landmask)) > 1) {lm <- terra::ifel(landmask >=0, NA, 1)
      #       } else(lm <- landmask)

      m <- landmask
      terra::values(m)[coords_adj$idx] <- -1
      d <- terra::costDist(m, target = -1, overwrite = TRUE)
      terra::varnames(d) <- 'distance'
      terra::longnames(d) <- 'Distance (m)'
      terra::units(d) <- 'm'

      if (dist_unit[1] == 'km') {
        d <- d/1000
        terra::longnames(d) <- 'Distance (km)'
        terra::units(d) <- 'km'
      }

      if (plot == TRUE) {
        terra::plot(d, main = coords_adj$location)
        terra::points(coords_adj[, c('x','y')], pch = 19, col = 'red')
      }

      terra::writeCDF(d,
                      out_file,
                      overwrite = overwrite,
                      varname = coords_adj$location,
                      longname = ifelse(dist_unit[1] == 'km', 'Overwater distance (km)', 'Overwater distance (m)'),
                      unit = ifelse(dist_unit[1] == 'km', 'km', 'm')
                      )
    } else print(paste('Skipping:', coords_adj$location))
  }
}

