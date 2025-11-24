
# get_flocks --------------------------------------------------------------

#' Create flocks from census observations
#'
#' @param dat An sf points object of survey observations. Must include a column named 'total' that records the number of birds seen per observation. Must be in a projected coordinate system.
#' @param flock_spacing Standard deviation in m for generating birds around previous observation
#' @param search_dist Standard deviation in m for generating birds when trying to move birds over land onto water
#' @param np Number of random points to generate when trying to move points on land onto water.
#' @param coast An sf polygon object of land, which will be used to ensure that birds are always lcoated in water. These polygons could be configured to exclude birds from any areas (e.g. exclude land birds from ocean).
#'
#' @details
#' In aerial census surveys of seaducks, a single point is recorded for each flock observed.
#' This function is the first step in a workflow for identifying hotspots from census surveys.
#' The function generates a single row for every individual in the flock. Because
#' very large flocks would have occupied more geographic space than small flocks, the
#' function adds individuals around the first observation using a random walk with the
#' flock_spacing argument determining the distribution of distances between subsequent points.
#'
#' The coast argument allows users to specify areas (land) that are inaccessible to the
#' focal species. Any observations that end up overlapping the coast layer are moved to the closest
#' randomly generated point within search_distance of the overlapping point. Currently there is no
#' option to run this function without a 'coast' layer. The example uses the rnaturalearth::ne_countries()
#' layer, however in mnay applications for coastal species this shoreline may be too coarse and the
#' user should consider a layer wih higher resolution.
#'
#' This function can be quite slow for large data sets. User should crop their coast layer to
#' focus on the area of interest (see example)
#' #'
#' @return An sf POINT object with one row for each bird in observed flocks. If
#' some birds are still overlapping with the coast polygons a warning will be given.
#' The output will include all fields from dat, plus a flock_id field indicating
#' which row the observation is associated with.

#'
#' @examples
#'
#' # format data
#' d <- coei
#' d <- sf::st_as_sf(d, coords = c('longitude', 'latitude'), crs = 4326)
#' d <- sf::st_transform(d,
#'                       crs = '+proj=laea +lat_0=50 +lon_0=-65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
#'
#' # get coast layer
#' coast <- rnaturalearth::ne_countries(scale = 10, return = 'sf')
#' coast <- sf::st_transform(coast, crs = sf::st_crs(d))
#' coast_crop <- sf::st_crop(coast, sf::st_buffer(d, dist = 10000))
#'
#' # generate flocks from observation data
#' f <- get_flocks(dat = d,
#'                 flock_spacing = 10,
#'                 search_dist = 100,
#'                 np = 100,
#'                 coast = coast_crop)
#'
#' sum(d$total) == nrow(f)
#' f
#'
#' @export
get_flocks <- function(dat,
                       flock_spacing,
                       search_dist = 100,
                       np = 100,
                       coast) {

  if (class(dat)[1] != "sf") stop('dat must be an sf object', call. = F)
  if (('total' %in% names(dat)) == F)  stop('dat must include a column named total, indicating the number of birds observed', call. = F)
  if (class(dat$total) != 'numeric')  stop('the total field must be numeric', call. = F)
  if (sf::st_is_longlat(dat)) stop('dat must be in a projected coordinate system', call. = F)
  if (class(coast)[1] != "sf") stop('coast must be an sf object', call. = F)


  dat <- dat |> dplyr::filter(total >0)

  cc <- suppressWarnings(sf::st_crop(coast, dat))

  coords <- sf::st_coordinates(dat)

  out <- lapply(1:nrow(dat), function(i) {

    # generate np random points around flock location within maxdist
    # this is to move flocks on land to the closest location over water
    pp <-  data.frame(
      X = rnorm(np, coords[i,1], search_dist),
      Y = rnorm(np, coords[i,2], search_dist)
    ) |>
      sf::st_as_sf(coords = c('X',"Y"), crs = sf::st_crs(dat))


    # select points that are not on land
    pp <- pp[lengths(sf::st_intersects(pp, cc, sparse = T)) == 0,]

    # only generate birds for flocks that are within maxdist of water
    if (nrow(pp) > 0) {

      # find the closest point to the original observation that is over water
      d <- sf::st_distance(pp, dat[i,])[,1]
      idx <- which.min(d)[1]
      temp <- data.frame(sf::st_coordinates(pp[idx,]))

      # generates birds within flock using a random walk
      k <- 2
      while (k <= dat$total[i]) {
        temp <- rbind(temp,
                      data.frame(
                        X = rnorm(1, temp$X[k - 1], flock_spacing),
                        Y = rnorm(1, temp$Y[k - 1], flock_spacing)
                      ))
        k <- k + 1
      }


      # Add flock index
      suppressWarnings(data.frame(
        sf::st_drop_geometry(dat[i,]),
        temp,
        flock_id = i
      ))
    }
  })

  out <- do.call(rbind, out)
  out <- sf::st_as_sf(out, coords = c('X','Y'), crs = sf::st_crs(dat))
  # checks for birds on land
  cr <- suppressWarnings(sf::st_crop(coast, out))
  idx <- which(lengths(sf::st_intersects(out, cr, sparse = T)) > 0)

  # Moves birds on land to water
  coords <- sf::st_coordinates(out)
  qq <- lapply(idx, function(i) {
    tt <- data.frame(
      X = rnorm(np, coords[i,1], search_dist),
      Y = rnorm(np, coords[i,2], search_dist)
    ) |>
      sf::st_as_sf(coords = c('X',"Y"), crs = sf::st_crs(out))
    tt <- tt[lengths(sf::st_intersects(tt, cr, sparse = T)) == 0,]
    d <- sf::st_distance(tt, out[i,])[,1]
    rep_idx <- which.min(d)[1]
    tt$geometry[rep_idx]
  })
  qq <- do.call(rbind, qq)
  out$geometry[idx] <- qq

  if (nrow(out) != sum(dat$total)) warning('Some observations still overlapped with coast layer, consider increasing search_dist')

  out

}


# get_ud_spatstat --------------------------------------------------------------

get_ud_spatstat <- function(data,
                            h,
                            res = 100,
                            ext = NULL,
                            factor = 1,
                            coast = coast
) {

  if (is.null(ext)) {
    my_ext = terra::ext(data)[1:4] # extent for UD
  } else my_ext <- ext
  br = terra::rast(ext = my_ext, res = res, crs = terra::crs(data)) # template raster
  terra::values(br) <- 1
  br <- terra::mask(br, sf::st_buffer(coast, -(res * factor)), inverse=TRUE, touches = T)
  my_win <- spatstat.geom::owin(mask = data.frame(terra::crds(br, na.rm = T)))

  cc <- sf::st_coordinates(data)
  my_ppp <- spatstat.geom::ppp(x = cc[,1],
                               y = cc[,2],
                               window=my_win)
  dd <- spatstat.explore::density.ppp(my_ppp, sigma = h, positive=FALSE)
  dd <- terra::rast(dd)
  terra::crs(dd) <- terra::crs(data)
  names(dd) <- 'ud'
  dd
}

# classify and smooth raster ----------------------------------------------
# function to classify rasters and make smoothed polygons for each level that includes the higher levels polygon

hs_class <- function(r, lvl, method = 'ksmooth', donuts = T) {
  names(r)<- 'class'
  mc <- lapply(lvl, function(x) {
    tt <- terra::ifel(r >= x,x,NA)
    tt <- terra::as.polygons(tt)
    tt <- smoothr::smooth(sf::st_as_sf(tt), method = method) |>
      sf::st_as_sf() |>
      sf::st_cast('MULTIPOLYGON')
    tt
  })
  mc <- do.call(rbind, mc)
  mc$class <- factor(mc$class)
  mc <- mc |>
    sf::st_cast('POLYGON')

  # add pretty labels
  lev_inf<- c(lvl, Inf)
  labs <- paste0(lev_inf, '-', dplyr::lead(lev_inf))[1:length(lvl)]
  labs[length(labs)] <- paste0('>',lvl[length(labs)])
  mc$class <- factor(mc$class, labels = labs[1:length(unique(mc$class))])
  row.names(mc) <- NULL

  # erase insides of overlapping polygon levels
  if (donuts == T) {
    mc <- lapply(1:length(labs), function(f) {
      if (f < length(labs)) {
        suppressWarnings(erase_poly(mc |> dplyr::filter(class == labs[f]),
                                    mc |> dplyr::filter(class == labs[f + 1])))
      } else {
        mc |> dplyr::filter(class ==  labs[f])
      }
    })
    mc <- do.call(rbind, mc)
  }

  mc
}


# erase_poly --------------------------------------------------------------
#' Erase parts of polygons
#'
#' @param x An sf POLYGONS object
#' @param y An sf POLYGONS object
#'
#' @description A helper function that errases any part of x that overlaps with y
#'
#' @return An sf POLYGONS object
#' @export
#' @examples
#'
#' library(ggplot2)
#' d <- coei
#' d <- sf::st_as_sf(d, coords = c('longitude', 'latitude'), crs = 4326)
#' d <- sf::st_transform(d,
#'                       crs = '+proj=laea +lat_0=50 +lon_0=-65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
#'
#' # get coast layer
#' coast <- rnaturalearth::ne_countries(scale = 10, return = 'sf')
#' coast <- sf::st_transform(coast, crs = sf::st_crs(d))
#' coast_crop <- sf::st_crop(coast, sf::st_buffer(d, dist = 10000))
#'
#' bbox_poly <-  d |>
#' sf::st_buffer(dist = 10000) |>
#'   sf::st_bbox() |>
#'   sf::st_as_sfc()
#'
#' p <- erase_poly(bbox_poly, coast_crop)

erase_poly <- function(x, y) {sf::st_difference(x, sf::st_union(sf::st_combine(y)))}


# get_hotspots ---------------------------------------------------------

#' Identify hot spots from aerial census data
#'
#' @param obs An sf points object of survey observations with one row per individual, must include a column named 'date' of class Date or POSIX.
#' @param regions An optional sf polygon object to split flocks into regions for faster analysis.
#' @param coast An sf polygons object of land areas to mask from utilization distribution.
#' @param time_group A date abbreviation to indicate how surveys should be grouped for analysis. One of %Y, %m, %d.
#' @param res Spatial raster resolution for utilization distribution in meters.
#' @param h Smoothing parameter for utliization distribution in meters.
#' @param sum_stat Statistic for summarizing distribution across survey events, one of 'max', 'mean', or 'median'.
#' @param lev List of breakpoints for summarizing the raster distribution into polygons. Values below first level will be dropped.
#' @param buffer_points Distance to buffer points when creating study area polygon if regions argument is NULL.
#' @param area_unit Unit of measurement for area calculations, one of 'km' or 'm'.
#' @param donuts Logical. Should higher level polygons be erased from lower level polygons to create donuts.
#' @param summaries Logical. Should summaries of bird observations within each polygon be returned.
#' @param output Should the function return summary polygons, summary rasters or all weighted utlization rasters, one of 'polygons', 'summary_rast', or 'all_rast'.
#'
#' @details This function calculates a utilization distribution (UD) for each unique level in time_group.
#' The observations can be broken up into regions to reduce the size of the raster on which the UD is calculated.
#' Each raster is then weighted by the sum of all observations within the time_group+region combination. Rasters
#' are mosaiced together and summarized using the function option provided in sum_stat.
#'
#' Continuous outputs from the rasters can be converted into polygons by classifying
#' observations according to the breakpoints in lev. This creates a smoothed polygon
#' representing areas with weighted utlization distribution values within each class.
#' Polygons can be converted to donuts with higher levels of expected distribution
#' removed from lower level polygons. Summaries are provided for the size of the polygons, and the
#' mean and max counts observed within polygons and the number of surveys conducted of that area.
#'
#'
#'
#' @return An sf POLYGONS object showing hotspots
#' @export
#'
#' @examples
#' # format data
#' d <- coei
#' d <- sf::st_as_sf(d, coords = c('longitude', 'latitude'), crs = 4326)
#' d <- sf::st_transform(d,
#'                       crs = '+proj=laea +lat_0=50 +lon_0=-65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
#'
#' # get coast layer
#' coast <- rnaturalearth::ne_countries(scale = 10, return = 'sf')
#' coast <- sf::st_transform(coast, crs = sf::st_crs(d))
#' coast_crop <- sf::st_crop(coast, sf::st_buffer(d, dist = 10000))
#'
#' # generate flocks from observation data
#' f <- get_flocks(dat = d,
#'                 flock_spacing = 10,
#'                 search_dist = 100,
#'                 np = 100,
#'                 coast = coast_crop)
#'
#' reg <- f |>
#'   sf::st_buffer(dist = 10000) |>
#'   sf::st_bbox() |>
#'   sf::st_as_sfc() |>
#'   sf::st_as_sf()
#'
#' # generate hoptpots
#' hs <- get_hotspots(
#'   obs = f,
#'   regions = reg,
#'   coast = coast_crop,
#'   time_group = c('%Y', '%m', '%d')[1],
#'   res = 250,
#'   h = 250,
#'   sum_stat = c('max','mean','median')[1],
#'   lev = c(1, 10, 100, 1000, 10000),
#'   buffer_points = 10000,
#'   area_unit = c('km','m')[1],
#'   donuts = T,
#'   summaries = T
#' )
#' hs
#'

get_hotspots <- function(
    obs,
    regions = NULL,
    coast,
    time_group = c('%Y', '%m', '%d')[1],
    res = 250,
    h = 250,
    sum_stat = c('max','mean','median')[1],
    lev = c(1, 10, 100, 1000, 10000),
    buffer_points = 10000,
    area_unit = c('km','m')[1],
    donuts = T,
    summaries = T,
    output = c('polygons','summary_rast','all_rast')[1]
){

  if (!('date' %in% names(obs))) stop('obs must include a column named date with class Date or POSIXct', call. = F)
  if (!(class(obs$date)[1] %in% c("Date",'POSIXct','POSIXlt'))) stop('obs must include a column named date with class Date or POSIXct', call. = F)
  if (class(obs)[1] != "sf") stop('obs must be an sf object', call. = F)
  if (class(coast)[1] != "sf") stop('coast must be an sf object', call. = F)
  if (!(time_group %in% c('%Y', '%m', '%d'))) stop('time_group must be one of %Y, %m, %d', call. = F)
  if (!(area_unit %in% c('km', 'm'))) stop('area_unit must be one of km, m', call. = F)
  if (!(sum_stat %in% c('max','mean','median'))) stop('sum_stat must be one of max, mean, median', call. = F)


  # create regions if missing
  if (is.null(regions)) {
    regions <-  obs |>
      sf::st_buffer(dist = buffer_points) |>
      sf::st_bbox() |>
      sf::st_as_sfc() |>
      sf::st_as_sf()
  }

  # create time grouping variable
  if (is.null(time_group)) {
    obs$tg <- 1
  } else {
    obs$tg <- strftime(obs$date, time_group)
  }

  out <- lapply(1:nrow(regions), function(i) {

    p <- regions[i,]
    tr <- suppressWarnings(obs |> sf::st_intersection(p))

    my_ext = terra::ext(tr)[1:4] # extent for UD
    br = terra::rast(ext = my_ext, res = res, crs = terra::crs(obs)) # template raster

    rd <- lapply(unique(tr$tg), function(y) {

      tt <- tr |> dplyr::filter(tg == y)

      my_hr <- get_ud_spatstat(data = tt,
                               ext = terra::ext(p)[1:4],
                               h = h,
                               res = res,
                               factor = 1,
                               coast = coast)
      (my_hr$ud * nrow(tt))/terra::cellSize(my_hr, unit = area_unit)

    })

    rd <- terra::rast(rd)
    names(rd) <- unique(tr$tg)
    rd <- terra::mask(rd, p)
    rm <- terra::app(rd, sum_stat)
    rm
  })

  m <- terra::sprc(out)

  if (output == 'all_rast') return(m)
  m <- terra::mosaic(m, fun = sum_stat)
  if (output == 'summary_rast') return(m)

  #convert to  polygon classes
  mc <-suppressWarnings(hs_class(r = m, lvl = lev, donuts = donuts))

  if (summaries == T) {
    my_ext <- terra::ext(obs)[1:4] # extent for UD
    br <- terra::rast(ext = my_ext, res = res, crs = terra::crs(obs)) # template raster
    tt <- terra::rasterize(obs, br, fun = 'sum', by = 'tg', background = 0)
    ex <- terra::extract(tt, mc, fun = 'sum', na.rm = T, exact = T, ID = F)

    #mc$surveys <- ncol(ex)
    mc$max_cnt <- round(apply(ex, 1, max))
    mc$mean_cnt <- round(apply(ex, 1, mean))
    mc$surv_cnt <- round(apply(ex, 1, function(x) sum(x>0)))
    mc$area_km2 <- as.numeric(units::set_units(sf::st_area(mc), km^2))


    mc <- mc |>
      dplyr::select(class, area_km2, mean_cnt, max_cnt, yrs_obs)
  }
  return(mc)
}

