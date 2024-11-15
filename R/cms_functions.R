
# make_cms_envir ----------------------------------------------------------

# Create a python environment with the copernicusmarine toolbox installed

#' Create a Python virtual environment with copernicusmarine toolbox installed
#'
#' @param my_env The name of, or path to, a Python virtual environment. If this
#' name contains any slashes, the name will be interpreted as a path; if the name
#' does not contain slashes, it will be treated as a virtual environment within virtualenv_root().
#' When NULL, the virtual environment as specified by the RETICULATE_PYTHON_ENV
#' environment variable will be used instead. To refer to a virtual environment
#' in the current working directory, you can prefix the path with ⁠./<name>⁠.
#' @param version Python version to use. String. Default is '3.11'.
#' @param user User name for copernicus marine service. String.
#' @param password Password for copernicus marine service. String.
#'
#' @details This function uses reticulate to create a Python virtual environment
#' and install the copernicusmarine toolbox. This environment can then be used to download
#' oceanographic data from Copernicus Marine Service.
#'
#' If the function is run including the user name and password then login details
#' will be saved to the virtual environment on creation.
#'
#' @return Creates a python virtual environment with the copernicusmarine toolbox installed on the path specified by my_env
#' @export
#'
#' @examples
#'
#' # not run
#' if (FALSE) {
#' cm_env <- './tmp/cms_env'
#' make_cms_envir(my_env = cm_env)
#' reticulate::use_virtualenv(paste0(my_env,'/'), required = TRUE)
#' cm <- reticulate::import("copernicusmarine")
#' }
#'
make_cms_envir <- function(my_env, version = '3.11', user = NULL, password = NULL) {
  if (reticulate::virtualenv_exists(my_env) == F) {
    reticulate::install_python(version = version, list = FALSE, optimized = TRUE)
    reticulate::virtualenv_create(envname = my_env, ignore_installed = TRUE)
    reticulate::virtualenv_install(envname = my_env, packages = c("copernicusmarine"))

  } else message(paste('Python environment already exists at:', my_env))

  if (!is.null(user) & !is.null(password)) {
    reticulate::use_virtualenv(paste0(my_env,'/'), required = TRUE)
    cm <- reticulate::import("copernicusmarine")
    cm$login(user, password, skip_if_user_logged_in = T)
  }
}

# get_cms_data ------------------------------------------------------------

#' Download subsetted data from Copernicus Marine Service (CMS)
#'
#' @param my_env Path to a virtual python environment containing the Copernicus Marine Toolbox
#' @param out_dir File path to the directory where raster will be saved. String.
#' @param var Variable name to be downloaded from CMS
#' @param prod CMS product ID
#' @param freq Frequency of the dataset to download. Currently supports: daily, monthly, and static.
#' @param region Bounding box for region as: xmin, xmax, ymin, ymax.
#' @param date_min The earliest date to include, formatted as a date class.
#' @param date_max The latest date to include, formatted as a date class.
#' @param depth_min For datasets with a depth component, the minimum depth value to include. Default is 0.
#' @param depth_max For datasets with a depth component, the maximum depth value to include. Default is 0.
#' @param overwrite TRUE/FALSE. Should existing files be overwritten? Default is FALSE.
#'
#' @details  This is a wrapper function to use the subset function of the Copernicus
#' Marine Toolbox directly from R with reticulate. CMS provides access to many useful
#' oceanographic datasets. This function can be used to download certain commonly
#' used datasets programmatically. Because the full data files cover the entire globe,
#' often at multiple depth layers making it impractical to download the full original files
#' for large time ranges. The subsetting tool allows us to request custom spatial, temporal, and depth
#' extents.
#'
#' This function relies on a data table that is pre-populated with required fields for the
#' datasets that are most commonly used in species distribution modelling: Global
#' Ocean Physics and Global Ocean Biogeochemistry. Other datasets could be added on request.
#' The hindcasast nd forecast models have different dataset IDs and cover different date ranges,
#' so this function requests data from the correct datasets for the temporal extent requested.
#' If the request covers multiple datasets, then a netcdf will be saved for each dataset.
#'
#' Before running this function for the first time you will need to create a Python
#' virtual environment with the copernicus marine toolbox installed. This can be
#' done using the get_cms_env() function. This environment could be saved outside
#' of an R project and called by any project that requires CMS data.
#'
#'
#' @return Netcdf file(s) saved to out_dir.
#' @export
#'
#' @seealso make_cms_envir()
#'
#' @examples
#'
#' data('cms_table')
#' cms_table
#'
#' # not run
#' if (FALSE) {
#' get_cms_data(
#'    my_env = cm_env,
#'    out_dir = 'tmp',
#'    var = 'so',
#'    freq = 'monthly',
#'    region = c(-68, -63, 44, 46),
#'    date_min = Sys.Date() + -31,
#'    date_max = Sys.Date() + 7,
#'    depth_min = 0,
#'    depth_max = 0
#'    )
#'}


get_cms_data <- function(
    my_env,
    out_dir,
    var = 'thetao',
    prod = 'cmems_mod_glo_phy_hindcast',
    freq = c('monthly', 'daily', 'static'),
    region = c(-68, -63, 44, 46),
    date_min = Sys.Date() - 1,
    date_max = Sys.Date(),
    depth_min = 0,
    depth_max = 0,
    overwrite = FALSE
) {

  reticulate::use_virtualenv(paste0(my_env,'/'), required = F)
  cm <- reticulate::import("copernicusmarine")

  if (!(prod %in% cms_table$prod)) stop(paste(prod, 'not in the included datasets, check cms_table'))

  vv <- cms_table[cms_table$prod == prod & cms_table$variable == var & cms_table$frequency == freq,]
  vv$max_date[is.na(vv$max_date)] <- Sys.Date() + 14
  vv <- vv[vv$max_date >=date_min,]
  vv <- vv[vv$min_date <=date_max,]

  if (nrow(vv) == 0) stop('No datasets matched requrested criteria for prod, var, and freq. Check cms_table')

  for (i in 1:nrow(vv)) {
    of <- paste0(out_dir,'/',vv$frequency[i],'/',vv$variable[i],'/',vv$dataset_name[i],
                 '-',vv$frequency[i],'-',vv$variable[i]
    )
    on <- paste0(vv$dataset_name[i],
                 '-',vv$frequency[i],'-',vv$variable[i]
    )

    if (file.exists(of) == FALSE | overwrite == T) {
      cm$subset(
        dataset_id=vv$dataset_id[i],
        variables=list(vv$variable[i]),
        minimum_longitude=region[1],
        maximum_longitude=region[2],
        minimum_latitude=region[3],
        maximum_latitude=region[4],
        start_datetime=paste0(as.Date(max(c(date_min, vv$min_date[i]), na.rm = T)),"T00:00:00"),
        end_datetime=paste0(as.Date(min(c(date_max, vv$max_date[i]), na.rm = T)),"T00:00:00"),
        minimum_depth=depth_min,
        maximum_depth=depth_max,
        output_filename = on,
        output_directory = paste0(out_dir,'/',vv$frequency[i],'/',vv$variable[i]),
        service = 'arco-geo-series',
        force_download = TRUE,
        overwrite_output_data = TRUE
      )

    } else print(paste('Skipping:', of))
  }

  rm(cm)
  rm(my_env)
}


# cms_table ---------------------------------------------------------------

#' Copernicus Marine Service product table
#'
#' Metadata information about Copernicus Marine Service datasets that can be downloaded using get_cms_data()
#'
#' @format A data frame
#' \describe{
#'   \item{variable}{CMS name of the variable that will be downloaded}
#'   \item{frequency}{Temporal frequency of the data to be downloaded}
#'   \item{product_id}{CMS product ID}
#'   \item{dataset_id}{CMS dataset ID}
#'   \item{dataset_name}{Name used when saving outputs from this dataset}
#'   \item{min_date}{Earliest date available in this dataset}
#'   \item{max_date}{Latest date available in this dataset}
#' }
#' @author Allison Patterson \email{allison.patterson@ec.gc.ca}
"cms_table"
