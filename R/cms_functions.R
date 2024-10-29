
# make_cms_envir ----------------------------------------------------------

# Create a python environment with the copernicusmarine toolbox installed

#' Create a Python vertual environment with copernicusmarine toolbox installed
#'
#' @param my_env The name of, or path to, a Python virtual environment. If this
#' name contains any slashes, the name will be interpreted as a path; if the name
#' does not contain slashes, it will be treated as a virtual environment within virtualenv_root().
#' When NULL, the virtual environment as specified by the RETICULATE_PYTHON_ENV
#' environment variable will be used instead. To refer to a virtual environment
#' in the current working directory, you can prefix the path with ⁠./<name>⁠.
#' @param user User name for copernicus marine service. String.
#' @param password Password for copernicus marine service. String.
#' @param force Force recreating the environment specified by envname, even if
#' it already exists. If TRUE, the pre-existing environment is first deleted and
#' then recreated. Otherwise, if FALSE (the default), the path to the existing environment is returned.
#'
#' @details This function uses reticulate to create a Python virtual environment
#' and install the copernicusmarine toolbox. This toolbox can be used to download
#' oceanographic data from copernicus marine service.
#'
#' If the function is run including the user name and password then login details
#' will be saved to the virtual environment on creation.
#'
#' @return Creates a python virtual environment with the copernicusmarine toolbox installed on the path specified by my_env
#' @export
#'
#' @examplesIf interactive()
#'
#' cm_env <- 'tmp/cms_env'
#' make_cms_envir(my_env = cm_env, force = FALSE)
#' reticulate::use_virtualenv(paste0(my_env,'/'), required = TRUE)
#' cm <- reticulate::import("copernicusmarine")
#' cm$login()
#'
make_cms_envir <- function(my_env, force = FALSE, user = NULL, password = NULL) {
  if (reticulate::virtualenv_exists(my_env) == F) {
    reticulate::install_python(version = "3.10:latest", list = FALSE, force = TRUE, optimized = TRUE)
    reticulate::virtualenv_create(envname = my_env,packages = c("numpy", 'copernicusmarine'),
                                  ignore_installed = FALSE,
                                  force = FALSE)
    reticulate::virtualenv_install(envname = my_env, packages = c("numpy", "copernicusmarine"))

  } else message(paste('Python environment already exists at:', my_env))

  if (!is.null(user) & !is.null(password)) {
    reticulate::use_virtualenv(paste0(my_env,'/'), required = TRUE)
    cm <- reticulate::import("copernicusmarine")
    cm$login(user, password, skip_if_user_logged_in = T)
  }
}

