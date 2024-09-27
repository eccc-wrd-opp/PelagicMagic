


test_that("bathymetry files download", {

  suppressMessages(get_bathymetry(out_dir = 'tmp',
                 region = c(-68, -63, 44, 46),
                 overwrite = FALSE,
                 stride = 1,
                 vars = c('slope','coastdist','shelfdist'),
                 neighbors = 8,
                 slope_unit = c('degrees','radians'),
                 land_elev = 0,
                 dist_unit = 'km',
                 shelf_elev = -50))

  expect_equal(list.files('tmp'), c("bathymetry.nc", "coastdist.nc",  "shelfdist.nc",  "slope.nc"))
  unlink('tmp', recursive = T)
})

test_that("bathymetry files have not changed", {

  suppressMessages(get_bathymetry(out_dir = 'tmp',
                                  region = c(-68, -63, 44, 46),
                                  overwrite = FALSE,
                                  stride = 1,
                                  vars = NULL))

  r <- terra::rast('tmp/bathymetry.nc')
  expect_equal(terra::ncell(r), 577681)
  expect_equal(min(terra::values(r, na.rm = T)), -251)
  expect_equal(max(terra::values(r, na.rm = T)), 414)

  unlink('tmp', recursive = T)
})
