library(ggplot2)

root <- normalizePath(".", mustWork = TRUE)
expressions <- parse(file.path(root, "script", "3_plots.R"))
definition <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("island_satellite"))
}, expressions)
stopifnot(length(definition) == 1L)
eval(definition[[1]])

local({
  work <- tempfile("island-satellite-test-")
  cache <- file.path(work, "data", ".eox_2025_cache")
  dir.create(cache, recursive = TRUE)
  original <- getwd()
  on.exit({
    setwd(original)
    unlink(file.path(cache, c("lon179_lat28.tile", "lon180_lat28.tile",
                              "lon181_lat28.tile")))
    unlink(cache)
    unlink(file.path(work, "data"))
    unlink(work)
  })

  image <- array(0, dim = c(2048, 2048, 4))
  image[1:1024, , 1] <- 1
  image[, 1:1024, 2] <- 1
  image[, , 4] <- 1
  image[900, 200, 4] <- 0
  png::writePNG(image, file.path(cache, "lon180_lat28.tile"))
  jpeg::writeJPEG(image[, , 1:3], file.path(cache, "lon179_lat28.tile"))
  jpeg::writeJPEG(image[seq(1, 2048, 8), seq(1, 2048, 8), 1:3],
                  file.path(cache, "lon181_lat28.tile"))
  setwd(work)

  warnings <- character()
  layers <- withCallingHandlers(
    island_satellite(c(179.75, 180.25), c(28.25, 28.625)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(length(layers) == 2L,
            any(grepl("missing coverage", warnings, fixed = TRUE)))

  raster <- layers[[2]]$geom_params$raster
  stopifnot(inherits(raster, "nativeRaster"),
            identical(dim(raster), c(768L, 512L)))
  original_image <- png::readPNG(file.path(cache, "lon180_lat28.tile"), native = TRUE)
  pixel <- function(x, row, col) as.vector(x)[(row - 1L) * ncol(x) + col]
  stopifnot(identical(pixel(raster, 100, 100), pixel(original_image, 868, 100)),
            identical(pixel(raster, 500, 300), pixel(original_image, 1268, 300)))
  warnings <- character()
  low_resolution <- withCallingHandlers(
    island_satellite(c(181.25, 181.75), c(28.25, 28.625)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(identical(dim(low_resolution[[1]]$geom_params$raster), c(96L, 128L)),
            any(grepl("lower resolution", warnings, fixed = TRUE)))
  stopifnot(inherits(try(island_satellite(c(180, 180), c(28, 29)), silent = TRUE),
                     "try-error"))
  cat("Island satellite tile cropping and dateline checks passed.\n")
})
