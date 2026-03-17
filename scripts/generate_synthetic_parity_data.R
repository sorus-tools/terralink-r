library(terra)
library(sf)

script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep('^--file=', commandArgs(trailingOnly = FALSE))], mustWork = FALSE)
script_path <- sub('^--file=', '', script_path)
root <- normalizePath(file.path(dirname(script_path), '..'), mustWork = TRUE)
out_dir <- file.path(root, 'inst', 'extdata')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mk_rect <- function(xmin, ymin, xmax, ymax) {
  st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax), c(xmin, ymin))))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Raster fixture: 35x35 pixels, five habitat patches, no impassables.
# This layout is deliberately chosen because TerraLink 1.7 in QGIS and the R
# package produce closely aligned selections across the 1.7 strategy set.
mat <- matrix(0, 35, 35)
mat[5:9, 5:9] <- 1
mat[5:9, 23:27] <- 1
mat[15:19, 14:18] <- 1
mat[25:29, 6:10] <- 1
mat[25:29, 24:28] <- 1
r <- rast(nrows = 35, ncols = 35, xmin = 0, xmax = 35, ymin = 0, ymax = 35, crs = 'EPSG:3857')
values(r) <- as.vector(t(mat))
writeRaster(r, file.path(out_dir, 'synthetic_habitat.tif'), overwrite = TRUE)

# Vector fixture: same general layout in metres plus a narrow obstacle.
patches <- st_sf(
  patch_name = c('north_west', 'north_east', 'midland', 'south_east'),
  quality = c(1.0, 0.8, 1.2, 1.0),
  geometry = st_sfc(
    mk_rect(100, 1600, 300, 1800),
    mk_rect(800, 1600, 1000, 1800),
    mk_rect(400, 900, 600, 1100),
    mk_rect(1100, 300, 1300, 500)
  ),
  crs = 3857
)
obstacles <- st_sf(
  obstacle = 'offsite_barrier',
  geometry = st_sfc(mk_rect(1700, 1700, 1800, 1800)),
  crs = 3857
)

patch_path <- file.path(out_dir, 'synthetic_patches.gpkg')
obs_path <- file.path(out_dir, 'synthetic_impassable.gpkg')
if (file.exists(patch_path)) file.remove(patch_path)
if (file.exists(obs_path)) file.remove(obs_path)
st_write(patches, patch_path, layer = 'patches', delete_layer = TRUE, quiet = TRUE)
st_write(obstacles, obs_path, layer = 'impassable', delete_layer = TRUE, quiet = TRUE)

cat('Wrote:\n')
cat(' -', file.path(out_dir, 'synthetic_habitat.tif'), '\n')
cat(' -', patch_path, '\n')
cat(' -', obs_path, '\n')
