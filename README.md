# terralink

Standalone R implementation of the TerraLink corridor optimization workflows (raster + vector).

This package mirrors the QGIS TerraLink plugin capabilities while keeping the workflow scriptable and clear for R users.

## Project links

- TerraLink online tool: https://sorusconsultingllc.com/tools/terralink-online
- QGIS plugin source: https://github.com/sorus-tools/TerraLink
- Contact: info@sorusconsultingllc.com

## Install

```r
install.packages("remotes")
remotes::install_local("/path/to/terralink")
library(terralink)
```

## Quick start

```r
library(terralink)

# Raster workflow
res_r <- terralink_raster(
  raster = "habitat.tif",
  patch_values = 1,
  budget = 50,
  min_patch_size = 2,
  min_corridor_width = 3,
  max_search_distance = 600,
  units = "pixels",
  strategy = "circuit_utility"
)
plot(res_r)

# Vector workflow
res_v <- terralink_vector(
  patches = "patches.gpkg",
  budget = 1,
  min_patch_size = 1,
  min_corridor_width = 20,
  max_search_distance = 5000,
  units = "metric",
  strategy = "circuit_utility"
)
plot(res_v)
```

## Scripts

Packaged scripts are discoverable from R:

```r
terralink_examples()
terralink_examples("raster")
terralink_examples("vector")
```

Main scripts:

- `inst/scripts/run_raster.R`
- `inst/scripts/run_vector.R`
- `inst/scripts/example_raster_barbados.R`
- `inst/scripts/example_vector_barbados.R`

The two `example_*_barbados.R` scripts use the exact parameter sets validated during raster/vector parity testing.

## Outputs

Raster:

- `corridor_raster`: corridor cells labeled by connected area score.
- `contiguous_raster`: habitat + corridors labeled by component size.
- `patches`: patch labels.
- CSV summaries + metrics report text.

Vector:

- `corridors` (sf): buffered corridor polygons with metrics.
- `networks` (sf): dissolved connected networks.
- CSV summary + metrics report text.

## Obstacle-aware routing (optional)

If you want corridors to route around impassable features, install the optional packages and enable obstacles:

```r
# Optional (for shortest-path routing around obstacles)
install.packages(c("gdistance", "raster", "sp"))

result <- terralink_vector(
  patches = "patches.gpkg",
  budget = 50,
  min_corridor_width = 200,
  max_search_distance = 5000,
  units = "metric",
  obstacle_layers = "impassable.gpkg",
  obstacle_resolution = 50
)
```

Raster mode supports impassable values/ranges without extra packages; the optional packages improve routing quality:

```r
result <- terralink_raster(
  raster = "habitat.tif",
  patch_values = 1,
  obstacle_values = c(9, 10),
  budget = 500,
  units = "pixels"
)
```
