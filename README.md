# terralink

Standalone R implementation of the TerraLink corridor optimization workflows (raster + vector).

This package mirrors the QGIS TerraLink plugin capabilities while keeping the workflow scriptable and clear for R users.

Raster inputs are polygonized and analyzed through the vector corridor pipeline, matching the current QGIS plugin workflow.

Canonical TerraLink 1.7 strategies:

- `most_connected_habitat`
- `largest_single_network`
- `reachable_habitat_advanced`
- `landscape_fluidity`

## What TerraLink is for

TerraLink is designed for scenario-based corridor planning rather than one-shot "best answer" mapping.

A typical workflow is:

1. load habitat data
2. choose a strategy that matches the ecological outcome you care about
3. set approximate planning constraints such as budget, corridor width, and search distance
4. run a few plausible scenarios
5. compare the PRE/POST metrics and mapped outputs

This is the same framing used in the QGIS plugin: TerraLink is most useful when you compare alternatives and decide which tradeoff is best for the landscape, not when you expect a single universal score to answer the problem.

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

paths <- terralink_sample_data()
paths
# $raster   -> .../extdata/habitat.tif
# $vector   -> .../extdata/patches.gpkg
# $obstacle -> .../extdata/impassable.gpkg

# Raster workflow
res_r <- terralink_raster(
  raster = paths["raster"],
  patch_values = 1,
  obstacle_values = 2,
  budget = 220,
  min_patch_size = 10,
  min_corridor_width = 3,
  max_search_distance = 30,
  corridor_cell_assignment = "sum_total_network_area",
  units = "pixels",
  strategy = "most_connected_habitat"
)
plot(res_r)

# Vector workflow
res_v <- terralink_vector(
  patches = paths["vector"],
  budget = 2.5,
  min_patch_size = 0.1,
  min_corridor_width = 20,
  max_search_distance = 500,
  units = "metric",
  strategy = "most_connected_habitat"
)
plot(res_v)
```

## Optimization modes

The TerraLink 1.7 strategies are not interchangeable labels. Each one rewards a different kind of connectivity improvement.

- `largest_single_network`: build one dominant connected backbone.
- `most_connected_habitat`: maximize structurally connected habitat across the whole map.
- `landscape_fluidity`: improve movement quality, alternative routing, and bottleneck relief.
- `reachable_habitat_advanced`: maximize functionally reachable habitat under a dispersal rule or kernel.

In practice:

- use `largest_single_network` when you want one strong consolidated system
- use `most_connected_habitat` when you want broad structural integration
- use `landscape_fluidity` when movement redundancy and detour relief matter
- use `reachable_habitat_advanced` when you have a focal species or explicit dispersal assumptions

That distinction matters because structural metrics and functional accessibility metrics capture different ecological processes, so mode choice should follow the planning question rather than habit.

## Scripts

Packaged scripts are discoverable from R:

```r
terralink_examples()
terralink_examples("raster")
terralink_examples("vector")
terralink_sample_data()
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

Raster runs are delegated through the vector engine after polygonizing habitat and impassable areas, so raster and vector workflows share the same corridor-selection core.

Vector:

- `corridors` (sf): buffered corridor polygons with metrics.
- `networks` (sf): dissolved connected networks.
- CSV summary + metrics report text.

## Interpreting metrics

The package writes PRE/POST metrics so you can compare scenarios on the same landscape.

The most useful habit is comparative interpretation:

- compare scenario A vs scenario B on the same data
- compare PRE vs POST for the same scenario
- avoid treating one metric value as a context-free "true" connectivity score

Some metrics are mainly structural, such as connected habitat area and largest network area. Others are closer to functional accessibility or whole-network movement behavior, such as habitat availability and fluidity-oriented summaries. That distinction is important when deciding whether your goal is consolidation, broad structural gain, or species-relevant reachability.

Raster and vector workflows can produce similar directional results without being numerically identical. They use different spatial representations and corridor construction steps, so the safest comparison is always scenario-vs-scenario within the same pipeline.

## Obstacle-aware routing (optional)

If you want corridors to route around impassable features, install the optional packages and enable obstacles:

```r
# Optional (for shortest-path routing around obstacles)
install.packages(c("gdistance", "raster", "sp"))

result <- terralink_vector(
  patches = terralink_sample_data("vector"),
  budget = 2.5,
  min_corridor_width = 20,
  max_search_distance = 500,
  units = "metric",
  obstacle_layers = terralink_sample_data("obstacle"),
  obstacle_resolution = 10
)
```

Raster mode supports impassable values/ranges directly from the raster and funnels them through the same vector corridor engine:

```r
result <- terralink_raster(
  raster = terralink_sample_data("raster"),
  patch_values = 1,
  obstacle_values = 2,
  budget = 220,
  units = "pixels"
)
```
