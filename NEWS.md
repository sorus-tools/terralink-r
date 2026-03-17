# terralink 1.7.0

- Updated the R package to TerraLink 1.7 strategy names and aliases:
  `most_connected_habitat`, `largest_single_network`,
  `reachable_habitat_advanced`, and `landscape_fluidity`.
- Added graph-based habitat availability and landscape fluidity selectors for
  raster and vector workflows, plus expanded pre/post connectivity reporting.
- Improved vector `most_connected_habitat` selection to stay closer to TerraLink
  1.7 QGIS results on the packaged parity fixture.
- Added packaged synthetic raster/vector parity fixtures and an
  synthetic parity fixture set validated against QGIS.
- Expanded README/vignette guidance with optimization-mode selection,
  scenario-comparison advice, and raster-vs-vector interpretation notes taken
  from the TerraLink 1.7 QGIS documentation.
- Added methods references in `DESCRIPTION` and removed the interactive wrapper
  from `terralink_raster()` examples for CRAN submission readiness.

# terralink 1.5.6

- Added an explicit `sf` availability check in `inst/scripts/run_vector.R`
  before calling `sf::st_read()`.
- Reduced `terralink_raster()` example workload to stay below CRAN example
  timing thresholds on slower builders.

- Added packaged sample data (`inst/extdata`) plus `terralink_sample_data()` for
  copy/paste runnable examples.
- Updated vignette, README, and example scripts to use packaged sample files.
- Aligned non-debug behavior with TerraLink 1.5.3 strategy/selector updates:
  `most_connectivity` defaults, enhanced overlap/distance guards, raster
  corridor assignment modes, and vector local efficiency scoring.

# terralink 1.5.1

- CRAN submission cleanup: metadata, URLs, and package build exclusions.
- Resubmission release to address incoming pretest notes.

# terralink 1.5.0

- First official R release matching TerraLink 1.5 features.
- Raster + vector workflows, metrics, and output writers.
- Optional obstacle-aware routing via gdistance.
