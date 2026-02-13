# ═══════════════════════════════════════════════════════════════════════════
# 🌲 TerraLink Analysis Runner
# ═══════════════════════════════════════════════════════════════════════════
# HOW TO USE:
# 1. Save this script as "terralink_runner.R"
# 2. In RStudio console, type: source("path/to/terralink_runner.R")
# 3. Then type either: run_raster() or run_vector()
# 4. A file picker will appear
# ═══════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(terralink)
  library(terra)
  library(sf)
})


# ═══════════════════════════════════════════════════════════════════════════
# ⚙️  DEFAULT PARAMETERS (edit these if needed)
# ═══════════════════════════════════════════════════════════════════════════

.RASTER_PARAMS <- list(
  patch_values          = 1,
  patch_ranges          = NULL,
  obstacle_values       = NULL,
  obstacle_ranges       = NULL,
  budget                = 500,
  min_patch_size        = 20,
  min_corridor_width    = 3,
  max_search_distance   = 300,
  units                 = "pixels",
  allow_bottlenecks     = FALSE,
  strategy              = "circuit_utility",
  verbose               = 1
)

.VECTOR_PARAMS <- list(
  budget                = 5,
  min_patch_size        = 1.0,
  min_corridor_width    = 20,
  max_search_distance   = 5000,
  units                 = "metric",
  return_crs            = "input",
  obstacle_resolution   = NULL,
  strategy              = "circuit_utility",
  verbose               = 1
)


# ═══════════════════════════════════════════════════════════════════════════
# 🚀 MAIN FUNCTIONS - Call these from the console
# ═══════════════════════════════════════════════════════════════════════════

run_raster <- function(params = .RASTER_PARAMS) {
  .banner("🗺️  RASTER ANALYSIS")

  # Pick input file
  .msg("Opening file picker...")
  input_file <- .pick_file("Select Raster Input (.tif, .img, etc.)")

  if (is.null(input_file)) {
    .msg("❌ No file selected. Cancelled.")
    return(invisible(NULL))
  }

  output_dir <- dirname(input_file)

  .msg("Input:  %s", basename(input_file))
  .msg("Output: %s", output_dir)
  .msg("Running analysis...")

  # Run analysis
  result <- do.call(terralink_raster, c(
    list(raster = input_file, progress = FALSE),
    params
  ))

  # Display
  print(result)
  .msg("Plotting results...")
  plot(result)

  # Save
  .msg("Saving outputs...")
  prefix <- sprintf("%s_terralink", tools::file_path_sans_ext(basename(input_file)))
  outputs <- write_terralink_raster_outputs(
    result, output_dir, prefix,
    overwrite = TRUE,
    output_paths = list()
  )

  .show_results(outputs)
  invisible(list(result = result, outputs = outputs))
}

run_vector <- function(params = .VECTOR_PARAMS, ask_for_obstacles = FALSE) {
  .banner("📍 VECTOR ANALYSIS")

  # Pick input file
  .msg("Opening file picker...")
  input_file <- .pick_file("Select Vector Patch File (.shp, .gpkg, etc.)")

  if (is.null(input_file)) {
    .msg("❌ No file selected. Cancelled.")
    return(invisible(NULL))
  }

  output_dir <- dirname(input_file)
  .msg("Patches: %s", basename(input_file))

  # Handle optional obstacles
  obstacle_file <- NULL
  if (ask_for_obstacles) {
    .msg("Opening file picker for obstacles (cancel to skip)...")
    obstacle_file <- .pick_file("Select Obstacle Layer (Optional)")
    if (!is.null(obstacle_file)) {
      .msg("Obstacles: %s", basename(obstacle_file))
    }
  }

  .msg("Output: %s", output_dir)
  .msg("Running analysis...")

  # Build arguments
  args <- c(
    list(patches = input_file, progress = FALSE),
    params
  )
  if (!is.null(obstacle_file)) {
    args$obstacle_layers <- obstacle_file
  }

  # Run analysis
  result <- do.call(terralink_vector, args)

  # Display
  print(result)
  .msg("Plotting results...")
  plot(result)

  # Save
  .msg("Saving outputs...")
  prefix <- sprintf("%s_terralink", tools::file_path_sans_ext(basename(input_file)))
  outputs <- write_terralink_vector_outputs(
    result, output_dir, prefix,
    overwrite = TRUE,
    output_paths = list()
  )

  .show_results(outputs)
  invisible(list(result = result, outputs = outputs))
}


# ═══════════════════════════════════════════════════════════════════════════
# 🛠️  HELPER FUNCTIONS (Internal use only)
# ═══════════════════════════════════════════════════════════════════════════

.pick_file <- function(caption = "Select a file") {
  # Try RStudio API
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    result <- tryCatch({
      path <- rstudioapi::selectFile(caption = caption, existing = TRUE)
      if (is.character(path) && length(path) > 0 && nzchar(path)) {
        return(normalizePath(path, mustWork = TRUE))
      }
      return(NULL)
    }, error = function(e) NULL)

    if (!is.null(result)) return(result)
  }

  # Fallback to base R file picker
  result <- tryCatch({
    path <- file.choose()
    if (is.character(path) && length(path) > 0 && nzchar(path) && file.exists(path)) {
      return(normalizePath(path, mustWork = TRUE))
    }
    return(NULL)
  }, error = function(e) NULL)

  return(result)
}

.banner <- function(text) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  ", text, "\n", sep = "")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("\n")
}

.msg <- function(fmt, ...) {
  cat(sprintf(paste0("  ", fmt, "\n"), ...))
}

.show_results <- function(outputs) {
  cat("\n")
  cat("  ✅ Analysis complete!\n\n")
  cat("  📦 Generated files:\n")
  for (name in names(outputs)) {
    cat(sprintf("     • %s: %s\n", name, basename(outputs[[name]])))
  }
  cat("\n")
}


# ═══════════════════════════════════════════════════════════════════════════
# 📋 INSTRUCTIONS
# ═══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  TerraLink Runner Loaded Successfully!\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")
cat("  Available commands:\n")
cat("    • run_raster()  - Run raster analysis\n")
cat("    • run_vector()  - Run vector analysis\n")
cat("    • run_vector(ask_for_obstacles = TRUE)  - Include obstacles\n")
cat("\n")
cat("  Just type one of these commands in the console and press Enter!\n")
cat("\n")
