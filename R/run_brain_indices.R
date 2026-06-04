#!/usr/bin/env Rscript

# Load required libraries silently
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(fs)
  library(dplyr)
  library(geometry) # for convhulln
  library(entropy)  # for grid entropy
  library(fields)   # for rdist
})

# ==========================================
# 1. Catch Arguments from Bash Shell Script
# ==========================================
option_list <- list(
  make_option(c("-o", "--output_dir"), type = "character", default = "output", 
              help = "Path to the main output directory (Default: output)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ------------------------------------------
# 🛠️ RSTUDIO TESTING AREA 🛠️
# UNCOMMENT ONLY THE ONE LINE BELOW to test in RStudio. 
#
# opt <- list(output_dir = "output")
# ------------------------------------------

# Set up paths based on the main output directory provided by Bash
output_dir <- opt$output_dir
input_mds_dir <- path(output_dir, "MDSpoint")

if (!dir_exists(input_mds_dir)) {
  stop(sprintf("MDS directory not found at: %s. Run MDS analysis first.", input_mds_dir), call. = FALSE)
}

# ==========================================
# 2. Define Helper Functions
# ==========================================
get_convex_hull_volume <- function(coords_mat) {
  coords_mat <- as.matrix(coords_mat)
  n_dim <- ncol(coords_mat)
  
  tryCatch({
    # FA option computes generalized area/volume
    hull <- convhulln(coords_mat, options = "FA")
    raw_vol <- hull$vol
    # Normalization: D-th root converts hypervolume to linear equivalent
    norm_vol <- raw_vol^(1/n_dim)
    return(list(raw = raw_vol, norm = norm_vol))
  }, error = function(e) {
    # Failsafe if points are perfectly coplanar/degenerate
    return(list(raw = NA, norm = NA))
  })
}

get_grid_entropy <- function(coords_mat, bins_per_dim = 5) {
  n_points <- nrow(coords_mat)
  n_dim <- ncol(coords_mat) # Extract the dimension count
  
  # Normalize coordinates to 0-1 range
  norm_coords <- apply(coords_mat, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Bin the data
  binned_data <- apply(norm_coords, 2, function(x) {
    cut(x, breaks = bins_per_dim, labels = FALSE, include.lowest = TRUE)
  })
  
  # Identify unique grid states
  grid_states <- apply(binned_data, 1, paste, collapse = "_")
  counts <- table(grid_states)
  
  # Calculate raw Shannon entropy
  raw_entropy <- entropy(counts, unit = "log")
  
  # NORMALIZATION: Max entropy is bounded by either the number of TRs 
  # or the total number of available bins in that specific D-dimensional space.
  total_possible_bins <- bins_per_dim^n_dim
  max_possible_states <- min(n_points, total_possible_bins)
  max_entropy <- log(max_possible_states)
  
  # Calculate normalized entropy (0 to 1)
  norm_entropy <- raw_entropy / max_entropy
  
  return(list(raw = raw_entropy, norm = norm_entropy))
}

get_laminarity <- function(coords_mat, target_RR = 0.05, min_vert_line = 2) {
  coords_mat <- scale(coords_mat) 
  n_points <- nrow(coords_mat)
  dist_mat <- rdist(coords_mat)
  
  dist_values <- dist_mat[lower.tri(dist_mat)]
  radius <- quantile(dist_values, probs = target_RR)
  
  rec_mat <- ifelse(dist_mat <= radius, 1, 0)
  total_recurrent_points <- sum(rec_mat) - n_points 
  vertical_line_points <- 0
  
  for (j in 1:ncol(rec_mat)) {
    rle_col <- rle(rec_mat[, j])
    vertical_lines <- rle_col$lengths[rle_col$values == 1]
    valid_lines <- vertical_lines[vertical_lines >= min_vert_line]
    vertical_line_points <- vertical_line_points + sum(valid_lines)
  }
  
  if (total_recurrent_points == 0) return(0)
  return(vertical_line_points / total_recurrent_points)
}

# ==========================================
# 3. File Setup & Initialization
# ==========================================
file_list <- dir_ls(input_mds_dir, regexp = "_MDSpoints.csv$")

if (length(file_list) == 0) {
  stop("No MDS point files found. Make sure step 1 completed successfully.", call. = FALSE)
}

results <- data.frame(
  Subject  = character(length(file_list)),
  CHA_Raw  = numeric(length(file_list)),
  CHA_Norm = numeric(length(file_list)),
  GE_Raw   = numeric(length(file_list)),
  GE_Norm  = numeric(length(file_list)),
  LAM_Norm = numeric(length(file_list)),
  stringsAsFactors = FALSE
)

# ==========================================
# 4. Process Each Subject
# ==========================================
cat("\n=================================================\n")
cat(sprintf("Calculating Brain Indices (CHA, GE, LAM) for %d subjects...\n", length(file_list)))
cat("=================================================\n")

for (i in seq_along(file_list)) {
  df <- read_csv(file_list[i], show_col_types = FALSE)
  
  # Clean subject name robustly
  subject_name <- path_file(file_list[i]) |> sub("_optimal_MDSpoints.csv|_MDSpoints.csv", "", x = _)
  
  # Execute functions
  cha_res <- get_convex_hull_volume(df)
  ge_res  <- get_grid_entropy(df)
  lam_val <- get_laminarity(df)
  
  # Assign to the dataframe
  results$Subject[i]  <- subject_name
  results$CHA_Raw[i]  <- cha_res$raw
  results$CHA_Norm[i] <- cha_res$norm
  results$GE_Raw[i]   <- ge_res$raw
  results$GE_Norm[i]  <- ge_res$norm
  results$LAM_Norm[i] <- lam_val
  
  cat(sprintf("  Processed: %s | CHA: %.3f | GE: %.3f | LAM: %.3f\n", 
              subject_name, cha_res$norm, ge_res$norm, lam_val))
}

# ==========================================
# 5. Save Results
# ==========================================
output_file <- path(output_dir, "Brain_Indices_Summary.csv")
write_csv(results, file = output_file)

cat("\n✅ Success! All indices calculated and saved to:", output_file, "\n")