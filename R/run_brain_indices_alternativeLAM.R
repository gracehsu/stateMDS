#!/usr/bin/env Rscript

# ==============================================================================
# BATCH CALCULATION: CHA, GE, and LAM Sweep from MDS Coordinates
# * MATHEMATICALLY RIGOROUS & PIPELINE READY *
# The "Parameter Sweep" (AUC) Method for LAM 
##You run it again at $0.02$, $0.03$, $0.04$, and $0.05$.You plot these 5 LAM scores on a curve 
##and calculate the Area Under the Curve (AUC), or simply take the Mean LAM.By doing this, 
# ==============================================================================

# Load required libraries silently
suppressPackageStartupMessages({
  library(optparse)
  library(geometry)
  library(entropy)
  library(fields)
  library(dplyr)
  library(readr)
  library(fs)
})

# ==========================================
# 1. Define Mathematical Functions
# ==========================================

# --- A. Convex Hull Volume (CHA) ---
get_convex_hull_volume <- function(coords_mat) {
  coords_mat <- as.matrix(coords_mat)
  n_dim <- ncol(coords_mat)
  
  tryCatch({
    hull <- convhulln(coords_mat, options = "FA")
    raw_vol <- hull$vol
    
    # NORMALIZATION: Take the D-th root of the hypervolume
    norm_vol <- raw_vol^(1/n_dim)
    
    return(list(raw = raw_vol, norm = norm_vol))
    
  }, error = function(e) {
    return(list(raw = NA, norm = NA))
  })
}

# --- B. Grid Entropy (GE) ---
get_grid_entropy <- function(coords_mat, bins_per_dim = 5) {
  n_points <- nrow(coords_mat)
  n_dim <- ncol(coords_mat) 
  
  norm_coords <- apply(coords_mat, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  binned_data <- apply(norm_coords, 2, function(x) {
    cut(x, breaks = bins_per_dim, labels = FALSE, include.lowest = TRUE)
  })
  grid_states <- apply(binned_data, 1, paste, collapse = "_")
  counts <- table(grid_states)
  
  raw_entropy <- entropy(counts, unit = "log")
  
  # NORMALIZATION FIX: Max entropy is bounded by the TRs OR the available bins
  total_possible_bins <- bins_per_dim^n_dim
  max_possible_states <- min(n_points, total_possible_bins)
  max_entropy <- log(max_possible_states)
  
  norm_entropy <- raw_entropy / max_entropy
  
  return(list(raw = raw_entropy, norm = norm_entropy))
}

# --- C. Laminarity (Single Threshold) ---
get_laminarity <- function(coords_mat, target_RR, min_vert_line = 2) {
  coords_mat <- scale(coords_mat) 
  n_points <- nrow(coords_mat)
  dist_mat <- rdist(coords_mat)
  
  dist_values <- dist_mat[lower.tri(dist_mat)]
  radius <- quantile(dist_values, probs = target_RR, na.rm = TRUE)
  
  rec_mat <- ifelse(dist_mat <= radius, 1, 0)
  total_recurrent_points <- sum(rec_mat) - n_points 
  vertical_line_points <- 0
  
  for (j in 1:ncol(rec_mat)) {
    rle_col <- rle(rec_mat[, j])
    vertical_lines <- rle_col$lengths[rle_col$values == 1]
    valid_lines <- vertical_lines[vertical_lines >= min_vert_line]
    vertical_line_points <- vertical_line_points + sum(valid_lines)
  }
  
  if (total_recurrent_points <= 0) return(0)
  return(vertical_line_points / total_recurrent_points)
}

# --- D. Laminarity Parameter Sweep (AUC) ---
get_laminarity_sweep <- function(coords_mat, min_vert_line = 2) {
  # We test across the physiologically relevant 1% to 5% density range
  rr_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  lam_scores <- c()
  
  for (target_RR in rr_thresholds) {
    current_lam <- get_laminarity(coords_mat, target_RR, min_vert_line)
    lam_scores <- c(lam_scores, current_lam)
  }
  
  # Return the Area Under the Curve (Mean)
  return(mean(lam_scores, na.rm = TRUE))
}

# ==========================================
# 2. Catch Arguments from Bash Shell Script
# ==========================================
option_list <- list(
  make_option(c("-o", "--output_dir"), type = "character", default = "output", 
              help = "Path to the main output directory")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

output_dir <- opt$output_dir
input_mds_dir <- path(output_dir, "MDSpoint")

if (!dir_exists(input_mds_dir)) {
  stop(sprintf("MDS directory not found at: %s. Run MDS analysis first.", input_mds_dir), call. = FALSE)
}

file_list <- dir_ls(input_mds_dir, regexp = "_optimal_MDSpoints\\.csv$|_MDSpoints\\.csv$")

if (length(file_list) == 0) {
  stop("No MDS points files found.", call. = FALSE)
}

cat(sprintf("\n🔍 Calculating CHA, GE, and LAM (Sweep) for %d subjects...\n", length(file_list)))

# ==========================================
# 3. Main Processing Loop
# ==========================================
results_list <- list()

for (i in seq_along(file_list)) {
  file_path <- file_list[i]
  subject_name <- path_file(file_path) |> sub("_optimal_MDSpoints\\.csv|_MDSpoints\\.csv", "", x = _)
  
  df_points <- read_csv(file_path, show_col_types = FALSE)
  mds_matrix <- as.matrix(df_points)
  n_dim <- ncol(mds_matrix)
  
  # Calculate Indices
  cha_res <- get_convex_hull_volume(mds_matrix)
  ge_res  <- get_grid_entropy(mds_matrix, bins_per_dim = 5) 
  lam_auc <- get_laminarity_sweep(mds_matrix, min_vert_line = 2)
  
  # Store Results
  results_list[[i]] <- data.frame(
    Subject = subject_name,
    Dimensions = n_dim,
    CHA_Raw = cha_res$raw,
    CHA_Norm = cha_res$norm,
    GE_Raw = ge_res$raw,
    GE_Norm = ge_res$norm,
    LAM_Norm = lam_auc, # This is now the robust AUC score!
    stringsAsFactors = FALSE
  )
}

# ==========================================
# 4. Save Master Summary
# ==========================================
final_summary_df <- bind_rows(results_list)
summary_file <- path(output_dir, "Brain_Indices_Summary.csv")

write_csv(final_summary_df, summary_file)
cat("✅ Master Brain Indices summary saved to:", summary_file, "\n")