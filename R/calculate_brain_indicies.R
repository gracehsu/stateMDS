# load libraries
library(readr)
library(here)
library(fs)
library(dplyr)
library(geometry) # for convhulln
library(entropy)  # for grid entropy
library(fields)   # for rdist

# 1. --- Define Helper Functions ---

get_convex_hull_volume <- function(coords_mat) {
  coords_mat <- as.matrix(coords_mat)
  n_dim <- ncol(coords_mat)
  
  tryCatch({
    hull <- convhulln(coords_mat, options = "FA")
    raw_vol <- hull$vol
    # Normalization: D-th root converts hypervolume to linear equivalent
    norm_vol <- raw_vol^(1/n_dim)
    return(list(raw = raw_vol, norm = norm_vol))
  }, error = function(e) {
    return(list(raw = NA, norm = NA))
  })
}

get_grid_entropy <- function(coords_mat, bins_per_dim = 5) {
  n_points <- nrow(coords_mat)
  norm_coords <- apply(coords_mat, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  binned_data <- apply(norm_coords, 2, function(x) {
    cut(x, breaks = bins_per_dim, labels = FALSE, include.lowest = TRUE)
  })
  grid_states <- apply(binned_data, 1, paste, collapse = "_")
  counts <- table(grid_states)
  
  raw_entropy <- entropy(counts, unit = "log")
  # Normalization: relative to max possible entropy log(N)
  norm_entropy <- raw_entropy / log(n_points)
  
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
  # LAM is naturally a ratio between 0 and 1
  return(vertical_line_points / total_recurrent_points)
}

# 2. --- File Setup ---
input_dir <- here("output", "MDSpoint")
output_dir <- here("output")
file_list <- dir_ls(input_dir, regexp = "_MDSpoints.csv$")

# Initialize results table with your specific column requirements
results <- data.frame(
  Subject  = character(length(file_list)),
  CHA_Raw  = numeric(length(file_list)),
  CHA_Norm = numeric(length(file_list)),
  GE_Raw   = numeric(length(file_list)),
  GE_Norm  = numeric(length(file_list)),
  LAM_Norm = numeric(length(file_list)),
  stringsAsFactors = FALSE
)

# 3. --- Process Each Subject ---
message("Calculating Brain Indices (CHA, GE, LAM)...")

for (i in seq_along(file_list)) {
  df <- read_csv(file_list[i], show_col_types = FALSE)
  subject_name <- path_file(file_list[i]) |> sub("_MDSpoints.csv", "", x = _)
  
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
  
  message(paste("Processed:", subject_name))
}

# 4. --- Save Results ---
# Saving as CSV for better GitHub/R compatibility
write_csv(results, file = path(output_dir, "Brain_Indices_Summary.csv"))
message("Success! Summary saved to output/Brain_Indices_Summary.csv")