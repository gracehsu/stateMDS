# load libraries
library(vegan)
library(ggplot2)
library(stats)
library(readr)
library(here)
library(fs)

# Define input/output folders using relative paths
input_dir <- here("data", "voxels")
output_arrowdis <- here("output", "arrowdis")
output_mdspoint <- here("output", "MDSpoint")

# Create output directories if they don’t exist
dir_create(output_arrowdis)
dir_create(output_mdspoint)

# Get voxel CSV files
file_names <- dir_ls(input_dir, regexp = "*voxels.csv")

# Set MDS dimension
dimension <- 2 # change to 3 or higher as needed

#Set how many TRs for MDS analysis 
TR <- 240 # Adjust this to analyze fewer or more time points per subject

# Read voxel data for each subject
values <- lapply(file_names, read_csv, col_names = FALSE)

# Matrix to save velocity info
m <- data.frame(matrix(nrow = length(values), ncol = 5))
colnames(m) <- c("Subject", "Velocity_Total", "Velocity_Mean", "Stress", "Converged")

# Main loop
for (i in seq_along(file_names)) {
  DMN <- t(values[[i]][1:TR,])
  DMN_XX <- crossprod(DMN)
  
  # Standardize the Gram matrix
  DMN_XX_standardized <- decostand(DMN_XX, method = "standardize")
  
  # MDS
  mds_result <- metaMDS(DMN_XX_standardized, autotransform = FALSE, k = dimension, distance = "euclidean", trymax = 50)
  mds_points <- mds_result[["points"]]
  
  # Save MDS points
  subject_name <- sub("_voxels$", "", path_ext_remove(path_file(file_names[i])))
  write_csv(as.data.frame(mds_points), file = path(output_mdspoint, paste0(subject_name, "_MDSpoints.csv")))
  
  # Arrow distances
  n_steps <- nrow(DMN_XX) - 1
  distances <- numeric(n_steps)
  for (k in 1:n_steps) {
    distances[k] <- sqrt(sum((mds_points[k, ] - mds_points[k + 1, ])^2))
  }
  
  # Save arrow distances
  write_csv(data.frame(Distance = distances), file = path(output_arrowdis, paste0(subject_name, "_arrow_distance.csv")))
  
  # Summary metrics
  total_dist <- sum(distances)
  mean_dist <- mean(distances)
  stress <- mds_result[["stress"]]
  converged <- mds_result[["converged"]]
  
  m[i, ] <- c(subject_name, total_dist, mean_dist, stress, converged)
}

# Save summary results
write_csv(m, file = path(output_arrowdis, "MDS_Velocity_Summary.csv"))
