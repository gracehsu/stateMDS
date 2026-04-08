# load libraries
library(vegan)
library(ggplot2)
library(stats)
library(readr)
library(here)
library(fs)

# 1. Define input/output folders using relative paths
input_dir       <- here("data", "voxels")
output_arrowdis <- here("output", "arrowdis")
output_mdspoint <- here("output", "MDSpoint")

# Create output directories if they don’t exist
dir_create(output_arrowdis)
dir_create(output_mdspoint)

# 2. Configuration
file_names <- dir_ls(input_dir, regexp = "\\.csv$")
dimension  <- 2    # MDS dimensions
max_TR     <- 240  # Maximum TRs to analyze per subject

# 3. Initialize Summary Table
# Pre-defining types prevents numeric data from being forced into characters
m <- data.frame(
  Subject       = character(length(file_names)),
  Velocity_Total = numeric(length(file_names)),
  Velocity_Mean  = numeric(length(file_names)),
  Stress        = numeric(length(file_names)),
  Converged     = logical(length(file_names)),
  stringsAsFactors = FALSE
)

# 4. Main loop
for (i in seq_along(file_names)) {
  
  # Load data
  raw_data <- read_csv(file_names[i], col_names = FALSE, show_col_types = FALSE)
  
  # Ensure we don't exceed available TRs
  actual_TR <- min(nrow(raw_data), max_TR)
  DMN       <- t(raw_data[1:actual_TR, ])
  DMN_XX    <- crossprod(DMN)
  
  # Standardize the Gram matrix
  DMN_XX_standardized <- decostand(DMN_XX, method = "standardize")
  
  # MDS Analysis
  # Note: trymax=50 is good; added trace=0 to keep the console clean
  mds_result <- metaMDS(DMN_XX_standardized, autotransform = FALSE, 
                        k = dimension, distance = "euclidean", 
                        trymax = 50, trace = 0)
  
  mds_points <- mds_result[["points"]]
  
  # Robust Subject Naming
  subject_name <- path_file(file_names[i]) |> path_ext_remove() |> sub("_voxels", "", x = _)
  
  # Save MDS points
  write_csv(as.data.frame(mds_points), 
            file = path(output_mdspoint, paste0(subject_name, "_MDSpoints.csv")))
  
  # 5. Vectorized Arrow Distances (Faster than a k-loop)
  # Calculates Euclidean distance between row i and row i+1
  diffs <- diff(mds_points) # Matrix of differences between consecutive TRs
  distances <- sqrt(rowSums(diffs^2))
  
  # Save arrow distances
  write_csv(data.frame(TR_Transition = 1:length(distances), Distance = distances), 
            file = path(output_arrowdis, paste0(subject_name, "_arrow_distance.csv")))
  
  # 6. Record Summary Metrics
  m$Subject[i]        <- subject_name
  m$Velocity_Total[i] <- sum(distances)
  m$Velocity_Mean[i]  <- mean(distances)
  m$Stress[i]         <- mds_result[["stress"]]
  m$Converged[i]      <- mds_result[["converged"]]
  
  message(paste("Successfully processed:", subject_name))
}

# 7. Save summary results
write_csv(m, file = path(output_arrowdis, "MDS_Velocity_Summary.csv"))
message("Analysis Complete. Results saved in output/ folder.")