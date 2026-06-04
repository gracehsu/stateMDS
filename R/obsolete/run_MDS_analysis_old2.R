#!/usr/bin/env Rscript

# Load required libraries silently
suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(readr)
  library(fs)
  library(dplyr)
})

# ==========================================
# 1. Catch Arguments from Bash Shell Script
# ==========================================
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = "data/voxels", 
              help = "Path to the directory containing TSV/CSV matrix files (Default: data/voxels)"),
  make_option(c("-o", "--output_dir"), type = "character", default = "output", 
              help = "Path to the output directory (Default: output)"),
  make_option(c("-t", "--max_tr"), type = "integer", default = 180, 
              help = "Maximum number of TRs to analyze per subject (Default: 180)"),
  make_option(c("-s", "--stress"), type = "numeric", default = 0.15, 
              help = "Maximum acceptable stress value for dimension selection (Default: 0.15)"),
  make_option(c("-k", "--max_dim"), type = "integer", default = 10, 
              help = "Maximum dimension to test before giving up (Default: 10)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ------------------------------------------
# 🛠️ RSTUDIO TESTING AREA 🛠️
# UNCOMMENT ONLY THE ONE LINE BELOW to test in RStudio. 
#
# opt <- list(input_dir = "data/voxels", output_dir = "output", max_tr = 180, stress = 0.15, max_dim = 10)
# ------------------------------------------

input_dir <- opt$input_dir
output_dir <- opt$output_dir
user_max_tr <- opt$max_tr
user_stress <- opt$stress
user_max_dim <- opt$max_dim

if (!dir_exists(input_dir)) {
  stop(sprintf("Input directory does not exist: %s", input_dir), call. = FALSE)
}

# ==========================================
# 2. Setup Output Directories
# ==========================================
out_mds <- path(output_dir, "MDSpoint")
out_vel <- path(output_dir, "arrowdis")
dir_create(out_mds)
dir_create(out_vel)

# ==========================================
# 3. Find All Files
# ==========================================
target_files <- dir_ls(input_dir, regexp = "\\.(csv|tsv)$")

if (length(target_files) == 0) {
  stop("No CSV or TSV files found in the input directory.", call. = FALSE)
}

cat(sprintf("Found %d files in %s. Starting batch processing...\n", length(target_files), input_dir))

summary_results <- list()

# ==========================================
# 4. Main Batch Loop
# ==========================================
for (file_path in target_files) {
  
  subject_name <- path_file(file_path) |> 
    path_ext_remove() |> 
    sub("_timeseries.*|_voxels.*", "", x = _)
  
  cat("\n-------------------------------------------------\n")
  cat("Processing Subject:", subject_name, "\n")
  
  raw_data <- read_delim(file_path, col_names = FALSE, show_col_types = FALSE)
  
  max_TR <- min(nrow(raw_data), user_max_tr)
  DMN <- t(raw_data[1:max_TR, ])
  DMN_XX <- crossprod(DMN)
  DMN_XX_standardized <- decostand(DMN_XX, method = "standardize")
  
  optimal_k <- NA
  optimal_stress <- NA
  optimal_points <- NULL
  
  cat(sprintf("Finding optimal dimension (stress < %.2f, testing up to k = %d)...\n", user_stress, user_max_dim))
  
  # The loop now dynamically searches from 2 up to your user_max_dim limit
  for (k in 2:user_max_dim) {
    
    # noshare=FALSE and wascores=FALSE are explicitly set to silence vegan's warnings
    mds_res <- suppressWarnings(
      metaMDS(DMN_XX_standardized, 
              autotransform = FALSE, 
              noshare = FALSE, 
              wascores = FALSE, 
              k = k, distance = "euclidean", trymax = 50, trace = 0)
    )
    
    stress_val <- mds_res$stress
    cat(sprintf("  k = %d | Stress: %.4f\n", k, stress_val))
    
    if (stress_val < user_stress) {
      cat(sprintf("  --> Success! Selected dimension k = %d\n", k))
      optimal_k <- k
      optimal_stress <- stress_val
      optimal_points <- mds_res$points
      break 
    }
  }
  
  # If it tests all the way up to max_dim and still fails, default to the max_dim
  if (is.null(optimal_points)) {
    cat(sprintf("  --> Warning: Stress never reached < %.2f. Defaulting to k = %d.\n", user_stress, user_max_dim))
    optimal_k <- user_max_dim
    optimal_stress <- stress_val
    optimal_points <- mds_res$points
  }
  
  # Calculate Velocity
  diffs <- diff(optimal_points) 
  distances <- sqrt(rowSums(diffs^2))
  
  total_velocity <- sum(distances)
  mean_velocity <- mean(distances)
  
  # Save Individual Subject Data
  write_csv(as.data.frame(optimal_points), 
            file = path(out_mds, paste0(subject_name, "_optimal_MDSpoints.csv")))
  
  write_csv(data.frame(TR_Transition = 1:length(distances), Distance = distances), 
            file = path(out_vel, paste0(subject_name, "_optimal_arrow_distance.csv")))
  
  # Store summary data in the list
  summary_results[[subject_name]] <- data.frame(
    Subject = subject_name,
    Selected_Dimension = optimal_k,
    Stress_Value = optimal_stress,
    Velocity_Total = total_velocity,
    Velocity_Mean = mean_velocity,
    stringsAsFactors = FALSE
  )
}

# ==========================================
# 5. Save Master Summary
# ==========================================
cat("\n=================================================\n")
cat("All subjects processed. Generating master summary...\n")

final_summary_df <- bind_rows(summary_results)

summary_file <- path(output_dir, "Master_Velocity_Summary.csv")
write_csv(final_summary_df, summary_file)

cat("✅ Master summary saved to:", summary_file, "\n")