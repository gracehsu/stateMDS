#!/usr/bin/env Rscript

# Load required libraries silently
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(fs)
  library(fields)
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
# opt <- list(output_dir = "output")
# ------------------------------------------

output_dir    <- opt$output_dir
input_mds_dir <- path(output_dir, "MDSpoint")
summary_file  <- path(output_dir, "Brain_Indices_Summary.csv")

# Setup Visualization Directories
dir_traj <- path(output_dir, "plots", "2d_trajectories")
dir_ge   <- path(output_dir, "plots", "ge_heatmaps")
dir_lam  <- path(output_dir, "plots", "LAM_plots")
dir_create(dir_traj)
dir_create(dir_ge)
dir_create(dir_lam)

# Verify inputs exist
if (!dir_exists(input_mds_dir)) stop("MDS directory not found. Run Step 1.", call. = FALSE)
if (!file_exists(summary_file)) stop("Brain_Indices_Summary.csv not found. Run Step 2.", call. = FALSE)

indices_df <- read_csv(summary_file, show_col_types = FALSE)
file_list  <- dir_ls(input_mds_dir, regexp = "_MDSpoints\\.csv$")

if (length(file_list) == 0) stop("No MDS points files found.", call. = FALSE)

# ==========================================
# 2. Define Plotting Functions
# ==========================================

# Function A: 2D Trajectory
plot_2d_trajectory <- function(df, subject_id) {
  df_2d <- df[, 1:2]
  colnames(df_2d) <- c("MDS1", "MDS2")
  df_2d$TR <- 1:nrow(df_2d)
  
  ggplot(df_2d, aes(x = MDS1, y = MDS2, color = TR)) +
    geom_path(alpha = 0.5, linewidth = 0.8) + 
    geom_point(size = 1.5) + 
    scale_color_viridis_c(option = "plasma") + 
    theme_minimal() +
    labs(title = paste("2D State Space Trajectory:", subject_id),
         subtitle = "Visualized across NMDS1 and NMDS2",
         x = "MDS Dimension 1", y = "MDS Dimension 2", color = "Time (TR)")
}

# Function B: GE Heatmap
plot_ge_heatmap <- function(df, subject_id, ge_raw, ge_norm, bins = 25) {
  df_2d <- df[, 1:2]
  colnames(df_2d) <- c("MDS1", "MDS2")
  
  h_raw_str  <- format(round(ge_raw, 4), nsmall = 4)
  h_norm_str <- format(round(ge_norm, 4), nsmall = 4)
  
  ggplot(df_2d, aes(x = MDS1, y = MDS2)) +
    geom_bin2d(bins = bins) + 
    scale_fill_gradientn(
      colors = c("#FFD700", "#FF8C00", "#FF4500", "#B22222", "#8B0000"),
      name = "TR Count", limits = c(1, NA)
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "#FFFFE0", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.2),
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste0("Subject ", subject_id, " (H_raw = ", h_raw_str, ", H_norm = ", h_norm_str, ")"),
         subtitle = "Grid Occupancy (NMDS1 vs NMDS2)", x = "NMDS1", y = "NMDS2")
}

# Function C: LAM Recurrence Plot (Updated for Subject-Specific Threshold)
plot_LAM_pro <- function(df, subject_id, lam_val, target_RR = 0.05) {
  coords_mat <- scale(as.matrix(df))
  n_points   <- nrow(coords_mat)
  dist_mat   <- rdist(coords_mat)
  
  # Dynamically calculate the 5% threshold specifically for this subject
  dist_values <- dist_mat[lower.tri(dist_mat)]
  subj_radius <- quantile(dist_values, probs = target_RR, na.rm = TRUE)
  
  rec_mat     <- dist_mat <= subj_radius
  rec_indices <- which(rec_mat, arr.ind = TRUE)
  rec_df      <- as.data.frame(rec_indices)
  colnames(rec_df) <- c("Time_A", "Time_B")
  
  # Calculate actual RR (excluding the Line of Identity to match standard RQA)
  actual_RR <- (nrow(rec_df) - n_points) / (n_points^2 - n_points)
  
  rr_str    <- format(round(actual_RR, 4), nsmall = 4)
  lam_str   <- format(round(lam_val, 4), nsmall = 4)
  eps_str   <- format(round(subj_radius, 4), nsmall = 4) 
  
  ggplot(rec_df, aes(x = Time_A, y = Time_B)) +
    geom_point(size = 0.2, shape = 15) + 
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid = element_blank()
    ) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    labs(title = paste0("LAM Plot (Subject ", subject_id, ")"),
         subtitle = paste0("RR = ", rr_str, " | LAM = ", lam_str, " | eps (\u03b5) = ", eps_str),
         x = "Time (TR)", y = "Time (TR)")
}

# ==========================================
# 3. Main Batch Loop
# ==========================================
cat("\n=================================================\n")
cat("🎨 Generating Trajectory, GE, and LAM Plots...\n")
cat("=================================================\n")

for (i in seq_along(file_list)) {
  subject_name <- path_file(file_list[i]) |> sub("_optimal_MDSpoints\\.csv|_MDSpoints\\.csv", "", x = _)
  
  subj_indices <- indices_df %>% filter(Subject == subject_name)
  if (nrow(subj_indices) == 0) {
    cat(sprintf("  --> Skipping %s: No summary data found.\n", subject_name))
    next
  }
  
  df_points <- read_csv(file_list[i], show_col_types = FALSE)
  
  # A. Trajectory Plot
  p_traj <- plot_2d_trajectory(df_points, subject_name)
  ggsave(path(dir_traj, paste0(subject_name, "_2D_Trajectory.png")), p_traj, width = 7, height = 5, dpi = 300, bg = "white")
  
  # B. GE Heatmap
  p_ge <- plot_ge_heatmap(df_points, subject_name, subj_indices$GE_Raw[1], subj_indices$GE_Norm[1])
  ggsave(path(dir_ge, paste0(subject_name, "_GE_Heatmap.png")), p_ge, width = 6.5, height = 6, dpi = 300)
  
  # C. LAM Plot (Target RR set to 5%)
  p_lam <- plot_LAM_pro(df_points, subject_name, subj_indices$LAM_Norm[1], target_RR = 0.05)
  ggsave(path(dir_lam, paste0(subject_name, "_LAMPlot_Pro.png")), p_lam, width = 6, height = 6.2, dpi = 300)
  
  cat(sprintf("  ✅ Completed visuals for: %s\n", subject_name))
}

cat("\n=================================================\n")
cat("🎉 Pipeline Complete! All visuals saved to:", path(output_dir, "plots"), "\n")