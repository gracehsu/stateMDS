# load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(here)
library(fs) # Added for easier file and directory management

# 1. --- Setup Paths & Directories ---
input_dir      <- here("output", "MDSpoint")
output_viz_dir <- here("output", "plots", "2d_trajectories")

# Automatically create the output directory if it doesn't exist
dir_create(output_viz_dir)

# Get a list of all MDS points files in the folder
file_list <- dir_ls(input_dir, regexp = "_MDSpoints\\.csv$")

if (length(file_list) == 0) {
  stop("No MDS points files found. Did you run the analysis script first?")
}

message("🚀 Generating and saving 2D Trajectory Plots for all subjects...")

# 2. --- Main Loop ---
for (file_path in file_list) {
  
  # Dynamically extract subject ID from the filename (e.g., "subject1" from "subject1_MDSpoints.csv")
  subject_id <- path_file(file_path) |> sub("_MDSpoints\\.csv$", "", x = _)
  
  # Load and Prepare Data
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Skip if the file is empty
  if(nrow(df) == 0) next
  
  colnames(df) <- c("MDS1", "MDS2")
  df$TR <- 1:nrow(df) # Add time column
  
  # ---------------------------------------------------------
  # 3. Create 2D Trajectory Plot (Static)
  # ---------------------------------------------------------
  plot_2d <- ggplot(df, aes(x = MDS1, y = MDS2, color = TR)) +
    geom_path(alpha = 0.5, linewidth = 0.8) +  # The trajectory line (using linewidth for newer ggplot2)
    geom_point(size = 1.5) +                   # The state points
    scale_color_viridis_c(option = "plasma") + 
    theme_minimal() +
    labs(title = paste("2D State Space Trajectory:", subject_id),
         subtitle = "Distance between points reflects state dissimilarity",
         x = "MDS Dimension 1",
         y = "MDS Dimension 2",
         color = "Time (TR)")
  
  # ---------------------------------------------------------
  # 4. Save the Figure
  # ---------------------------------------------------------
  plot_filename <- path(output_viz_dir, paste0(subject_id, "_2D_Trajectory.png"))
  
  # ggsave automatically saves the last generated plot (or the one explicitly passed to it)
  # bg = "white" ensures the background isn't transparent, which helps if viewing in dark mode
  ggsave(plot_filename, plot = plot_2d, width = 7, height = 5, dpi = 300, bg = "white")
  
  message(paste("✅ Saved plot for:", subject_id))
}

message("🎉 Done! All 2D plots have been successfully saved to: ", output_viz_dir)