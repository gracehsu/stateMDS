# load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(here)
library(fs)

# 1. --- Input Setup ---
input_dir      <- here("output", "MDSpoint")
summary_file   <- here("output", "Brain_Indices_Summary.csv")
output_viz_dir <- here("output", "plots", "ge_heatmaps")
dir_create(output_viz_dir)

# Load the summary results to get GE values
if (!file.exists(summary_file)) {
  stop("Brain_Indices_Summary.csv not found. Please run Step 4 first.")
}
indices_df <- read_csv(summary_file, show_col_types = FALSE)

file_list <- dir_ls(input_dir, regexp = "_MDSpoints.csv$")

# 2. --- Define Final Plotting Function ---
plot_ge_heatmap <- function(df, subject_id, ge_raw, ge_norm, bins = 25) {
  
  colnames(df)[1:2] <- c("MDS1", "MDS2")
  background_color  <- "#FFFFE0"
  warm_colors       <- c("#FFD700", "#FF8C00", "#FF4500", "#B22222", "#8B0000")
  
  # Format both values to 4 decimal places
  h_raw_str  <- format(round(ge_raw, 4), nsmall = 4)
  h_norm_str <- format(round(ge_norm, 4), nsmall = 4)
  
  ggplot(df, aes(x = MDS1, y = MDS2)) +
    geom_bin2d(bins = bins) + 
    scale_fill_gradientn(
      colors = warm_colors,
      name = "TR Count",
      limits = c(1, NA) 
    ) +
    theme_minimal() +
    theme(aspect.ratio = 1) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = background_color, color = NA),
      panel.grid.major = element_line(color = "white", size = 0.2),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0("Subject ", subject_id, " (H_raw = ", h_raw_str, ", H_norm = ", h_norm_str, ")"),
      subtitle = "Grid Occupancy (NMDS1 vs NMDS2)",
      x = "NMDS1",
      y = "NMDS2"
    )
}

# 3. --- Main Loop ---
message("Generating Final GE Heatmaps with Raw and Norm values...")

for (i in seq_along(file_list)) {
  subject_name <- path_file(file_list[i]) |> sub("_MDSpoints.csv", "", x = _)
  
  # Filter for the correct row in indices_df
  subj_indices <- indices_df %>% filter(Subject == subject_name)
  
  if (nrow(subj_indices) == 0) {
    warning(paste("No indices found for", subject_name, "- Skipping."))
    next
  }
  
  # Pull both Raw and Norm values
  ge_raw_val  <- subj_indices$GE_Raw[1]
  ge_norm_val <- subj_indices$GE_Norm[1]
  
  df <- read_csv(file_list[i], show_col_types = FALSE)
  
  # Create plot
  p <- plot_ge_heatmap(df, subject_name, ge_raw_val, ge_norm_val)
  
  # Save
  plot_filename <- path(output_viz_dir, paste0(subject_name, "_GE_Heatmap.png"))
  ggsave(plot_filename, p, width = 6.5, height = 6, dpi = 300)
  
  message(paste("Saved:", subject_name))
}

message("Success! Figures saved in: output/plots/ge_heatmaps/")