# load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(here)
library(fs)
library(fields)

# 1. --- Input Setup ---
input_dir      <- here("output", "MDSpoint")
summary_file   <- here("output", "Brain_Indices_Summary.csv")
output_viz_dir <- here("output", "plots", "recurrence_plots")
dir_create(output_viz_dir)

if (!file.exists(summary_file)) {
  stop("Brain_Indices_Summary.csv not found. Please run calculate_brain_indices.R first.")
}
indices_df <- read_csv(summary_file, show_col_types = FALSE)
file_list  <- dir_ls(input_dir, regexp = "_MDSpoints.csv$")


# 2. --- PRE-COMPUTATION: Find a single Global Epsilon ---
message("🔍 Pre-scanning dataset to calculate global epsilon...")
all_eps_values <- c()

for (f in file_list) {
  df_temp <- read_csv(f, show_col_types = FALSE)
  # Standardize to make distances comparable across subjects
  coords_mat <- scale(as.matrix(df_temp)) 
  dist_mat <- rdist(coords_mat)
  dist_values <- dist_mat[lower.tri(dist_mat)]
  
  # Find what the 5% threshold would be for this specific subject
  subj_eps <- quantile(dist_values, probs = 0.05, na.rm = TRUE)
  all_eps_values <- c(all_eps_values, subj_eps)
}

# Calculate the average threshold across everyone to use as our global standard
global_epsilon <- mean(all_eps_values)
message(paste("✅ Global epsilon calculated as:", round(global_epsilon, 4)))


# 3. --- Define Updated Plotting Function (Fixed Epsilon) ---
# Notice we removed target_RR and replaced it with a fixed epsilon parameter
plot_recurrence_pro <- function(df, subject_id, lam_val, fixed_eps) {
  
  # Standardize
  coords_mat <- scale(as.matrix(df))
  n_points   <- nrow(coords_mat)
  
  # Calculate distance matrix
  dist_mat    <- rdist(coords_mat)
  
  # Create Recurrence Data using the GLOBAL threshold!
  rec_mat     <- dist_mat <= fixed_eps
  rec_indices <- which(rec_mat, arr.ind = TRUE)
  rec_df      <- as.data.frame(rec_indices)
  colnames(rec_df) <- c("Time_A", "Time_B")
  
  # Calculate Actual RR for the title (This will now vary beautifully between subjects!)
  actual_RR <- nrow(rec_df) / (n_points^2)
  
  # Formatting strings for the title
  rr_str  <- format(round(actual_RR, 4), nsmall = 4)
  lam_str <- format(round(lam_val, 4), nsmall = 4)
  eps_str <- format(round(fixed_eps, 4), nsmall = 4) # Updated to show the fixed global eps
  
  ggplot(rec_df, aes(x = Time_A, y = Time_B)) +
    geom_point(size = 0.2, shape = 15) + 
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid = element_blank()
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(
      title = paste0("Recurrence Plot (Subject ", subject_id, ")"),
      subtitle = paste0("RR = ", rr_str, " | LAM = ", lam_str, " | eps (\u03b5) = ", eps_str),
      x = "Time (TR)",
      y = "Time (TR)"
    )
}


# 4. --- Main Loop ---
message("🎨 Generating Professional Recurrence Plots...")

for (i in seq_along(file_list)) {
  subject_name <- path_file(file_list[i]) |> sub("_MDSpoints.csv", "", x = _)
  
  subj_indices <- indices_df %>% filter(Subject == subject_name)
  if (nrow(subj_indices) == 0) next
  
  lam_value <- subj_indices$LAM_Norm[1]
  df_points <- read_csv(file_list[i], show_col_types = FALSE)
  
  # Create plot passing the global_epsilon!
  p <- plot_recurrence_pro(df_points, subject_name, lam_value, global_epsilon)
  
  # Save
  plot_filename <- path(output_viz_dir, paste0(subject_name, "_RecurrencePlot_Pro.png"))
  ggsave(plot_filename, p, width = 6, height = 6.2, dpi = 300)
}

message("🎉 Done! Check output/plots/recurrence_plots/")