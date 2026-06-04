# load libraries
library(ggplot2)
library(plotly)
library(dplyr)
library(readr)
library(here)

# 1. Setup paths
# Change "subject1" to match one of your actual subject IDs in the folder
subject_id <- "subject1" 
input_file <- here("output", "MDSpoint", paste0(subject_id, "_MDSpoints.csv"))

if (!file.exists(input_file)) {
  stop("MDS points file not found. Did you run the analysis script first?")
}

# 2. Load and Prepare Data
df <- read_csv(input_file, show_col_types = FALSE)
colnames(df) <- c("MDS1", "MDS2")
df$TR <- 1:nrow(df) # Add time column

# ---------------------------------------------------------
# 3. 2D Trajectory Plot (Static)
# ---------------------------------------------------------
plot_2d <- ggplot(df, aes(x = MDS1, y = MDS2, color = TR)) +
  geom_path(alpha = 0.5, size = 0.8) +  # The trajectory line
  geom_point(size = 1.5) +            # The state points
  scale_color_viridis_c(option = "plasma") + 
  theme_minimal() +
  labs(title = paste("2D State Space Trajectory:", subject_id),
       subtitle = "Distance between points reflects state dissimilarity",
       x = "MDS Dimension 1",
       y = "MDS Dimension 2",
       color = "Time (TR)")

print(plot_2d)

# ---------------------------------------------------------
# 4. 3D Rotating Trajectory (MDS1, MDS2 + Time)
# ---------------------------------------------------------
# We use TR as the Z-axis to show the "flow" of time
plot_3d <- plot_ly(df, 
                   x = ~MDS1, 
                   y = ~MDS2, 
                   z = ~TR, 
                   type = 'scatter3d', 
                   mode = 'lines+markers',
                   line = list(width = 4, color = ~TR, colorscale = 'Viridis'),
                   marker = list(size = 3, color = ~TR, colorscale = 'Viridis', opacity = 0.8)) %>%
  layout(title = paste("3D Spatiotemporal Trajectory:", subject_id),
         scene = list(
           xaxis = list(title = 'MDS 1'),
           yaxis = list(title = 'MDS 2'),
           zaxis = list(title = 'Time (TR)')
         ))

# This will open in your RStudio Viewer or Default Browser
plot_3d
