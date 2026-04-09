# stateMDS: Functional State Dynamics via Multidimensional Scaling

[![R-Project](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Matlab](https://img.shields.io/badge/Language-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🧠 Overview
**stateMDS** is a framework designed to quantify and visualize the dynamic trajectories of brain network states using resting-state fMRI (rsfMRI) data.

### The Concept
* **Functional State:** Defined as the spatial pattern of multi-voxel activity (BOLD signal intensities) captured at a single time point ($TR$) within a predefined network or ROI.
* **State Transition:** A temporal change in this multi-voxel pattern. By quantifying frame-by-frame changes, we capture the dynamic trajectory of a network’s activity.
* **Dimensionality Reduction:** Using Multidimensional Scaling (MDS), high-dimensional voxel-wise patterns are projected into a lower-dimensional space (e.g., 2D, 3D, or 4D). 

This method is atlas-agnostic and can be applied to any network of interest (e.g., Default Mode Network, Salience Network) defined by standard atlases like **AAL3**, **DiFuMo**, or custom ROI sets.

---

## 📁 Repository Structure

```text

├── data/               
│   ├── raw/            # Sample 4D EPI volume (.nii) & Sample mask (.nii)
│   └── voxels/         # Extracted voxel-wise fMRI time series (CSV per subject)
│
├── output/             # Generated analysis results
│   ├── arrowdis/       # Stepwise Euclidean "arrow" distances
│   ├── MDSpoint/       # Low-dimensional MDS coordinates
│   └── plots/          # Visualization outputs
│       ├── ge_heatmaps/   # Grid Entropy heatmaps
│       └── LAM_plots/     # Recurrence plots for Laminarity
│
├── matlab/             # Preprocessing scripts
│   └── catCarryingVoxel.m           # Step 1: Voxel extraction
│  
├── R/                  # Core R analysis and visualization scripts
│   ├── run_mds_analysis.R           # Step 2: Main MDS pipeline
│   ├── calculate_brain_indices.R    # Step 3: Compute CHA, GE, and LAM
│   ├── visualize_trajectories.R     # Step 4a: 2D/3D Trajectory plotting
│   ├── plot_state_space_heatmap.R   # Step 4b: GE Grid occupancy heatmaps
│   └── plot_LAM_plots.R             # Step 4c: Recurrence plots
│
├── stateMDS.Rproj      # RStudio project file
└── README.md           # Project documentation

```

## 📦 Requirements

MATLAB (Extraction Pipeline)
To run the extraction using catCarryingVoxel.m, you will need:

MATLAB (R2016b or later).

SPM12 (Tested on revision 6906). The spm_vol, spm_read_vols, and spm_get_data functions must be in your MATLAB path.

Install the following R packages:

```r
install.packages(c("vegan", "ggplot2", "readr", "here", "fs", "dplyr", "plotly", "geometry", "entropy", "fields"))

▶️ How to Run the Analysis

Step 1: Voxel Extraction (MATLAB)
Use catCarryingVoxel.m to extract values from your 4D EPI volumes using a mask (e.g., AAL3 or DiFuMo).

1. Open MATLAB and add SPM12 to your path.

2. Navigate to the matlab/ folder.

3. Call the function using your mask and data:
% Example: mode 2 (base on data space), threshold 0.5
[meanval, voxelval, voxmni, voxcor] = catCarryingVoxel('my_mask.nii', 'data/raw/4D_EPI_volume.nii', 2, 0.5);

4. Save the resulting voxelval as a CSV in data/voxels/ (e.g., subject1_voxels.csv).

Step 2: MDS Analysis (R)
1. Open Project: Launch brainMDS.Rproj in RStudio.

2. Execute Script: Run the following command in the R console:

source("R/run_mds_analysis.R")

The Pipeline Will:
1. Load and preprocess voxel-level time series.

2. Perform Multidimensional Scaling (MDS) for each subject.

3. Calculate frame-to-frame "arrow distances" (velocity) across the scan duration.

4. Export coordinates and summary statistics to the output/ folder.

Step 3: Calculate Brain Indices (R)
Compute high-level geometric and dynamic metrics:

Step 4: Visualization (R)
Choose a script based on your visualization needs:

Trajectories: Run R/visualize_trajectories.R for interactive 3D "corkscrew" views (MDS1, MDS2, Time) and 2D paths.
See a [sample 3D trajectory visualization here](assets/subject1_3D_screenshot.png).

Grid Entropy: Run R/plot_state_space_heatmap.R to generate warm-colored grid occupancy heatmaps in output/plots/ge_heatmaps/.

Laminarity: Run R/plot_LAM_plots.R to generate Laminarity plots in output/plots/LAM_plots/.


📤 Outputs
data/voxels/
For each subject: a CSV file containing voxel x time matrices extracted via MATLAB.

output/MDSpoint/
For each subject: a CSV with 2D/3D/4D MDS coordinates over time.

output/arrowdis/
For each subject: a CSV containing step-by-step displacement distances.

MDS_Velocity_Summary.csv: Group-level summary including total distance, mean distance, MDS stress, and convergence metrics.

output/Brain_Indices_Summary.csv: CHA_Raw/Norm, GE_Raw/Norm, LAM_Norm.

output/plots/
PNG files for GE heatmaps, and LAM recurrence plots.

📌 Implementation Notes
Extraction Modes: catCarryingVoxel provides three modes. Use Mode 1 if your mask has a smaller voxel size than your data, or Mode 2 if you prefer to retain every voxel in the data space.

NIfTI Data: A sample 4D EPI volume.nii is included in data/raw/ for testing your extraction parameters.

Pathing: Always use the .Rproj file to maintain relative pathing via the here package.


