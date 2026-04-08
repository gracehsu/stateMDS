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
.
├── data/               
│   ├── raw/            # Sample 4D EPI volume (.nii) & Sample mask (.nii)
│   └── voxels/         # Extracted voxel-wise fMRI time series (CSV per subject)
│
├── output/             # Generated analysis results
│   ├── arrowdis/       # Stepwise Euclidean "arrow" distances
│   └── MDSpoint/       # Low-dimensional MDS coordinates
│
├── R/                  # Core R analysis scripts
│   └── run_mds_analysis.R         # run MDS analysis 
│   └── visualize_trajectories.R   # Visualize the trajectories 
│  
├── matlab/             # Preprocessing scripts
│   └── catCarryingVoxel.m  # Extracts voxel activity from 4D NIfTI to CSV
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
install.packages(c("vegan", "ggplot2", "readr", "here", "fs", "dplyr"))


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

Step 3: Visualize the Trajectories (R)
1. Open the Visualization Script: In RStudio, open `R/viz/plot_trajectories.R`.
2. Select Subject: Change the `subject_id` variable at the top of the script to the ID of the subject you wish to visualize (e.g., `"subject1"`).
3. Generate Plots:
   2D View: Run the first half of the script to generate a static `ggplot2` path. This is ideal for checking the overall "territory" the brain network covered.
   3D Interactive View: Run the `plotly` section to generate an interactive 3D spatiotemporal trajectory.
4. Interact: Use your mouse to click and drag the 3D plot. By placing **Time (TR)** on the Z-axis, you can visualize the "flow" of brain states; vertical stacks indicate state stability, while wide horizontal leaps indicate rapid state transitions.
5. Export: Save the 2D plot using the "Export" button in the Plots pane.
   Save the 3D plot as an interactive `.html` file using `htmlwidgets::saveWidget()` or the "Save as Web Page" option in the Viewer pane.


📤 Outputs
data/voxels/
For each subject: a CSV file containing voxel x time matrices extracted via MATLAB.

output/MDSpoint/
For each subject: a CSV with 2D/3D/4D MDS coordinates over time.

output/arrowdis/
For each subject: a CSV containing step-by-step displacement distances.

MDS_Velocity_Summary.csv: Group-level summary including total distance, mean distance, MDS stress, and convergence metrics.

📌 Implementation Notes
Extraction Modes: catCarryingVoxel provides three modes. Use Mode 1 if your mask has a smaller voxel size than your data, or Mode 2 if you prefer to retain every voxel in the data space.

NIfTI Data: A sample 4D EPI volume.nii is included in data/raw/ for testing your extraction parameters.

Pathing: Always use the .Rproj file to maintain relative pathing via the here package.


