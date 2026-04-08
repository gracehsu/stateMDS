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
│   ├── raw/            # Sample 4D EPI volume (.nii)
│   └── voxels/         # Extracted voxel-wise fMRI time series (CSV per subject)
│
├── output/             # Generated analysis results
│   ├── arrowdis/       # Stepwise Euclidean "arrow" distances
│   └── MDSpoint/       # Low-dimensional MDS coordinates
│
├── R/                  # Core R analysis scripts
│   └── run_mds_analysis.R
│
├── matlab/             # Preprocessing scripts
│   └── catCarryingVoxel.m  # Extracts voxel activity from 4D NIfTI to CSV
│
├── brainMDS.Rproj      # RStudio project file
└── README.md           # Project documentation

```


## 📦 Requirements

MATLAB (Preprocessing)
To run catCarryingVoxel.m, you will need:

MATLAB (R2020b or later recommended)

A NIfTI toolbox (e.g., SPM12 or the NIfTI toolset) added to your MATLAB path.

Install the following R packages:

```r
install.packages(c("vegan", "ggplot2", "readr", "here", "fs", "dplyr"))


▶️ How to Run the Analysis

Step 1: Voxel Extraction (MATLAB)
Run catCarryingVoxel.m to convert your 4D EPI volumes into the required CSV format.

Open MATLAB and navigate to the matlab/ folder.

Ensure your .nii files are in data/raw/.

Run the script to generate subject_voxels.csv files inside data/voxels/.

Step 2: MDS Analysis (R)
Open Project: Launch brainMDS.Rproj in RStudio.

Execute Script: Run the following command in the R console:

source("R/run_mds_analysis.R")

The Pipeline Will:
Load and preprocess voxel-level time series.

Perform Multidimensional Scaling (MDS) for each subject.

Calculate frame-to-frame "arrow distances" (velocity) across the scan duration.

Export coordinates and summary statistics to the output/ folder.



📤 Outputs
data/voxels/
For each subject: a CSV extracted from 4D volumes containing voxel x time matrices.

output/MDSpoint/
For each subject: a CSV with 2D/3D/4D MDS coordinates over time.

output/arrowdis/
For each subject: a CSV containing step-by-step displacement distances.

MDS_Velocity_Summary.csv: Group-level summary including total distance, mean distance, MDS stress, and convergence metrics.

📌 Notes
Paths are handled via the here package, so the project must be run from within the project root or with the .Rproj open in RStudio.

