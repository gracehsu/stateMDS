# stateMDS: Functional State Dynamics via Multidimensional Scaling

[![R-Project](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🧠 Overview
**stateMDS** is an R-based framework designed to quantify and visualize the dynamic trajectories of brain network states using resting-state fMRI (rsfMRI) data.

### The Concept
* **Functional State:** Defined as the spatial pattern of multi-voxel activity (BOLD signal intensities) captured at a single time point ($TR$) within a predefined network or ROI.
* **State Transition:** A temporal change in this multi-voxel pattern. By quantifying frame-by-frame changes, we capture the dynamic trajectory of a network’s activity.
* **Dimensionality Reduction:** Using Multidimensional Scaling (MDS), high-dimensional voxel-wise patterns are projected into a lower-dimensional space (e.g., 2D, 3D, or 4D). 
* **Analysis:** Distances between points in this low-dimensional space reflect the dissimilarity between functional states, allowing for the assessment of **network stability**, **state transitions**, and their relation to individual differences such as **psychological/cognitive resilience**.

This method is atlas-agnostic and can be applied to any network of interest (e.g., Default Mode Network, Salience Network) defined by standard atlases like **AAL3**, **DiFuMo**, or custom ROI sets.

---

## 📁 Repository Structure

```text
.
├── data/               # Input data
│   └── voxels/         # Voxel-wise fMRI time series (one CSV per subject)
│       ├── subject1_voxels.csv
│       └── subject2_voxels.csv
│
├── output/             # Generated analysis results
│   ├── arrowdis/       # Stepwise Euclidean "arrow" distances per subject
│   └── MDSpoint/       # Low-dimensional MDS coordinates per subject
│
├── R/                  # Core analysis scripts
│   └── run_mds_analysis.R
│
├── brainMDS.Rproj      # RStudio project file (defines project root)
└── README.md           # Project documentation

```


## 📦 Requirements

Install the following R packages:

```r
install.packages(c("vegan", "ggplot2", "readr", "here", "fs"))


▶️ How to Run the Analysis
Open the DMN_MDS_Analysis.Rproj in RStudio.

Make sure the input voxel data (*_voxels.csv) are in the data/voxels/ folder.

Run the script:
source("run_MDS_analysis.R")

This script will:

Load and preprocess all voxel data

Run MDS on each subject

Compute arrow distances across time

Save outputs into output/arrowdis/ and output/MDSpoint/


📤 Outputs
output/MDSpoint/
For each subject: a CSV with 2D/3D/4D MDS coordinates over time.

output/arrowdis/
For each subject: step-by-step arrow distances across time.

MDS_Velocity_Summary.csv: summary for all subjects (total distance, mean distance, stress, convergence).

📌 Notes
Paths are handled via the here package, so the project must be run from within the project root or with the .Rproj open in RStudio.

