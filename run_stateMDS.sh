#!/bin/bash

# Exit immediately if a command crashes
set -e

# ==========================================
# 1. HELP MENU FUNCTION
# ==========================================
show_help() {
    echo "Usage: ./run_stateMDS.sh [Options]"
    echo ""
    echo "Options:"
    echo "  -d  Input directory containing TSV/CSV files (Default: data/voxels)"
    echo "  -o  Main output directory (Default: output)"
    echo "  -t  Maximum TRs to analyze per subject (Default: 180)"
    echo "  -s  Maximum acceptable stress value (Default: 0.15)"
    echo "  -k  Maximum dimension to test (Default: 10)"
    echo "  -h  Show this help message"
    exit 0
}

# ==========================================
# 2. DEFAULT VARIABLES
# ==========================================
INPUT_DIR="data/voxels"
OUTPUT_DIR="output"
MAX_TR=180
STRESS=0.15
MAX_DIM=10

# ==========================================
# 3. PARSE COMMAND LINE ARGUMENTS
# ==========================================
while getopts "d:o:t:s:k:h" opt; do
    case "$opt" in
        d) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) MAX_TR="$OPTARG" ;;
        s) STRESS="$OPTARG" ;;
        k) MAX_DIM="$OPTARG" ;;
        h) show_help ;;
        *) show_help ;;
    esac
done

# ==========================================
# 4. PIPELINE EXECUTION
# ==========================================
echo "========================================================="
echo " 🚀 Starting stateMDS Pipeline"
echo "========================================================="
echo " Input Directory: $INPUT_DIR"
echo " Output Folder:   $OUTPUT_DIR"
echo " Settings:        Max TRs: $MAX_TR | Target Stress: <$STRESS | Max Dim: $MAX_DIM"
echo "========================================================="

# Run Step 1: MDS Analysis
echo -e "\n[Step 1/3] Running Batch Multidimensional Scaling Analysis..."
Rscript R/run_mds_analysis.R \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --max_tr "$MAX_TR" \
    --stress "$STRESS" \
    --max_dim "$MAX_DIM"

# Run Step 2: Brain Indices
echo -e "\n[Step 2/3] Calculating Advanced Brain Indices (CHA, GE, LAM)..."
Rscript R/run_brain_indices.R --output_dir "$OUTPUT_DIR"

# Run Step 3: Visualization
echo -e "\n[Step 3/3] Generating Visualizations and Trajectory Plots..."
Rscript R/visualize_trajectories.R --output_dir "$OUTPUT_DIR"

echo -e "\n========================================================="
echo " ✅ stateMDS Pipeline Completed Successfully!"
echo " All results saved to: $OUTPUT_DIR"
echo "========================================================="