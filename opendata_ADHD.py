import os
import shutil
import numpy as np
import pandas as pd
import subprocess
import nibabel as nib
from datetime import datetime
from nilearn import datasets, image
from nilearn.maskers import NiftiMasker

# ==========================================
# CONFIGURATION BLOCK
# ==========================================
ATLAS = "schaefer" 
TARGET_ROI = "Default" 
INPUT_DIR = "data/voxels"
OUTPUT_DIR = "output_adhd"
TARGET_TR = 2.0  

# ==========================================
# STEP 0: Intelligent Folder Archiving (Q2)
# ==========================================
if os.path.exists(INPUT_DIR) and len(os.listdir(INPUT_DIR)) > 0:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_dir = f"data/voxels_archive_{timestamp}"
    os.rename(INPUT_DIR, archive_dir)
    print(f"📁 Archived old voxel files to: {archive_dir}")

os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# STEP 1: Fetch and Clean Phenotypic Data
# ==========================================
print("\nFetching ADHD-200 Dataset pool (40 subjects)...")
dataset = datasets.fetch_adhd(n_subjects=40)

if isinstance(dataset.phenotypic, list) and len(dataset.phenotypic) > 0 and isinstance(dataset.phenotypic[0], str):
    pheno_df = pd.read_csv(dataset.phenotypic[0])
elif isinstance(dataset.phenotypic, str):
    pheno_df = pd.read_csv(dataset.phenotypic)
else:
    pheno_df = pd.DataFrame(dataset.phenotypic)

pheno_df.columns = [str(c).lower() for c in pheno_df.columns]

for col in pheno_df.columns:
    if pheno_df[col].dtype == object:
        pheno_df[col] = pheno_df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else str(x))

if 'subject' in pheno_df.columns:
    pheno_df['subject'] = pheno_df['subject'].astype(str).str.zfill(7)

# ==========================================
# STEP 2: Filter by TR and Find Global Minimum Volumes (Q3)
# ==========================================
valid_controls = []
valid_patients = []
global_min_volumes = float('inf') # Start with infinity, reduce as we find smaller scans
all_fetched_indices = list(range(len(dataset.func)))

print(f"Scanning headers for TR = {TARGET_TR}s and standardizing scan length...")

for idx in all_fetched_indices:
    func_file = dataset.func[idx]
    img_header = nib.load(func_file).header
    
    # get_zooms()[3] is the TR (in seconds), get_data_shape()[3] is the number of volumes
    actual_tr = img_header.get_zooms()[3]
    num_volumes = img_header.get_data_shape()[3]
    
    if round(actual_tr, 2) == TARGET_TR:
        sub_id = os.path.basename(func_file).split('_')[0]
        
        if 'subject' in pheno_df.columns:
            subject_row = pheno_df[pheno_df['subject'] == sub_id]
            if subject_row.empty: continue 
            diagnosis_code = int(subject_row.iloc[0]['adhd'])
        else:
            diagnosis_code = int(pheno_df.iloc[idx].get('adhd', 0))
        
        # Add to cohort if we need them
        if diagnosis_code == 0 and len(valid_controls) < 6:
            valid_controls.append((idx, sub_id, "CTRL", func_file))
            if num_volumes < global_min_volumes: global_min_volumes = num_volumes
                
        elif diagnosis_code == 1 and len(valid_patients) < 6:
            valid_patients.append((idx, sub_id, "ADHD", func_file))
            if num_volumes < global_min_volumes: global_min_volumes = num_volumes
            
    if len(valid_controls) == 6 and len(valid_patients) == 6:
        break

cohort_list = valid_controls + valid_patients
cohort_indices = [item[0] for item in cohort_list]

print(f"✅ Selected {len(valid_controls)} Controls and {len(valid_patients)} ADHD subjects.")
print(f"⏱️  Standardizing all subjects to exactly {global_min_volumes} TRs.")

# ==========================================
# STEP 3: Automatic Garbage Collection (Q1)
# ==========================================
print("\n🗑️  Cleaning up hard drive (deleting unused subjects from Nilearn cache)...")
deleted_count = 0
for idx in all_fetched_indices:
    if idx not in cohort_indices:
        try:
            os.remove(dataset.func[idx])
            os.remove(dataset.confounds[idx])
            deleted_count += 1
        except OSError:
            pass # File might already be deleted or locked
print(f"Cleared {deleted_count} unused subjects from your storage.")

# ==========================================
# STEP 4: Fetch Atlas and Build Mask
# ==========================================
print(f"\nLoading the {ATLAS.upper()} atlas...")

if ATLAS == "schaefer":
    atlas = datasets.fetch_atlas_schaefer_2018(n_rois=100, yeo_networks=7, resolution_mm=2)
    atlas_img, labels = atlas.maps, atlas.labels
    target_indices = [
        i + 1 for i, label in enumerate(labels) 
        if TARGET_ROI in (label.decode('utf-8') if isinstance(label, bytes) else str(label))
    ]
elif ATLAS == "aal":
    atlas = datasets.fetch_atlas_aal()
    atlas_img, labels = atlas.maps, atlas.labels
    target_indices = [
        int(atlas.indices[i]) for i, label in enumerate(labels) 
        if TARGET_ROI in (label.decode('utf-8') if isinstance(label, bytes) else str(label))
    ]

roi_mask = image.math_img(f"np.isin(img, {target_indices})", img=atlas_img)

# ==========================================
# STEP 5: Process and Truncate the Matrices
# ==========================================
for file_idx, sub_id, diagnosis, func_file in cohort_list:
    file_prefix = f"sub-{sub_id}_{diagnosis}"
    print(f"\nProcessing {file_prefix}...")
    
    confounds_file = dataset.confounds[file_idx] 
    
    masker = NiftiMasker(
        mask_img=roi_mask,
        standardize='zscore_sample',
        detrend=True,
        high_pass=0.01,
        low_pass=0.08,
        t_r=TARGET_TR,
        memory='nilearn_cache',
        verbose=0 
    )

    time_series = masker.fit_transform(func_file, confounds=confounds_file)
    
# THE FIX: Slice the matrix to enforce identical lengths across the cohort
    time_series_truncated = time_series[:global_min_volumes, :]
    
    print(f" -> Extracted Shape: {time_series.shape[0]} TRs × {time_series.shape[1]} Voxels | Truncated to: {time_series_truncated.shape[0]} TRs × {time_series_truncated.shape[1]} Voxels")

    output_csv = os.path.join(INPUT_DIR, f"{file_prefix}_timeseries.csv")
    np.savetxt(output_csv, time_series_truncated, delimiter=",")

# ==========================================
# STEP 6: Trigger stateMDS Bash Script
# ==========================================
print("\nFiring stateMDS pipeline on the standardized cohort...")
try:
    # Pass the global minimum as the -t parameter
    bash_command = f"./run_stateMDS.sh -d {INPUT_DIR} -o {OUTPUT_DIR} -t {global_min_volumes}"
    result = subprocess.run(bash_command, shell=True, check=True, text=True, capture_output=True)
    print(result.stdout)
    
except subprocess.CalledProcessError as e:
    print("\n❌ stateMDS script hit an error:")
    print(e.stderr)
    if e.stdout:
        print(e.stdout)