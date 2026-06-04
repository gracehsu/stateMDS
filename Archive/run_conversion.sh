#!/bin/bash

# Define your directories
RAW_DIR="."
BIDS_DIR="stateMDS_BIDS"
CONFIG="bids_config.json"

# Loop through all ADMS folders in the current directory
for SUBJ_DIR in ${RAW_DIR}/ADMS*; do
  
  # Extract just the subject ID (e.g., "ADMS024")
  SUBJ_ID=$(basename ${SUBJ_DIR})
  
  echo "----------------------------------------"
  echo "🔄 Converting sub-${SUBJ_ID} to BIDS..."
  echo "----------------------------------------"
  
  # Run dcm2bids
  dcm2bids -d ${SUBJ_DIR} \
           -p ${SUBJ_ID} \
           -c ${CONFIG} \
           -o ${BIDS_DIR} \
           --auto_extract_entities
           
done

echo "✅ All subjects successfully converted to BIDS!"
