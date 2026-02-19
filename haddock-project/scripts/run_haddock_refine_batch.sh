#!/bin/bash
#
# Stage 6: Run HADDOCK 3 Refinement Batch
#
# Executes HADDOCK 3 refinement for all prepared runs.
# Requires haddock3 conda environment.
#
# Usage: ./run_haddock_refine_batch.sh
#

set -e

# Configuration — resolve to absolute paths before cd-ing into run dirs
RUN_LIST="$(cd "$(dirname "haddock-project/outputs/haddock_refine/run_dirs.txt")" && pwd)/$(basename "haddock-project/outputs/haddock_refine/run_dirs.txt")"
LOG_FILE="$(mkdir -p haddock-project/logs && cd haddock-project/logs && pwd)/haddock_refine.log"

# Activate HADDOCK 3 environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate haddock3

# Create log directory
mkdir -p "$(dirname "$LOG_FILE")"

echo "Starting HADDOCK 3 refinement batch..." | tee -a "$LOG_FILE"
echo "Run list: $RUN_LIST" | tee -a "$LOG_FILE"
echo "---" | tee -a "$LOG_FILE"

# Process each run directory
while IFS= read -r run_dir; do
    [[ -z "$run_dir" ]] && continue

    echo "Processing: $run_dir" | tee -a "$LOG_FILE"

    # Check if already completed (HADDOCK 3 creates numbered module dirs)
    if [[ -d "$run_dir/run" ]] && ls "$run_dir/run/"*_caprieval 1>/dev/null 2>&1; then
        echo "  Already completed, skipping." | tee -a "$LOG_FILE"
        continue
    fi

    # Run HADDOCK 3
    cd "$run_dir"

    echo "  Running refinement..." | tee -a "$LOG_FILE"
    haddock3 refine.cfg >> "$LOG_FILE" 2>&1 || true

    cd - > /dev/null

    # Check completion
    if ls "$run_dir/run/"*_caprieval 1>/dev/null 2>&1; then
        echo "  Completed successfully." | tee -a "$LOG_FILE"
    else
        echo "  Warning: run may not have completed." | tee -a "$LOG_FILE"
    fi

done < "$RUN_LIST"

echo "---" | tee -a "$LOG_FILE"
echo "Batch refinement complete." | tee -a "$LOG_FILE"
