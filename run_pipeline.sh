#!/usr/bin/env bash

set -euo pipefail

# -------------------------------
# Basic paths
# -------------------------------

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
CONFIG_FILE="${PROJECT_DIR}/config.yaml"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "${LOG_DIR}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/run_${TIMESTAMP}.log"

# -------------------------------
# Dry-run mode
# -------------------------------

DRY_RUN=false

if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo ">>> Running in DRY-RUN mode"
fi

# -------------------------------
# Simple conda activation
# -------------------------------

echo ">>> Activating conda environment"

# Load conda into bash
source ~/anaconda3/etc/profile.d/conda.sh

conda activate var-explorer 2>/dev/null || {
    echo "[ERROR] Could not activate conda environment: var-explorer"
    exit 1
}

# -------------------------------
# Function to run commands
# -------------------------------

run_cmd () {
    echo ">>> $*" | tee -a "${LOG_FILE}"

    if [[ "${DRY_RUN}" == false ]]; then
        "$@" 2>&1 | tee -a "${LOG_FILE}"
    fi
}

# -------------------------------
# Run pipeline
# -------------------------------

echo ">>> Starting Var-Explorer pipeline" | tee -a "${LOG_FILE}"

run_cmd python3 "${SCRIPTS_DIR}/summary.py" --config "${CONFIG_FILE}"
run_cmd python3 "${SCRIPTS_DIR}/report.py"  --config "${CONFIG_FILE}"

echo ">>> Pipeline finished successfully" | tee -a "${LOG_FILE}"
echo ">>> Log saved at: ${LOG_FILE}"

