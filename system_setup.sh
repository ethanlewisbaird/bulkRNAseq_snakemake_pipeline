#!/bin/bash

# Snakemake Pipeline Setup Script
# This script automates the setup of you system for the RNA-seq analysis pipeline
# It only needs to be run once per system

set -e  # Exit on error

# Configuration variables
CONDA_ENV_NAME="snakemake6"
SNAKEMAKE_VERSION="6.15.5"
SLURM_PROFILE_DIR="$HOME/snakemake/slurm"

# Function to print colored output
print_status() {
    echo -e "\033[1;34m==> $1\033[0m"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if conda is available
if ! command_exists conda; then
    echo "Error: Conda is not installed or not in PATH"
    exit 1
fi

# Create conda environment
print_status "Creating conda environment: $CONDA_ENV_NAME"
if conda env list | grep -q "$CONDA_ENV_NAME"; then
    echo "Conda environment $CONDA_ENV_NAME already exists"
else
    conda create -y -n "$CONDA_ENV_NAME" snakemake="$SNAKEMAKE_VERSION" -c bioconda -c conda-forge
fi

# Create SLURM profile directory
print_status "Setting up SLURM profile"
mkdir -p "$SLURM_PROFILE_DIR"
mkdir -p logs/slurm

# Create SLURM profile config
cat > "$SLURM_PROFILE_DIR/config.yaml" << 'EOF'
restart-times: 2
max-jobs-per-second: 10
max-status-checks-per-second: 30
local-cores: 1
latency-wait: 60
jobs: 5
conda-frontend: conda
use-conda: true
conda-prefix: ~/.snakemake/conda
scheduler: greedy
cluster-status: ~/.config/snakemake/slurm/slurm_status.py

cluster: "sbatch --parsable --cpus-per-task={threads} --mem={resources.mem_mb} --time={resources.time_min} --job-name=smk.{rule}.{wildcards} --output=logs/slurm/{rule}_{wildcards}.%j.out --error=logs/slurm/{rule}_{wildcards}.%j.err"
EOF

# Create slurm_status.py script
cat > "$SLURM_PROFILE_DIR/slurm_status.py" << 'EOF'
#!/usr/bin/env python3
import subprocess
import sys

def main():
    jobid = sys.argv[1]
    try:
        # Use scontrol which is faster than squeue
        output = subprocess.run(
            f"scontrol show job {jobid} | grep JobState",
            shell=True, 
            capture_output=True, 
            text=True,
            timeout=10  # Add timeout to prevent hanging
        )
        
        if output.returncode == 0 and output.stdout:
            state = output.stdout.split("=")[1].split()[0]
            # Map SLURM states to Snakemake states
            state_mapping = {
                "RUNNING": "running",
                "COMPLETED": "success", 
                "FAILED": "failure",
                "CANCELLED": "failure",
                "PENDING": "running",
                "CONFIGURING": "running"
            }
            print(state_mapping.get(state, "failed"))
        else:
            print("failed")
    except subprocess.TimeoutExpired:
        print("running")  # Assume still running if timeout
    except Exception as e:
        print(f"failed")  # Any other error

if __name__ == "__main__":
    main()
EOF

# Make the status script executable
chmod +x "$SLURM_PROFILE_DIR/slurm_status.py"

print_status "Setup complete!"