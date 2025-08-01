#!/bin/bash

# Script to copy events from parallel job directories to the correct directory structure for Delphes processing
# Based on the structure expected by 03a_read_delphes.py
# Handles parallel structure: signal_sm_1, signal_sm_2, etc. -> run_01, run_02, etc.

set -e  # Exit on any error

echo "Setting up directory structure for Delphes processing from parallel jobs..."

# Create base directory
FINAL_EVENTS_DIR="/vols/cms/us322/02b_final_events_14"
SOURCE_DIR="/vols/cms/us322/02_event_generation_14"

mkdir -p "$FINAL_EVENTS_DIR"

echo "Creating directory structure..."

# Function to find and copy signal_sm jobs
copy_signal_sm_jobs() {
    echo "Copying signal_sm events from parallel jobs..."
    mkdir -p "$FINAL_EVENTS_DIR/signal_sm/batch_0"
    
    # Find all signal_sm_* directories and copy them as run_01, run_02, etc.
    job_count=1
    for job_dir in "$SOURCE_DIR/mg_processes"/signal_sm_*; do
        if [ -d "$job_dir" ]; then
            echo "Copying from $job_dir to run_$(printf "%02d" $job_count)..."
            cp -r "$job_dir/Events/run_01" "$FINAL_EVENTS_DIR/signal_sm/batch_0/run_$(printf "%02d" $job_count)"
            cp -r "$job_dir/Events/run_01_decayed_1" "$FINAL_EVENTS_DIR/signal_sm/batch_0/run_$(printf "%02d" $job_count)_decayed_1"
            ((job_count++))
        fi
    done
    
    echo "Copied $((job_count-1)) signal_sm jobs"
}

# Function to find and copy background jobs
copy_background_jobs() {
    echo "Copying background events from parallel jobs..."
    mkdir -p "$FINAL_EVENTS_DIR/background_0/batch_0"
    
    # Find all background_* directories and copy them as run_01, run_02, etc.
    job_count=1
    for job_dir in "$SOURCE_DIR/mg_processes_2"/background_*; do
        if [ -d "$job_dir" ]; then
            echo "Copying from $job_dir to run_$(printf "%02d" $job_count)..."
            cp -r "$job_dir/Events/run_01" "$FINAL_EVENTS_DIR/background_0/batch_0/run_$(printf "%02d" $job_count)"
            ((job_count++))
        fi
    done
    
    echo "Copied $((job_count-1)) background jobs"
}

# Function to find and copy BSM jobs
copy_bsm_jobs() {
    echo "Copying signal_supp events from parallel jobs..."
    mkdir -p "$FINAL_EVENTS_DIR/signal_supp"
    
    # For each morphing basis vector (0-9), find all parallel jobs and copy them
    for mb_vector in {0..9}; do
        # Check if any morphing_basis_vector_${mb_vector}_* directories exist
        if ls "$SOURCE_DIR/mg_processes/signal_supp"/morphing_basis_vector_${mb_vector}_* 1> /dev/null 2>&1; then
            echo "Processing morphing_basis_vector_$mb_vector..."
            mkdir -p "$FINAL_EVENTS_DIR/signal_supp/mb_vector_$mb_vector/batch_0"
            
            job_count=1
            for job_dir in "$SOURCE_DIR/mg_processes/signal_supp"/morphing_basis_vector_${mb_vector}_*; do
                if [ -d "$job_dir" ]; then
                    echo "Copying from $job_dir to run_$(printf "%02d" $job_count)..."
                    cp -r "$job_dir/Events/run_01" "$FINAL_EVENTS_DIR/signal_supp/mb_vector_$mb_vector/batch_0/run_$(printf "%02d" $job_count)"
                    cp -r "$job_dir/Events/run_01_decayed_1" "$FINAL_EVENTS_DIR/signal_supp/mb_vector_$mb_vector/batch_0/run_$(printf "%02d" $job_count)_decayed_1"
                    ((job_count++))
                fi
            done
            
            echo "Copied $((job_count-1)) jobs for morphing_basis_vector_$mb_vector"
        else
            echo "No jobs found for morphing_basis_vector_$mb_vector, skipping..."
        fi
    done
}

# Execute the copy functions
copy_signal_sm_jobs
copy_background_jobs
copy_bsm_jobs

echo "Directory structure created successfully!"
echo "Final structure:"
echo "  $FINAL_EVENTS_DIR/signal_sm/batch_0/run_01, run_02, ..."
echo "  $FINAL_EVENTS_DIR/background_0/batch_0/run_01, run_02, ..."
echo "  $FINAL_EVENTS_DIR/signal_supp/mb_vector_X/batch_0/run_01, run_02, ..." 