#!/bin/bash

eval "$(/home/hep/us322/miniforge3/bin/conda shell.bash hook)"
conda activate madminer_env

# Parse command line arguments
GENERATION_TYPE=$1
SUPP_ID=$2

# Generate a unique job ID using Condor's PROCESS variable for sequential numbering
JOB_ID=$((PROCESS + 1))
echo "Generated job ID: $JOB_ID"

# Create a unique temporary directory for this job
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"

# Function to generate random seed
generate_random_seed() {
    # Generate a random number between 1 and 999999
    echo $((RANDOM % 999999 + 1))
}

# Function to modify run card with random seed
modify_run_card() {
    local run_card_path=$1
    local temp_run_card_path=$2
    local random_seed=$3
    
    # Copy the original run card to temp directory
    cp "$run_card_path" "$temp_run_card_path"
    
    # Replace XXX with the random seed
    sed -i "s/XXX/${random_seed}/g" "$temp_run_card_path"
    
    echo "Modified run card with seed: $random_seed"
}

# Function to clean up temporary directory
cleanup() {
    echo "Cleaning up temporary directory: $TEMP_DIR"
    rm -rf "$TEMP_DIR"
}

# Set up cleanup on exit
trap cleanup EXIT

# Run the event generation script based on the generation type
case $GENERATION_TYPE in
    "signal")
        echo "Running signal generation (-sm) with random seed"
        
        # Use 14TeV signal run card
        RUN_CARD="./cards/run_cards/run_card_signal_14TeV.dat"
        if [ ! -f "$RUN_CARD" ]; then
            echo "Error: Signal run card not found: $RUN_CARD"
            exit 1
        fi
        
        # Create temporary run card with random seed
        TEMP_RUN_CARD="$TEMP_DIR/run_card_signal_job${JOB_ID}.dat"
        RANDOM_SEED=$(generate_random_seed)
        modify_run_card "$RUN_CARD" "$TEMP_RUN_CARD" "$RANDOM_SEED"
        
        # Run with modified run card
        python 02_generate_events_parallel.py -sm -run_card "$TEMP_RUN_CARD" -job_id "$JOB_ID"
        ;;
        
    "background")
        echo "Running background generation (-b) with random seed"
        
        # Use 14TeV background run card
        RUN_CARD="./cards/run_cards/run_card_background_14TeV.dat"
        if [ ! -f "$RUN_CARD" ]; then
            echo "Error: Background run card not found: $RUN_CARD"
            exit 1
        fi
        
        # Create temporary run card with random seed
        TEMP_RUN_CARD="$TEMP_DIR/run_card_background_job${JOB_ID}.dat"
        RANDOM_SEED=$(generate_random_seed)
        modify_run_card "$RUN_CARD" "$TEMP_RUN_CARD" "$RANDOM_SEED"
        
        # Run with modified run card
        python 02_generate_events_parallel.py -b -run_card "$TEMP_RUN_CARD" -job_id "$JOB_ID"
        ;;
        
    "bsm")
        if [ -z "$SUPP_ID" ]; then
            echo "Error: BSM generation requires a supp_id (1-9)"
            exit 1
        fi
        echo "Running BSM generation with supp_id $SUPP_ID and random seed"
        
        # Use 14TeV signal run card for BSM (same as signal generation)
        RUN_CARD="./cards/run_cards/run_card_signal_14TeV.dat"
        if [ ! -f "$RUN_CARD" ]; then
            echo "Error: Signal run card not found for BSM generation: $RUN_CARD"
            exit 1
        fi
        
        # Create temporary run card with random seed
        TEMP_RUN_CARD="$TEMP_DIR/run_card_signal_job${JOB_ID}.dat"
        RANDOM_SEED=$(generate_random_seed)
        modify_run_card "$RUN_CARD" "$TEMP_RUN_CARD" "$RANDOM_SEED"
        
        # Run with modified run card
        python 02_generate_events_parallel.py -supp -supp_id $SUPP_ID -run_card "$TEMP_RUN_CARD" -job_id "$JOB_ID"
        ;;
        
    *)
        echo "Error: Unknown generation type: $GENERATION_TYPE"
        exit 1
        ;;
esac

echo "Job completed successfully" 