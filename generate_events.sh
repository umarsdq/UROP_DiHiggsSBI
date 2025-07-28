#!/bin/bash

eval "$(/home/hep/us322/miniforge3/bin/conda shell.bash hook)"
conda activate madminer_env
#export TMPDIR=/vols/cms/us322/tmp

# Run the event generation script
python 02_generate_events.py -sm