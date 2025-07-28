#!/bin/bash

# Script to clean ANSI escape codes from LHE files
# Usage: ./clean_lhe_files.sh [directory]

TARGET_DIR=${1:-"events/"}

echo "Scanning for LHE files with ANSI escape codes in: $TARGET_DIR"

find "$TARGET_DIR" -name "unweighted_events.lhe.gz" | while read file; do
    echo "Checking: $file"
    
    # Check if file contains ANSI escape codes
    if zcat "$file" | head -100 | grep -q $'\x1b\['; then
        echo "  -> Found ANSI codes, cleaning..."
        
        # Create backup
        cp "$file" "${file}.backup"
        
        # Clean and replace
        zcat "$file" | sed 's/\x1b\[[0-9;]*m//g' | gzip > "${file}.tmp"
        mv "${file}.tmp" "$file"
        
        echo "  -> Cleaned successfully (backup saved as ${file}.backup)"
    else
        echo "  -> Already clean"
    fi
done

echo "Done!" 