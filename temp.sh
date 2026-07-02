#!/bin/bash

# Configuration
SOURCE_DIR="harita-sims"
N=8

# Get the absolute path of your current working directory
TARGET_DIR=$(pwd)

# Move into the source directory safely
cd "$SOURCE_DIR" || { echo "Source directory '$SOURCE_DIR' not found!"; exit 1; }

echo "Starting copy process for N=$N..."
echo "Target directory: $TARGET_DIR"
echo "-----------------------------------"

# Loop i from 1 to N-1
for (( i=1; i<N; i++ )); do
    # Loop j from i+1 to N
    for (( j=i+1; j<=N; j++ )); do
        
        for ext in inp out; do
            # Since we are already inside the SOURCE_DIR, match the pattern locally
            PATTERN="*frag${i}_frag${j}.${ext}"
            
            # Use a nullglob so the loop doesn't literally search for the string '*frag...' 
            # if no files match
            shopt -s nullglob
            for file in $PATTERN; do
                echo "-> Copying to current dir: $file"
                cp "$file" "$TARGET_DIR/"
            done
            shopt -u nullglob # Turn off nullglob to keep environment clean
            
        done

    done
done

echo "-----------------------------------"
echo "Done! All valid pairs (up to N=$N) have been copied here."