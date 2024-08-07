#!/bin/bash

# Ensure conda is initialized
source ~/miniconda3/etc/profile.d/conda.sh

# TD2
conda deactivate
conda activate td2
pip install .

NUM_THREADS=32

BASE_COMMAND="/usr/bin/time -v TD2.LongOrfs -t data/chess.fa"
OUTPUT_DIR="results/chess_benchmark/td2/"

# Remove the output directory if it exists, then create it
if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"
fi
mkdir -p "$OUTPUT_DIR"

# File to store the results
RESULT_FILE="${OUTPUT_DIR}time_stats.txt"
echo "User Time (seconds), System Time (seconds), Elapsed Time (seconds), Maximum Memory (kbytes)" > $RESULT_FILE

# Construct the command
COMMAND="$BASE_COMMAND -O $OUTPUT_DIR -@ $NUM_THREADS > ${OUTPUT_DIR}time.log 2>&1"
eval $COMMAND

# Extract relevant statistics from the log file
user_time=$(grep "User time (seconds):" "${OUTPUT_DIR}time.log" | awk '{print $4}')
system_time=$(grep "System time (seconds):" "${OUTPUT_DIR}time.log" | awk '{print $4}')
elapsed_time=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" "${OUTPUT_DIR}time.log" | awk '{print $8}')
max_memory=$(grep "Maximum resident set size (kbytes):" "${OUTPUT_DIR}time.log" | awk '{print $6}')

# Convert elapsed time to seconds
IFS=: read -r min sec <<< "$elapsed_time"
elapsed_time_sec=$(echo "$min * 60 + $sec" | bc)

# Append individual run stats to the results file
echo "$user_time, $system_time, $elapsed_time_sec, $max_memory" >> $RESULT_FILE

# Output stats to console
echo "Results for TD2:"
echo "User Time (seconds): $user_time"
echo "System Time (seconds): $system_time"
echo "Elapsed Time (seconds): $elapsed_time_sec"
echo "Maximum Memory (kbytes): $max_memory"