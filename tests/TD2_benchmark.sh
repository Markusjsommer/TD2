#!/bin/bash

# conda activate td2
# pip install .

# Base command to be executed
BASE_COMMAND="/usr/bin/time -v TD2.LongOrfs -O results/TD2_benchmark/"

# Number of times to run the command
NUM_RUNS=5

# Array of partition files
PARTITION_FILES=("data/partition_20.fna" "data/partition_40.fna" "data/partition_60.fna" "data/partition_80.fna" "data/partition_100.fna")

# Iterate over each partition file
for PARTITION_FILE in "${PARTITION_FILES[@]}"; do
    # Extract partition name for result files
    PARTITION_NAME=$(basename "$PARTITION_FILE" .fna)
    
    # Directory for the current partition results
    OUTPUT_DIR="results/TD2_benchmark/${PARTITION_NAME}/"
    mkdir -p "$OUTPUT_DIR"
    
    # Initialize total variables for each stat
    total_user_time=0
    total_system_time=0
    total_elapsed_time=0
    total_max_memory=0

    # File to store the results
    RESULT_FILE="${OUTPUT_DIR}time_stats.txt"
    echo "Run, User Time (seconds), System Time (seconds), Elapsed Time (seconds), Maximum Memory (kbytes)" > $RESULT_FILE

    # Execute the command multiple times and collect stats
    for i in $(seq 1 $NUM_RUNS); do
        # Remove the output directory if it exists
        if [ -d "$OUTPUT_DIR" ]; then
            rm -rf "$OUTPUT_DIR"
            mkdir -p "$OUTPUT_DIR"
        fi

        # Construct the command with the current partition file
        COMMAND="$BASE_COMMAND${PARTITION_NAME}/ -t $PARTITION_FILE > ${OUTPUT_DIR}time.log 2>&1"

        eval $COMMAND

        # Extract relevant statistics from the log file
        user_time=$(grep "User time (seconds):" "${OUTPUT_DIR}time.log" | awk '{print $4}')
        system_time=$(grep "System time (seconds):" "${OUTPUT_DIR}time.log" | awk '{print $4}')
        elapsed_time=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" "${OUTPUT_DIR}time.log" | awk '{print $8}')
        max_memory=$(grep "Maximum resident set size (kbytes):" "${OUTPUT_DIR}time.log" | awk '{print $6}')

        # Convert elapsed time to seconds
        IFS=: read -r min sec <<< "$elapsed_time"
        elapsed_time_sec=$(echo "$min * 60 + $sec" | bc)

        # Accumulate totals
        total_user_time=$(echo "$total_user_time + $user_time" | bc)
        total_system_time=$(echo "$total_system_time + $system_time" | bc)
        total_elapsed_time=$(echo "$total_elapsed_time + $elapsed_time_sec" | bc)
        total_max_memory=$(echo "$total_max_memory + $max_memory" | bc)

        # Append individual run stats to the results file
        echo "$i, $user_time, $system_time, $elapsed_time_sec, $max_memory" >> $RESULT_FILE
    done

    # Calculate averages
    avg_user_time=$(echo "$total_user_time / $NUM_RUNS" | bc -l)
    avg_system_time=$(echo "$total_system_time / $NUM_RUNS" | bc -l)
    avg_elapsed_time=$(echo "$total_elapsed_time / $NUM_RUNS" | bc -l)
    avg_max_memory=$(echo "$total_max_memory / $NUM_RUNS" | bc -l)

    # Append averages to the results file
    echo "Averages, $avg_user_time, $avg_system_time, $avg_elapsed_time, $avg_max_memory" >> $RESULT_FILE

    # Output averages to console
    echo "Results for $PARTITION_NAME:"
    echo "Average User Time (seconds): $avg_user_time"
    echo "Average System Time (seconds): $avg_system_time"
    echo "Average Elapsed Time (seconds): $avg_elapsed_time"
    echo "Average Maximum Memory (kbytes): $avg_max_memory"
    echo "------------------------------------------------------"
done