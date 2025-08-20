#!/bin/bash

# Script to monitor R processes during searchlight execution

echo "========================================="
echo "R Process Monitor for Searchlight"
echo "========================================="
echo
echo "This script will monitor R processes every 0.1 seconds"
echo "Run your searchlight in another terminal and watch here"
echo
echo "Press Ctrl+C to stop monitoring"
echo
echo "Time       | R Processes | Delta"
echo "-----------|-------------|-------"

prev_count=0
while true; do
    # Count R processes
    count=$(ps aux | grep '[R] ' | grep -v grep | wc -l | tr -d ' ')
    
    # Calculate delta
    delta=$((count - prev_count))
    delta_str=""
    if [ $delta -gt 0 ]; then
        delta_str="+$delta"
    elif [ $delta -lt 0 ]; then
        delta_str="$delta"
    else
        delta_str="="
    fi
    
    # Print with timestamp
    printf "%s | %11s | %s\n" "$(date +%H:%M:%S.%N | cut -c1-11)" "$count" "$delta_str"
    
    # Update previous count
    prev_count=$count
    
    # Short sleep
    sleep 0.1
done