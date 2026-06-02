#!/bin/bash
# bench_test2.sh: Benchmark script for parallel weighted Jaccard clustering
# Usage: ./bench_test2.sh [fasta_file] [threshold] [runs]

FASTA=${1:-test2.fasta}
THRESHOLD=${2:-0.85}
RUNS=${3:-3}
EXE="./jaccard_cluster_test"

cd ~/expr2026_baseline/src

if [ ! -f "$FASTA" ]; then
    echo "Error: File $FASTA not found!"
    exit 1
fi

if [ ! -f "$EXE" ]; then
    echo "Error: Executable $EXE not found! Please compile wj.cpp first."
    exit 1
fi

echo "=========================================================="
echo "  Benchmarking weighted Jaccard clustering on $FASTA"
echo "  Threshold: $THRESHOLD | Runs: $RUNS"
echo "=========================================================="
printf "%-8s | " "Threads"
for ((r=1; r<=RUNS; r++)); do
    printf "%-10s | " "Run $r"
done
printf "%-10s | %-8s\n" "Average" "Speedup"
echo "----------------------------------------------------------"

# To store the 1-thread average time for speedup calculation
avg_1t=0.0

for t in 1 2 4 8 16 32 64; do
    export OMP_NUM_THREADS=$t
    sum=0.0
    
    printf "%-8d | " "$t"
    
    for ((r=1; r<=RUNS; r++)); do
        # Run and extract real time
        elapsed=$( { time -p $EXE $FASTA $THRESHOLD > /dev/null; } 2>&1 | grep '^real' | awk '{print $2}' )
        sum=$(echo "$sum $elapsed" | awk '{print $1 + $2}')
        printf "%-10.2f | " "$elapsed"
    done
    
    # Calculate average
    avg=$(echo "$sum $RUNS" | awk '{printf "%.3f", $1 / $2}')
    
    # Calculate speedup
    if [ "$t" -eq 1 ]; then
        avg_1t=$avg
        speedup="1.00x"
    else
        speedup=$(echo "$avg_1t $avg" | awk '{if ($2 > 0) printf "%.2fx", $1 / $2; else printf "0.00x"}')
    fi
    
    printf "%-10.3f | %-8s\n" "$avg" "$speedup"
done
echo "=========================================================="
