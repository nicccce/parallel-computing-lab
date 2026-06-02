#!/bin/bash
cd ~/expr2026_baseline/src

echo "===== BASELINE ====="
for t in 1 2 4 8 16 32 64; do
  export OMP_NUM_THREADS=$t
  elapsed=$( { time -p ./jaccard_baseline test.fasta 0.85 > /dev/null; } 2>&1 | grep '^real' | awk '{print $2}')
  echo "BASELINE Threads=$t Time=${elapsed}"
done

echo "===== PHASE1 ====="
for t in 1 2 4 8 16 32 64; do
  export OMP_NUM_THREADS=$t
  elapsed=$( { time -p ./jaccard_cluster_test $t test.fasta 0.85 > /dev/null; } 2>&1 | grep '^real' | awk '{print $2}')
  echo "PHASE1 Threads=$t Time=${elapsed}"
done
