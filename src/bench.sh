#!/bin/bash
cd ~/expr2026_baseline/src

for t in 1 2 4 8 16 32 64; do
  echo "=== Threads: $t ==="
  export OMP_NUM_THREADS=$t
  { time -p ./jaccard_cluster_test $t test.fasta 0.85 > test_output_t${t}.txt; } 2> time_t${t}.txt
  cat time_t${t}.txt
  head -c 64 test_output_t${t}.txt | md5sum
done

echo "--- REFERENCE ---"
md5sum result1.txt
