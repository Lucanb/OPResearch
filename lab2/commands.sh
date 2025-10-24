#!/bin/bash
set -e

export GUROBI_HOME="$HOME/Gurobi/gurobi1203/linux64"
export PATH="$PATH:$GUROBI_HOME/bin"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib:$LD_LIBRARY_PATH"
export GRB_LICENSE_FILE="$HOME/Gurobi/gurobi.lic"

OUT="results.txt"
DATE=$(date +"%Y-%m-%d %H:%M:%S")

echo "==============================" > "$OUT"
echo " TWO-PHASE SIMPLEX & GUROBI " >> "$OUT"
echo " Started at: $DATE " >> "$OUT"
echo "==============================" >> "$OUT"
echo "" >> "$OUT"

echo "[BUILD] g++ two_phase.cpp -O2 -std=c++17 -o two_phase" | tee -a "$OUT"
g++ two_phase.cpp -O2 -std=c++17 -o two_phase 2>&1 | tee -a "$OUT"
echo "" | tee -a "$OUT"

for t in a b c d; do
  echo "=== TEST ($t) ===" | tee -a "$OUT"
  ./two_phase < "$t.in" | tee -a "$OUT"
  echo "" >> "$OUT"
  echo "--- GUROBI ($t) ---" | tee -a "$OUT"
  gurobi_cl "problem_${t}.lp" | tee -a "$OUT"
  echo "" >> "$OUT"
  echo "--------------------------------" >> "$OUT"
done

DATE_END=$(date +"%Y-%m-%d %H:%M:%S")
echo "Completed at: $DATE_END" >> "$OUT"
