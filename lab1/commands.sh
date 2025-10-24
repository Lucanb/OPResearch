#!/bin/bash
set -e

SRC="simplex.cpp"
EXE="simplex"
OUT="results.txt"

export GUROBI_HOME="$HOME/Gurobi/gurobi1203/linux64"
export PATH="$PATH:$GUROBI_HOME/bin"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib:$LD_LIBRARY_PATH"
export GRB_LICENSE_FILE="$HOME/Gurobi/gurobi.lic"

DATE=$(date +"%Y-%m-%d %H:%M:%S")
echo "==============================" > "$OUT"
echo " SIMPLEX ALGORITHM + GUROBI " >> "$OUT"
echo " Started at: $DATE " >> "$OUT"
echo "==============================" >> "$OUT"
echo "" >> "$OUT"

echo "[BUILD] g++ $SRC -O2 -std=c++17 -o $EXE" | tee -a "$OUT"
g++ "$SRC" -O2 -std=c++17 -o "$EXE" 2>&1 | tee -a "$OUT"
echo "" >> "$OUT"

for infile in input_*.txt; do
  [ -e "$infile" ] || continue
  base="${infile#input_}"
  name="${base%.txt}"
  outfile="output_${name}.txt"
  lpfile="problem_${name}.lp"
  
  echo "=== RUNNING SIMPLEX ON $infile ===" | tee -a "$OUT"
  ./"$EXE" "$infile" "$outfile" | tee -a "$OUT"

  echo "" >> "$OUT"
  echo "--- SIMPLEX OUTPUT: $outfile ---" | tee -a "$OUT"
  cat "$outfile" | tee -a "$OUT"
  echo "" >> "$OUT"

  if [ -f "$lpfile" ]; then
    echo "--- GUROBI RESULT ($lpfile) ---" | tee -a "$OUT"
    gurobi_cl "$lpfile" | tee -a "$OUT"
  else
    echo "(Nu există fișier LP pentru $infile → căutat $lpfile)" | tee -a "$OUT"
  fi

  echo "--------------------------------" >> "$OUT"
  echo "" >> "$OUT"
done

DATE_END=$(date +"%Y-%m-%d %H:%M:%S")
echo "Completed at: $DATE_END" >> "$OUT"
echo "Rezultatele complete sunt salvate în $OUT"
