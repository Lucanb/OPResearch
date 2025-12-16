#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
OUT="results/results.out"

: > "$OUT"

echo "Problema: a.lp" >> "$OUT"
echo "----------------------------------------" >> "$OUT"
gurobi_cl ResultFile=results/a.sol a.lp >> "$OUT" 2>&1
echo "" >> "$OUT"

echo "Problema: b.lp" >> "$OUT"
echo "----------------------------------------" >> "$OUT"
gurobi_cl ResultFile=results/b.sol b.lp >> "$OUT" 2>&1
echo "" >> "$OUT"

echo "GATA. Output: $OUT"
echo "Solutii: results/a.sol si results/b.sol"
