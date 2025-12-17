#!/usr/bin/env bash
set -e

LP1="a.lp"
LP2="b.lp"

OUTDIR="results"
OUTFILE="$OUTDIR/results.out"

mkdir -p "$OUTDIR"
: > "$OUTFILE"

echo "Problema: $LP1" >> "$OUTFILE"
echo "----------------------------------------" >> "$OUTFILE"
gurobi_cl ResultFile="$OUTDIR/a.sol" "$LP1" >> "$OUTFILE" 2>&1
echo "" >> "$OUTFILE"

echo "Problema: $LP2" >> "$OUTFILE"
echo "----------------------------------------" >> "$OUTFILE"
gurobi_cl ResultFile="$OUTDIR/b.sol" "$LP2" >> "$OUTFILE" 2>&1
echo "" >> "$OUTFILE"

echo "GATA."
echo "Output: $OUTFILE"
echo "Solutii: $OUTDIR/a.sol si $OUTDIR/b.sol"
