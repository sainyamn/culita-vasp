#!/bin/bash
# Core substitution of X (Nb, Zr, ...) for one Cu atom in the core of the embedded Cu3Li precipitate.
#   bash scripts/run_core.sh [X ...]
# Uses the 459-atom cell shipped in core_inputs/ (relaxed DFT geometry). The host, Ta and every X are relaxed here with
# the same INCAR, so E_subs = E(X) - E(host) - mu_X + mu_Cu is self-contained and Ta is the like-for-like reference.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; source "$ROOT/config.sh"; source "$ROOT/scripts/vasp_step.sh"
ELS="$*"
if [ -z "$ELS" ]; then ELS="host Ta $CORE_ELEMENTS"; fi   # culita submits host, Ta, Nb, Zr as separate jobs
ELS=$(echo $ELS | tr ' ' '\n' | awk '!seen[$0]++' | tr '\n' ' ')
EXE="$VASP_GAM"; [ -z "$EXE" ] && EXE="$VASP_STD"       # Gamma-only cell
echo "core substitution runs: $ELS"
FAIL=0
for X in $ELS; do
  d="$ROOT/core/$X"; echo "=== core $X"
  if step_done "$d"; then echo "  already done"; continue; fi
  if [ ! -f "$d/POSCAR" ]; then
    [ -f "$ROOT/core_inputs/$X/POSCAR" ] || { echo "  no POSCAR for $X (core_inputs/$X)"; FAIL=1; continue; }
    mkdir -p "$d"; cp "$ROOT/core_inputs/$X/POSCAR" "$ROOT/core_inputs/INCAR" "$ROOT/core_inputs/KPOINTS" "$d/"
    python3 "$ROOT/scripts/make_potcar.py" "$d/POTCAR" $(sed -n 6p "$d/POSCAR") || { FAIL=1; continue; }
  elif [ -s "$d/CONTCAR" ]; then cp "$d/CONTCAR" "$d/POSCAR"; fi     # resume after a wall-time stop
  run_step "$d" "$EXE" || { echo "  core/$X needs attention: $(cat "$d/FAILED.txt" 2>/dev/null)"; FAIL=1; }
done
python3 "$ROOT/scripts/collect.py" > /dev/null
exit $FAIL
