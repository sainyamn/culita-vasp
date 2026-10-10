#!/bin/bash
# Core substitution (Nb, Zr, ...) in the same cell as the original Ta core run, with its INCAR/KPOINTS.
#   bash scripts/run_core.sh [X ...]
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; source "$ROOT/config.sh"; source "$ROOT/scripts/vasp_step.sh"
ELS="$*"; [ -z "$ELS" ] && ELS="$CORE_ELEMENTS"
NK=$(awk 'NR==4{print $1*$2*$3}' "$ORIG_TA_CORE_DIR/KPOINTS" 2>/dev/null)
EXE="$VASP_STD"; [ "$NK" = "1" ] && EXE="$VASP_GAM"
FAIL=0
for X in $ELS; do
  d="$ROOT/core/$X"; echo "=== core $X"
  if step_done "$d"; then echo "  already done"; continue; fi
  if [ ! -f "$d/POSCAR" ]; then
    python3 "$ROOT/scripts/prep_core.py" "$X" "$d" || { echo "  could not prepare core/$X"; FAIL=1; continue; }
  elif [ -s "$d/CONTCAR" ]; then cp "$d/CONTCAR" "$d/POSCAR"; fi     # resume after a wall-time stop
  run_step "$d" "$EXE" || { echo "  core/$X needs attention: $(cat "$d/FAILED.txt" 2>/dev/null)"; FAIL=1; }
done
python3 "$ROOT/scripts/collect.py" > /dev/null
exit $FAIL
