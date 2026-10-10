#!/bin/bash
# Core substitution of X (Nb, Zr, ...) for one Cu atom in the core of the embedded Cu3Li precipitate.
#   bash scripts/run_core.sh [X ...]
# Two modes:
#  * ORIGINAL (ORIG_TA_CORE_DIR set): X replaces the core Ta of our original Ta-core run; its INCAR/KPOINTS/geometry are
#    reused, and the original host energy E1 is used.
#  * PACKAGED (no original runs available): the 459-atom cell shipped in core_inputs/ (relaxed DFT geometry); the host,
#    Ta and every X are relaxed here with the same INCAR, so E_subs = E(X) - E(host) - mu_X + mu_Cu is self-contained and
#    the Ta value can be compared with the published +0.81 eV.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; source "$ROOT/config.sh"; source "$ROOT/scripts/vasp_step.sh"
ELS="$*"; [ -z "$ELS" ] && ELS="$CORE_ELEMENTS"
FAIL=0
if [ -n "$ORIG_TA_CORE_DIR" ] && [ -d "$ORIG_TA_CORE_DIR" ]; then
  MODE=original
  NK=$(awk 'NR==4{print $1*$2*$3}' "$ORIG_TA_CORE_DIR/KPOINTS" 2>/dev/null)
else
  MODE=packaged; NK=1
  [ -z "$*" ] && ELS="host Ta $ELS"        # host and Ta are always needed in this mode (culita submits them as separate jobs)
  ELS=$(echo $ELS | tr ' ' '\n' | awk '!seen[$0]++' | tr '\n' ' ')
fi
EXE="$VASP_STD"; [ "$NK" = "1" ] && EXE="$VASP_GAM"
echo "core substitution mode: $MODE  (runs: $ELS)"
for X in $ELS; do
  d="$ROOT/core/$X"; echo "=== core $X"
  if step_done "$d"; then echo "  already done"; continue; fi
  if [ ! -f "$d/POSCAR" ]; then
    if [ $MODE = original ]; then
      python3 "$ROOT/scripts/prep_core.py" "$X" "$d" || { echo "  could not prepare core/$X"; FAIL=1; continue; }
    else
      [ -f "$ROOT/core_inputs/$X/POSCAR" ] || { echo "  no packaged POSCAR for $X (core_inputs/$X)"; FAIL=1; continue; }
      mkdir -p "$d"; cp "$ROOT/core_inputs/$X/POSCAR" "$ROOT/core_inputs/INCAR" "$ROOT/core_inputs/KPOINTS" "$d/"
      python3 "$ROOT/scripts/make_potcar.py" "$d/POTCAR" $(sed -n 6p "$d/POSCAR") || { FAIL=1; continue; }
    fi
  elif [ -s "$d/CONTCAR" ]; then cp "$d/CONTCAR" "$d/POSCAR"; fi     # resume after a wall-time stop
  run_step "$d" "$EXE" || { echo "  core/$X needs attention: $(cat "$d/FAILED.txt" 2>/dev/null)"; FAIL=1; }
done
python3 "$ROOT/scripts/collect.py" > /dev/null
exit $FAIL
