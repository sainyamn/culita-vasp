#!/bin/bash
# Bulk chemical potentials: ISIF=3 relaxation twice, static (tetrahedron), and a denser-k static as the k-point test.
#   bash scripts/run_bulk.sh [El ...]
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; source "$ROOT/config.sh"; source "$ROOT/scripts/vasp_step.sh"
ELS="$*"; [ -z "$ELS" ] && ELS=$(ls "$ROOT/bulk")
FAIL=0
for el in $ELS; do
  d="$ROOT/bulk/$el"; echo "=== bulk $el"
  python3 "$ROOT/scripts/make_potcar.py" "$d/POTCAR" "$el" || { echo "  POTCAR for $el could not be built"; FAIL=1; continue; }
  awk 'NR==4{for(i=1;i<=3;i++) printf "%d ", int($i*1.35+0.999); print ""; next} {print}' "$d/KPOINTS" > "$d/KPOINTS.dense"
  for step in relax1 relax2 static static_dense; do
    w="$d/$step"
    if step_done "$w"; then echo "  $step: already done"; continue; fi
    mkdir -p "$w"; cp "$d/KPOINTS" "$d/POTCAR" "$w/"
    case $step in
      relax1) src="$d/POSCAR"; inc=INCAR.relax ;;
      relax2) src="$d/relax1/CONTCAR"; inc=INCAR.relax ;;
      static) src="$d/relax2/CONTCAR"; inc=INCAR.static ;;
      static_dense) src="$d/relax2/CONTCAR"; inc=INCAR.static; cp "$d/KPOINTS.dense" "$w/KPOINTS" ;;
    esac
    [ -f "$w/INCAR" ] || cp "$d/$inc" "$w/INCAR"          # keep auto-fixed INCAR on re-runs
    [ -s "$w/CONTCAR" ] && [ "$step" != static ] && [ "$step" != static_dense ] && cp "$w/CONTCAR" "$w/POSCAR" || cp "$src" "$w/POSCAR"
    run_step "$w" "$VASP_STD" || { echo "  $el/$step needs attention: $(cat "$w/FAILED.txt" 2>/dev/null)"; FAIL=1; break; }
  done
done
python3 "$ROOT/scripts/collect.py" > /dev/null
exit $FAIL
