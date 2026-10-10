#!/bin/bash
# run_step DIR "VASP COMMAND": run VASP in DIR, check the result, apply an automatic fix and retry (max 4 attempts).
# Returns 0 when the step finished correctly, 1 when it needs a human (reason in DIR/FAILED.txt and DIR/FIXES.txt).
run_step() {
  local d="$1" exe="$2" i
  rm -f "$d/FAILED.txt"
  for i in 1 2 3 4 5; do
    (cd "$d" && $exe > vasp.log 2>&1)
    if python3 "$ROOT/culita.py" _check "$d"; then
      echo "  $(basename "$d"): done  $(grep 'energy(sigma->0)' "$d/OUTCAR" | tail -1 | awk '{print "E0 =",$NF}')"
      rm -f "$d/FAILED.txt"; return 0
    fi
    python3 "$ROOT/culita.py" _fix "$d" || return 1
  done
  return 1
}
step_done() {   # finished earlier (re-runs skip it)
  python3 "$ROOT/culita.py" _check "$1" >/dev/null 2>&1
}
