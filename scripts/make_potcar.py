#!/usr/bin/env python3
"""Build a POTCAR for a list of species, choosing the PAW variant of each element in this order:
  1. the variant used in an original run (read from the TITEL lines of that run's POTCAR),
  2. the default below (standard Li, Cu, Ta as stated in the paper's supplement: 1, 11, 5 valence electrons).
Usage: make_potcar.py OUT_POTCAR El1 El2 ...      (species in POSCAR order)
Environment: POTCAR_LIB, ORIG_TA_CORE_DIR, ORIG_HOST_DIR, SCREEN_DIR_<El>
Writes a log line per species to stdout (element -> variant, source)."""
import os, re, sys

DEFAULT = {'Cu': 'Cu', 'Li': 'Li', 'Ta': 'Ta', 'O': 'O',
           # element-screen elements: Materials-Project-style choices, used only if the original runs are not available
           'Zr': 'Zr_sv', 'Nb': 'Nb_pv', 'V': 'V_pv', 'Mo': 'Mo_pv', 'W': 'W_pv', 'Ag': 'Ag'}

def variants_in(potcar):
    """Return {element: variant} from the TITEL lines of an existing POTCAR."""
    out = {}
    if not potcar or not os.path.isfile(potcar):
        return out
    for line in open(potcar, errors='replace'):
        m = re.match(r'\s*TITEL\s*=\s*PAW\S*\s+(\S+)', line)
        if m:
            var = m.group(1); el = re.match(r'[A-Z][a-z]?', var).group(0)
            out.setdefault(el, var)
    return out

def choose(el):
    for key in ('SCREEN_DIR_' + el, 'ORIG_TA_CORE_DIR', 'ORIG_HOST_DIR'):
        d = os.environ.get(key, '')
        v = variants_in(os.path.join(d, 'POTCAR')) if d else {}
        if el in v:
            return v[el], f'from {key}/POTCAR'
    return DEFAULT.get(el, el), 'default (original run not available - CHECK)'

def main():
    out, els = sys.argv[1], sys.argv[2:]
    lib = os.environ.get('POTCAR_LIB', '')
    blocks = []
    for el in els:
        var, src = choose(el)
        f = os.path.join(lib, var, 'POTCAR')
        if not os.path.isfile(f):
            sys.exit(f'ERROR: {f} not found (set POTCAR_LIB in config.sh)')
        blocks.append(open(f).read()); print(f'POTCAR  {el:2s} -> {var:8s} ({src})')
    open(out, 'w').write(''.join(blocks))

if __name__ == '__main__':
    main()
