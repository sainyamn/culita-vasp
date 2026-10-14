#!/usr/bin/env python3
"""Build a POTCAR for a list of species with fixed PAW variants:
  Cu, Li, Ta standard (as stated in the supplement: 11, 1, 5 valence electrons);
  Zr_sv, Nb_pv, V_pv, Mo_pv, W_pv, Ag for the element screen.
Usage: make_potcar.py OUT_POTCAR El1 El2 ...      (species in POSCAR order)
Environment: POTCAR_LIB"""
import os, sys

VARIANT = {'Cu': 'Cu', 'Li': 'Li', 'Ta': 'Ta', 'O': 'O',
           'Zr': 'Zr_sv', 'Nb': 'Nb_pv', 'V': 'V_pv', 'Mo': 'Mo_pv', 'W': 'W_pv', 'Ag': 'Ag'}

def main():
    out, els = sys.argv[1], sys.argv[2:]
    lib = os.environ.get('POTCAR_LIB', '')
    blocks = []
    for el in els:
        var = VARIANT.get(el, el)
        f = os.path.join(lib, var, 'POTCAR')
        if not os.path.isfile(f):
            sys.exit(f'ERROR: {f} not found (set POTCAR_LIB in config.sh)')
        blocks.append(open(f).read()); print(f'POTCAR  {el:2s} -> {var}')
    open(out, 'w').write(''.join(blocks))

if __name__ == '__main__':
    main()
