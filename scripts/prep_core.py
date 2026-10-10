#!/usr/bin/env python3
"""Make the Nb (or Zr) core-substitution run from the ORIGINAL Ta core run, so that the cell, settings
and starting geometry are identical and only the core atom changes.

  prep_core.py X OUTDIR [--index N]

- Reads CONTCAR (else POSCAR) of $ORIG_TA_CORE_DIR, finds the core Ta atom (the Ta atom with the fewest Ta
  neighbours within 3.3 A; bilayer Ta atoms have several), replaces it by X, and writes OUTDIR/POSCAR with the
  species regrouped (selective-dynamics flags are kept).
- Copies INCAR and KPOINTS unchanged from the original run.
- Builds the POTCAR with scripts/make_potcar.py (original variants for Cu/Li/Ta).
Use --index N (1-based atom number in the original file) to choose the atom explicitly if the automatic
choice is reported as ambiguous. Pure Python 3, no extra packages."""
import os, re, sys, shutil, subprocess

def read_poscar(fn):
    L = open(fn).read().splitlines()
    scale = float(L[1].split()[0])
    cell = [[float(x) * scale for x in L[i].split()[:3]] for i in (2, 3, 4)]
    if all(t.lstrip('-').isdigit() for t in L[5].split()):          # VASP 4 format: no element line
        pot = os.path.join(os.path.dirname(fn), 'POTCAR')
        species = re.findall(r'TITEL\s*=\s*\S+\s+([A-Z][a-z]?)', open(pot).read())
        counts = [int(x) for x in L[5].split()]; i = 6
        print('VASP 4 format: element names taken from POTCAR:', species)
    else:
        species = L[5].split(); counts = [int(x) for x in L[6].split()]; i = 7
    sd = False
    if L[i].strip()[0] in 'sS':
        sd = True; i += 1
    direct = L[i].strip()[0] in 'dD'; i += 1
    sym, pos, flags = [], [], []
    for s, n in zip(species, counts):
        for _ in range(n):
            t = L[i].split(); i += 1
            x = [float(v) for v in t[:3]]
            if not direct:   # convert Cartesian to fractional
                x = cart_to_frac(cell, [v * scale for v in x])
            sym.append(s); pos.append(x); flags.append(' '.join(t[3:6]) if sd else '')
    return L[0], cell, sym, pos, flags, sd

def inv3(m):
    a, b, c = m
    det = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
    return [[(b[1]*c[2]-b[2]*c[1])/det, (a[2]*c[1]-a[1]*c[2])/det, (a[1]*b[2]-a[2]*b[1])/det],
            [(b[2]*c[0]-b[0]*c[2])/det, (a[0]*c[2]-a[2]*c[0])/det, (a[2]*b[0]-a[0]*b[2])/det],
            [(b[0]*c[1]-b[1]*c[0])/det, (a[1]*c[0]-a[0]*c[1])/det, (a[0]*b[1]-a[1]*b[0])/det]]

def cart_to_frac(cell, r):
    inv = inv3(cell)   # rows of cell are lattice vectors: r = f @ cell  ->  f = r @ inv(cell)
    return [sum(r[k] * inv[k][j] for k in range(3)) for j in range(3)]

def expand(tokens):
    out = []
    for t in tokens:
        if '*' in t:
            n, v = t.split('*'); out += [v] * int(n)
        else:
            out.append(t)
    return out

def write_incar(src, dst, orig_sym, new_sym, order):
    """Copy the original INCAR unchanged except: (1) start from scratch (ISTART = 0, ICHARG = 2), since the
    original WAVECAR/CHGCAR are not copied; (2) re-map MAGMOM (if present) to the new atom order."""
    changes = []
    lines = open(src).read().splitlines(); out = []
    for line in lines:
        body = re.split(r'[!#]', line)[0]
        key = body.split('=')[0].strip().upper() if '=' in body else ''
        if key in ('ISTART', 'ICHARG'):
            changes.append(f'removed "{line.strip()}"'); continue
        if key == 'MAGMOM':
            vals = expand(body.split('=', 1)[1].split())
            if len(vals) != len(orig_sym):
                sys.exit(f'MAGMOM has {len(vals)} entries but the cell has {len(orig_sym)} atoms; edit by hand')
            new = [vals[i] for s in order for i in range(len(new_sym)) if new_sym[i] == s]
            comp = []; i = 0
            while i < len(new):
                j = i
                while j < len(new) and new[j] == new[i]: j += 1
                comp.append(f'{j-i}*{new[i]}'); i = j
            line = 'MAGMOM = ' + ' '.join(comp); changes.append('MAGMOM re-mapped to the new atom order')
        out.append(line)
    out += ['ISTART = 0', 'ICHARG = 2']
    open(dst, 'w').write('\n'.join(out) + '\n')
    print('INCAR copied from the original run;', '; '.join(changes) if changes else 'no changes', '+ ISTART = 0, ICHARG = 2')

def dist(cell, f1, f2):
    d = [f1[j] - f2[j] for j in range(3)]; d = [x - round(x) for x in d]
    v = [sum(d[k] * cell[k][j] for k in range(3)) for j in range(3)]
    return sum(x * x for x in v) ** 0.5

def main():
    X, out = sys.argv[1], sys.argv[2]
    idx = int(sys.argv[sys.argv.index('--index') + 1]) - 1 if '--index' in sys.argv else None
    src = os.environ['ORIG_TA_CORE_DIR']
    fn = os.path.join(src, 'CONTCAR') if os.path.isfile(os.path.join(src, 'CONTCAR')) else os.path.join(src, 'POSCAR')
    title, cell, sym, pos, flags, sd = read_poscar(fn)
    ta = [i for i, s in enumerate(sym) if s == 'Ta']
    if idx is None:
        nTa = {i: sum(1 for j in ta if j != i and dist(cell, pos[i], pos[j]) < 3.3) for i in ta}
        ranked = sorted(nTa, key=lambda i: nTa[i])
        idx = ranked[0]
        print(f'Ta atoms: {len(ta)}; core candidate = atom {idx+1} with {nTa[idx]} Ta neighbours '
              f'(next: atom {ranked[1]+1 if len(ranked) > 1 else "-"} with {nTa[ranked[1]] if len(ranked) > 1 else "-"})')
        if len(ranked) > 1 and nTa[ranked[1]] == nTa[idx]:
            sys.exit('AMBIGUOUS core Ta: two candidates with the same neighbour count. Re-run with --index N.')
    elif sym[idx] != 'Ta':
        sys.exit(f'atom {idx+1} is {sym[idx]}, not Ta')
    near = sorted((round(dist(cell, pos[idx], pos[j]), 2), sym[j]) for j in range(len(sym)) if j != idx)[:12]
    print('neighbours of the replaced atom (A):', near)
    orig_sym = list(sym)
    sym[idx] = X
    order = []
    for s in sym:
        if s not in order: order.append(s)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'POSCAR'), 'w') as f:
        f.write(f'{title} | core Ta (atom {idx+1}) -> {X}\n1.0\n')
        for v in cell: f.write('  %.10f %.10f %.10f\n' % tuple(v))
        f.write('  ' + ' '.join(order) + '\n  ' + ' '.join(str(sym.count(s)) for s in order) + '\n')
        if sd: f.write('Selective dynamics\n')
        f.write('Direct\n')
        for s in order:
            for i in range(len(sym)):
                if sym[i] == s:
                    f.write('  %.10f %.10f %.10f %s\n' % (*pos[i], flags[i]))
    shutil.copy(os.path.join(src, 'KPOINTS'), os.path.join(out, 'KPOINTS'))
    write_incar(os.path.join(src, 'INCAR'), os.path.join(out, 'INCAR'), orig_sym, sym, order)
    here = os.path.dirname(os.path.abspath(__file__))
    subprocess.check_call([sys.executable, os.path.join(here, 'make_potcar.py'), os.path.join(out, 'POTCAR')] + order)
    open(os.path.join(out, 'SOURCE.txt'), 'w').write(f'geometry from {fn}\nreplaced atom {idx+1} (Ta) by {X}\nINCAR/KPOINTS copied from {src}\n')
    print(f'{out}: ready')

if __name__ == '__main__':
    main()
