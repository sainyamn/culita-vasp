#!/usr/bin/env python3
"""Collect the results into results.txt, with convergence checks.
  E_subs(X) = E(core X) - E1(host) - mu_X + mu_Cu        (Eq. S6)

Checks written to results.txt (PASS / CHECK):
  bulk  - k-points : |mu(dense mesh) - mu(standard mesh)| < 1 meV/atom
        - cell     : |pressure| of the final static run < 2 kB (cell fully relaxed)
        - relax    : second ISIF=3 relaxation reached the force criterion
        - controls : Cu, Li, Ta within 10 meV/atom of the paper's mu
  core  - ionic    : 'reached required accuracy' and the largest remaining force
        - field    : which OUTCAR energy (energy(sigma->0) or TOTEN) reproduces the original Ta value of +0.81 eV
                     from the original OUTCARs; that field is then used for Nb and Zr.
Pure Python 3."""
import os, re, glob

PAPER_MU = {'Cu': -3.728, 'Li': -2.051, 'Ta': -11.865}
PAPER_TA_SUBS = 0.81

def read_outcar(outcar):
    r = dict(e0=None, toten=None, n=None, p=None, fmax=None, conv=False)
    if not os.path.isfile(outcar):
        return r
    L = open(outcar, errors='replace').read().splitlines()
    for i, line in enumerate(L):
        if 'NIONS' in line and r['n'] is None:
            r['n'] = int(line.split('=')[-1])
        m = re.search(r'energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)', line)
        if m: r['e0'] = float(m.group(1))
        m = re.search(r'free\s+energy\s+TOTEN\s*=\s*(-?\d+\.\d+)', line)
        if m: r['toten'] = float(m.group(1))
        m = re.search(r'external pressure\s*=\s*(-?\d+\.\d+)', line)
        if m: r['p'] = float(m.group(1))
        if 'TOTAL-FORCE' in line and r['n']:
            f = []
            for l2 in L[i + 2:i + 2 + r['n']]:
                t = l2.split()
                if len(t) >= 6:
                    f.append(sum(float(x) ** 2 for x in t[3:6]) ** 0.5)
            if f: r['fmax'] = max(f)
        if 'reached required accuracy' in line:
            r['conv'] = True
    return r

def flag(ok):
    return 'PASS ' if ok else 'CHECK'

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = ['BULK CHEMICAL POTENTIALS  (mu = E(sigma->0)/N of the final static run)',
       f'{"El":3s} {"mu (eV/atom)":>13s} {"dense-k diff":>13s} {"p (kB)":>8s} {"relax":>6s} {"paper mu":>9s}  checks / potentials']
mu = {}
for d in sorted(glob.glob(os.path.join(root, 'bulk', '*'))):
    el = os.path.basename(d)
    s = read_outcar(os.path.join(d, 'static', 'OUTCAR')); sd = read_outcar(os.path.join(d, 'static_dense', 'OUTCAR'))
    r2 = read_outcar(os.path.join(d, 'relax2', 'OUTCAR'))
    if s['e0'] is None:
        out.append(f'{el:3s} not run yet'); continue
    mu[el] = s['e0'] / s['n']
    dk = sd['e0'] / sd['n'] - mu[el] if sd['e0'] is not None else None
    checks = [f'k {flag(dk is not None and abs(dk) < 1e-3)}', f'cell {flag(s["p"] is not None and abs(s["p"]) < 2)}', f'relax {flag(r2["conv"])}']
    if el in PAPER_MU:
        checks.append(f'control {flag(abs(mu[el] - PAPER_MU[el]) < 0.010)}')
    pot = ' '.join(re.findall(r'TITEL\s*=\s*\S+\s+(\S+)', open(os.path.join(d, 'POTCAR')).read())) if os.path.isfile(os.path.join(d, 'POTCAR')) else ''
    out.append(f'{el:3s} {mu[el]:13.5f} {(dk*1000 if dk is not None else float("nan")):10.2f} meV {(s["p"] if s["p"] is not None else float("nan")):8.2f} '
               f'{str(r2["conv"]):>6s} {PAPER_MU.get(el, float("nan")):9.3f}  {", ".join(checks)}  [{pot}]')

out += ['', 'CORE SUBSTITUTION  E_subs = E(core X) - E1(host) - mu_X + mu_Cu']
h = read_outcar(os.path.join(os.environ.get('ORIG_HOST_DIR', ''), 'OUTCAR'))
t = read_outcar(os.path.join(os.environ.get('ORIG_TA_CORE_DIR', ''), 'OUTCAR'))
field = 'e0'
if h['e0'] is not None and t['e0'] is not None:
    best = None
    for fld, name in (('e0', 'energy(sigma->0)'), ('toten', 'TOTEN')):
        if h[fld] is None or t[fld] is None: continue
        v = t[fld] - h[fld] - PAPER_MU['Ta'] + PAPER_MU['Cu']
        out.append(f'CHECK original Ta with {name:16s}: E_subs(Ta) = {v:+.3f} eV (paper +{PAPER_TA_SUBS})')
        if best is None or abs(v - PAPER_TA_SUBS) < best[0]:
            best = (abs(v - PAPER_TA_SUBS), fld, name)
    field = best[1]
    out.append(f'  -> using {best[2]} for Nb/Zr  [{flag(best[0] < 0.02)} reproduces +0.81 within 0.02 eV]')
    if 'Ta' in mu and 'Cu' in mu:
        out.append(f'  original Ta with the NEW bulk mu: E_subs(Ta) = {t[field] - h[field] - mu["Ta"] + mu["Cu"]:+.3f} eV')
else:
    out.append('original host / Ta-core OUTCAR not found: check ORIG_HOST_DIR and ORIG_TA_CORE_DIR in config.sh')
for d in sorted(glob.glob(os.path.join(root, 'core', '*'))):
    X = os.path.basename(d); c = read_outcar(os.path.join(d, 'OUTCAR'))
    if c[field] is None or h[field] is None or X not in mu:
        out.append(f'{X}: not available yet'); continue
    fm = f'{c["fmax"]:.3f} eV/A' if c['fmax'] is not None else 'n/a'
    out.append(f'{X}: relaxation {flag(c["conv"])} (largest remaining force {fm}); E = {c[field]:.4f} eV')
    for lab, mcu in (('paper mu_Cu', PAPER_MU['Cu']), ('new mu_Cu', mu.get('Cu'))):
        if mcu is not None:
            out.append(f'    E_subs({X}) = {c[field] - h[field] - mu[X] + mcu:+.3f} eV   [{lab}, new mu_{X}]')
    if 'Ta' in mu and 'Cu' in mu:
        out.append(f'    E_subs({X}) - E_subs(Ta) = {(c[field] - mu[X]) - (t[field] - mu["Ta"]):+.3f} eV  (all new mu; independent of mu_Cu and E1)')
txt = '\n'.join(out) + '\n'
open(os.path.join(root, 'results.txt'), 'w').write(txt); print(txt)
