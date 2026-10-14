#!/usr/bin/env python3
"""culita - run the remaining VASP calculations for the Cu-Ta-Li bilayer-complexion letter.

  ./culita setup      guided configuration (auto-detects VASP, PAW library, scheduler)
  ./culita check      pre-flight checks only (nothing is submitted)
  ./culita run        submit the jobs (Slurm / PBS) or run them here   [--bulk | --core] [--local]
  ./culita status     progress of every calculation, with VASP errors explained
  ./culita results    write and show results.txt (chemical potentials, E_subs, convergence checks)
  ./culita pack       bundle results to send back (results.txt, CONTCARs, OSZICARs, bulk OUTCARs)
  ./culita selftest   end-to-end test with a fake VASP (no VASP or cluster time needed)

Only the Python 3 standard library is used (Python >= 3.6)."""
import configparser, datetime, glob, json, os, re, shutil, subprocess, sys, tarfile, tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
CFG = os.path.join(ROOT, 'culita.cfg')
STATE = os.path.join(ROOT, '.culita_state.json')
BULK = ['Zr', 'Nb', 'V', 'Mo', 'W', 'Ag', 'Cu', 'Li', 'Ta']
SCREEN = ['Zr', 'Nb', 'V', 'Mo', 'W', 'Ag']
VARIANTS = {'Cu': ['Cu'], 'Li': ['Li', 'Li_sv'], 'Ta': ['Ta', 'Ta_pv'], 'Zr': ['Zr_sv', 'Zr'], 'Nb': ['Nb_pv', 'Nb_sv'],
            'V': ['V_pv', 'V_sv', 'V'], 'Mo': ['Mo_pv', 'Mo_sv', 'Mo'], 'W': ['W_pv', 'W_sv', 'W'], 'Ag': ['Ag', 'Ag_pv']}

# ------------------------------------------------------------------ helpers
C = {'ok': '\033[32m', 'warn': '\033[33m', 'err': '\033[31m', 'b': '\033[1m', 'end': '\033[0m'}
if not sys.stdout.isatty():
    C = {k: '' for k in C}

def say(msg, kind=None):
    pre = {'ok': '[ OK ] ', 'warn': '[WARN] ', 'err': '[FAIL] ', None: ''}[kind]
    col = C.get(kind, '') if kind else ''
    print(col + pre + C['end'] * bool(col) + msg)

def die(msg, hint=None):
    say(msg, 'err')
    if hint:
        print('       -> ' + hint)
    sys.exit(1)

def ask(prompt, default=''):
    if not sys.stdin.isatty():
        return default
    v = input('%s%s%s [%s]: ' % (C['b'], prompt, C['end'], default)).strip()
    return v or default

def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              universal_newlines=True, timeout=60).stdout
    except Exception as e:
        return str(e)

def load_cfg():
    c = configparser.ConfigParser(); c.optionxform = str
    if os.path.isfile(CFG):
        c.read(CFG)
    for s in ('vasp', 'paths', 'scheduler', 'run'):
        if not c.has_section(s): c.add_section(s)
    return c

def save_cfg(c):
    with open(CFG, 'w') as f:
        f.write('# written by ./culita setup - edit by hand or re-run setup\n'); c.write(f)
    write_config_sh(c)

def write_config_sh(c):
    g = lambda s, k, d='': c.get(s, k, fallback=d)
    lines = ['#!/bin/bash', '# generated from culita.cfg by ./culita - do not edit (edit culita.cfg instead)']
    if g('vasp', 'modules'):
        lines.append(g('vasp', 'modules'))
    lines += ['export VASP_STD=%s' % sq(g('vasp', 'std')), 'export VASP_GAM=%s' % sq(g('vasp', 'gam') or g('vasp', 'std')),
              'export POTCAR_LIB=%s' % sq(g('paths', 'potcar_lib'))]
    lines.append('export CORE_ELEMENTS=%s' % sq(g('run', 'core_elements', 'Nb Zr')))
    open(os.path.join(ROOT, 'config.sh'), 'w').write('\n'.join(lines) + '\n')

def sq(s):
    return "'" + str(s).replace("'", "'\\''") + "'"

def env_from_cfg(c):
    e = dict(os.environ)
    e.update({'VASP_STD': c.get('vasp', 'std', fallback=''), 'VASP_GAM': c.get('vasp', 'gam', fallback=''),
              'POTCAR_LIB': c.get('paths', 'potcar_lib', fallback='')})
    return e

def core_jobs(c):
    xs = c.get('run', 'core_elements', fallback='Nb Zr').split()
    return ['host', 'Ta'] + [x for x in xs if x not in ('host', 'Ta')]

# ------------------------------------------------------------------ validation
def read_species_counts(fn):
    L = open(fn).read().splitlines()
    if all(t.lstrip('-').isdigit() for t in L[5].split()):
        pot = os.path.join(os.path.dirname(fn), 'POTCAR')
        sp = re.findall(r'TITEL\s*=\s*\S+\s+([A-Z][a-z]?)', open(pot).read()) if os.path.isfile(pot) else []
        return sp, [int(x) for x in L[5].split()]
    return L[5].split(), [int(x) for x in L[6].split()]

def outcar_energy(fn):
    e = None
    if os.path.isfile(fn):
        for line in open(fn, errors='replace'):
            m = re.search(r'energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)', line)
            if m: e = float(m.group(1))
    return e

def check_vasp(c):
    ok = True
    for key in ('std', 'gam'):
        cmd = c.get('vasp', key, fallback='')
        if not cmd:
            if key == 'std': say('VASP launch command (vasp std) not set', 'err'); ok = False
            continue
        exe = [t for t in cmd.split() if 'vasp' in os.path.basename(t)]
        if not exe:
            say('cannot find a vasp executable name in "%s"' % cmd, 'warn'); continue
        found = shutil.which(exe[-1]) or os.path.isfile(exe[-1])
        mod = c.get('vasp', 'modules', fallback='')
        if not found and mod:
            found = 'vasp' in sh('bash -lc %s' % sq(mod + ' >/dev/null 2>&1; command -v ' + exe[-1]))
        if found: say('VASP (%s): %s' % (key, cmd), 'ok')
        else:
            say('"%s" not found on PATH%s; jobs will fail unless it exists on the compute nodes' % (exe[-1], ' (after modules)' if mod else ''), 'warn')
    return ok

def check_potcars(c):
    lib = c.get('paths', 'potcar_lib', fallback='')
    if not lib or not os.path.isdir(lib):
        say('PAW library not found: %r' % lib, 'err'); return False
    ok = True
    for el, vs in VARIANTS.items():
        have = [v for v in vs if os.path.isfile(os.path.join(lib, v, 'POTCAR'))]
        if have: say('PAW %-2s: %s available' % (el, ', '.join(have)), 'ok')
        else:
            say('PAW %-2s: none of %s in %s' % (el, '/'.join(vs), lib), 'err'); ok = False
    return ok

def check_core_inputs(c):
    need = core_jobs(c)
    miss = [x for x in need if not os.path.isfile(os.path.join(ROOT, 'core_inputs', x, 'POSCAR'))]
    if miss:
        say('core cells missing for: %s (core_inputs/)' % ', '.join(miss), 'err'); return False
    say('core runs: 459-atom cell in core_inputs/ - %s will be relaxed with core_inputs/INCAR' % ', '.join(need), 'ok')
    return True

def preflight(c, quiet_header=False):
    if not quiet_header: print(C['b'] + 'Pre-flight checks' + C['end'])
    a = check_vasp(c); b = check_potcars(c); d = check_core_inputs(c)
    sched = c.get('scheduler', 'type', fallback='local')
    say('scheduler: %s' % sched, 'ok')
    return a and b and d

# ------------------------------------------------------------------ setup
def detect_potcar_lib():
    cands = []
    for var in ('VASP_PP_PATH', 'PMG_VASP_PSP_DIR', 'POTCAR_LIB'):
        v = os.environ.get(var)
        if v: cands += [v, os.path.join(v, 'potpaw_PBE'), os.path.join(v, 'potpaw_PBE.54'), os.path.join(v, 'POT_GGA_PAW_PBE')]
    for base in ('/opt', '/usr/local', '/apps', '/packages', '/software', os.path.expanduser('~')):
        cands += glob.glob(os.path.join(base, '*', 'potpaw_PBE*')) + glob.glob(os.path.join(base, '*', '*', 'potpaw_PBE*'))
    for p in cands:
        if os.path.isfile(os.path.join(p, 'Cu', 'POTCAR')) and os.path.isfile(os.path.join(p, 'Ta', 'POTCAR')):
            return p
    return ''

def cmd_setup(args):
    c = load_cfg()
    print(C['b'] + 'culita setup' + C['end'] + ' - press Enter to accept the value in [brackets].\n')
    # VASP
    mods = sh("bash -lc 'module -t avail 2>&1 | grep -i vasp'").split()
    mod_default = c.get('vasp', 'modules', fallback='') or (('module load ' + mods[0]) if mods else '')
    if mods: print('  VASP modules found: ' + ', '.join(mods[:6]))
    c.set('vasp', 'modules', ask('Command(s) to load VASP 5.4.4 (blank if not needed)', mod_default))
    np_ = ask('MPI ranks per job', c.get('run', 'ntasks', fallback='32')); c.set('run', 'ntasks', np_)
    std = shutil.which('vasp_std') or 'vasp_std'; gam = shutil.which('vasp_gam') or 'vasp_gam'
    c.set('vasp', 'std', ask('Command to run vasp_std', c.get('vasp', 'std', fallback='mpirun -np %s %s' % (np_, std))))
    c.set('vasp', 'gam', ask('Command to run vasp_gam (Gamma-only 459-atom cells; same as above if vasp_gam is not available)',
                             c.get('vasp', 'gam', fallback='mpirun -np %s %s' % (np_, gam))))
    # PAW
    c.set('paths', 'potcar_lib', os.path.abspath(os.path.expanduser(ask('PBE PAW library (folder with Cu/POTCAR, Ta/POTCAR, ...)',
                                                                      c.get('paths', 'potcar_lib', fallback='') or detect_potcar_lib()))))
    print('\n  Core substitution: the 459-atom cell in core_inputs/ is used; the host, Ta, Nb and Zr are all relaxed here.')
    c.set('run', 'core_elements', ask('Core substitutions to run', c.get('run', 'core_elements', fallback='Nb Zr')))
    # scheduler
    sched = 'slurm' if shutil.which('sbatch') else ('pbs' if shutil.which('qsub') else 'local')
    sched = ask('Scheduler (slurm / pbs / local)', c.get('scheduler', 'type', fallback=sched)); c.set('scheduler', 'type', sched)
    if sched in ('slurm', 'pbs'):
        c.set('scheduler', 'partition', ask('Partition / queue', c.get('scheduler', 'partition', fallback='')))
        c.set('scheduler', 'account', ask('Account / allocation (blank if none)', c.get('scheduler', 'account', fallback='')))
        c.set('scheduler', 'nodes', ask('Nodes per job', c.get('scheduler', 'nodes', fallback='1')))
        c.set('scheduler', 'time_bulk', ask('Wall time for the bulk job (HH:MM:SS)', c.get('scheduler', 'time_bulk', fallback='06:00:00')))
        c.set('scheduler', 'time_core', ask('Wall time per core-substitution job (HH:MM:SS)', c.get('scheduler', 'time_core', fallback='48:00:00')))
    save_cfg(c)
    say('saved culita.cfg and config.sh\n', 'ok')
    if preflight(c):
        say('all checks passed - next: ./culita run', 'ok')
    else:
        say('fix the items marked FAIL (re-run ./culita setup, or edit culita.cfg and run ./culita check)', 'err')

# ------------------------------------------------------------------ VASP output diagnosis
ERRORS = [
    (r'VERY BAD NEWS|internal error in subroutine IBZKPT|Tetrahedron method fails', 'k-point / symmetry problem', 'for bulk runs try ISYM = 0 in bulk/<El>/INCAR.*; report the message'),
    (r'ZBRENT: fatal error', 'ionic optimiser (ZBRENT) failure', 'usually harmless near convergence: re-run ./culita run (restarts from CONTCAR)'),
    (r'Error EDDDAV|EDDRMM: call to ZHEGV|Sub-Space-Matrix is not hermitian|ZPOTRF|PSSYEVX|CNORMN', 'electronic minimiser failure',
     'add ALGO = Normal (or ALGO = All) to the INCAR of that step and re-run'),
    (r'BRMIX: very serious problems', 'charge-density mixing failure', 're-run; if it repeats, add AMIX = 0.1, BMIX = 0.01 to that INCAR'),
    (r'POTCAR and POSCAR are incompatible|number of potentials on File POTCAR incompatible', 'POTCAR does not match POSCAR', 'check the PAW library path / variants in culita.cfg'),
    (r'out of memory|OOM|Killed|oom-kill|Cannot allocate memory', 'out of memory', 'use more nodes or fewer ranks per node (scheduler settings in culita.cfg)'),
    (r'segmentation fault|SIGSEGV|forrtl: severe', 'VASP crashed (segfault)', 'often stack size: add "ulimit -s unlimited" to the modules line in culita.cfg'),
    (r'command not found|No such file or directory.*vasp', 'VASP executable not found on the compute node', 'fix the modules / VASP command in culita.cfg'),
    (r'DUE TO TIME LIMIT|CANCELLED|walltime', 'job hit the wall-time limit', 'just re-run ./culita run: finished steps are kept, the core run restarts from CONTCAR'),
]

def diagnose(d, finished=False):
    """Return (state, message) for a VASP run folder: done / running / failed / incomplete / not started.
    finished=True: VASP has exited (called right after a run), so never report 'running'."""
    out = os.path.join(d, 'OUTCAR')
    if not os.path.isdir(d) or not os.path.isfile(out):
        return 'not started', ''
    txt = open(out, errors='replace').read()
    logs = ''
    p = os.path.join(d, 'vasp.log')
    if os.path.isfile(p): logs = open(p, errors='replace').read()[-20000:]
    normal = 'General timing and accounting' in txt
    for pat, what, fix in ERRORS:
        if re.search(pat, txt[-50000:] + logs):
            return 'failed', '%s -> %s' % (what, fix)
    if not normal:
        try:
            age = (datetime.datetime.now() - datetime.datetime.fromtimestamp(os.path.getmtime(out))).total_seconds()
        except OSError:
            age = 1e9
        return ('running', 'OUTCAR updated %d s ago' % age) if (age < 900 and not finished) else ('incomplete', 'stopped before normal termination (wall time?) -> re-run ./culita run')
    # electronic convergence of the last ionic step
    osz = os.path.join(d, 'OSZICAR'); nelm = 60
    inc = os.path.join(d, 'INCAR')
    if os.path.isfile(inc):
        m = re.search(r'^\s*NELM\s*=\s*(\d+)', open(inc).read(), re.M | re.I)
        if m: nelm = int(m.group(1))
    if os.path.isfile(osz):
        steps = [l for l in open(osz).read().splitlines() if re.match(r'^\s*(DAV|RMM|CG|SDA|DIA)\s*:', l)]
        last_block = []
        for l in open(osz).read().splitlines():
            if re.match(r'^\s*\d+\s+F=', l): last_block = []
            elif re.match(r'^\s*(DAV|RMM|CG|SDA|DIA)\s*:', l): last_block.append(l)
        if last_block and len(last_block) >= nelm:
            return 'failed', 'electronic SCF did not converge within NELM = %d -> add ALGO = Normal and NELM = 200 to that INCAR, then re-run' % nelm
    relax = re.search(r'^\s*NSW\s*=\s*([1-9]\d*)', open(inc).read(), re.M | re.I) if os.path.isfile(inc) else None
    if relax and 'reached required accuracy' not in txt:
        return 'incomplete', 'ionic relaxation not converged within NSW -> re-run ./culita run (restarts from CONTCAR)'
    return 'done', ''

def cmd_check_step(args):
    """internal: called by the run scripts after each VASP step; non-zero exit stops the chain."""
    d = args[0]; st, msg = diagnose(d, finished=True)
    if st == 'done':
        return 0
    open(os.path.join(d, 'FAILED.txt'), 'w').write('%s: %s\n' % (st, msg))
    say('%s: %s %s' % (os.path.relpath(d, ROOT), st, msg), 'err'); sys.exit(2)


def set_incar(d, settings, note):
    inc = os.path.join(d, 'INCAR'); L = open(inc).read().splitlines()
    keys = {k.upper() for k in settings}
    L = [l for l in L if not (('=' in l) and l.split('=')[0].strip().upper() in keys)]
    L += ['%s = %s   ! culita auto-fix: %s' % (k, v, note) for k, v in settings.items()]
    open(inc, 'w').write('\n'.join(L) + '\n')

def backup(d):
    n = 1
    while os.path.isdir(os.path.join(d, 'attempt%d' % n)): n += 1
    b = os.path.join(d, 'attempt%d' % n); os.makedirs(b)
    for f in ('INCAR', 'POSCAR', 'CONTCAR', 'OUTCAR', 'OSZICAR', 'vasp.log', 'FAILED.txt'):
        if os.path.isfile(os.path.join(d, f)): shutil.copy(os.path.join(d, f), b)
    return n

def restart_from_contcar(d):
    c = os.path.join(d, 'CONTCAR')
    if os.path.isfile(c) and os.path.getsize(c) > 0:
        shutil.copy(c, os.path.join(d, 'POSCAR')); return True
    return False

def cmd_fix_step(args):
    """internal: apply one automatic fix after a failed step (exit 0 = retry, 3 = needs a human)."""
    d = args[0]; st, msg = diagnose(d, finished=True)
    n = backup(d); log = open(os.path.join(d, 'FIXES.txt'), 'a')
    inc = open(os.path.join(d, 'INCAR')).read() if os.path.isfile(os.path.join(d, 'INCAR')) else ''
    def done(text):
        log.write('attempt %d: %s | %s\n' % (n, msg or st, text)); log.close(); say('%s: auto-fix: %s' % (os.path.relpath(d, ROOT), text), 'warn'); sys.exit(0)
    if n > 4:
        log.write('giving up after %d attempts\n' % (n - 1)); log.close(); say('%s: still failing after %d automatic attempts: %s' % (os.path.relpath(d, ROOT), n - 1, msg), 'err'); sys.exit(3)
    if st in ('incomplete', 'running') or 'ZBRENT' in msg or 'wall-time' in msg:
        done('restarted from CONTCAR' if restart_from_contcar(d) else 'restarted from the same POSCAR')
    if 'electronic' in msg or 'SCF did not converge' in msg:
        restart_from_contcar(d)
        if re.search(r'^\s*ALGO\s*=\s*N', inc, re.M | re.I):
            set_incar(d, {'ALGO': 'All', 'NELM': '200'}, 'electronic minimiser'); done('ALGO = All, NELM = 200')
        set_incar(d, {'ALGO': 'Normal', 'NELM': '200'}, 'electronic minimiser'); done('ALGO = Normal, NELM = 200')
    if 'mixing' in msg:
        restart_from_contcar(d); set_incar(d, {'AMIX': '0.1', 'BMIX': '0.01'}, 'charge mixing'); done('AMIX = 0.1, BMIX = 0.01')
    if 'symmetry' in msg:
        set_incar(d, {'ISYM': '0'}, 'k-point/symmetry'); done('ISYM = 0')
    log.write('attempt %d: %s | no automatic fix - needs a human\n' % (n, msg)); log.close()
    say('%s: %s (no automatic fix; see FIXES.txt)' % (os.path.relpath(d, ROOT), msg), 'err'); sys.exit(3)

# ------------------------------------------------------------------ run / status
def job_script(c, name, body, walltime):
    sched = c.get('scheduler', 'type', fallback='local')
    nt, nn = c.get('run', 'ntasks', fallback='32'), c.get('scheduler', 'nodes', fallback='1')
    part, acct = c.get('scheduler', 'partition', fallback=''), c.get('scheduler', 'account', fallback='')
    h = ['#!/bin/bash']
    if sched == 'slurm':
        h += ['#SBATCH -J %s' % name, '#SBATCH -N %s' % nn, '#SBATCH -n %s' % nt, '#SBATCH -t %s' % walltime,
              '#SBATCH -o %s/logs/%s.%%j.out' % (ROOT, name)]
        if part: h.append('#SBATCH -p %s' % part)
        if acct: h.append('#SBATCH -A %s' % acct)
    elif sched == 'pbs':
        ppn = max(1, int(nt) // max(1, int(nn)))
        h += ['#PBS -N %s' % name, '#PBS -l nodes=%s:ppn=%d' % (nn, ppn), '#PBS -l walltime=%s' % walltime,
              '#PBS -j oe', '#PBS -o %s/logs/%s.out' % (ROOT, name)]
        if part: h.append('#PBS -q %s' % part)
        if acct: h.append('#PBS -A %s' % acct)
    h += ['set -o pipefail', 'cd %s' % sq(ROOT), body,
          'python3 %s results >/dev/null 2>&1 || true' % sq(os.path.join(ROOT, 'culita.py'))]
    os.makedirs(os.path.join(ROOT, 'jobs'), exist_ok=True); os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    fn = os.path.join(ROOT, 'jobs', name + '.sh'); open(fn, 'w').write('\n'.join(h) + '\n'); os.chmod(fn, 0o755)
    return fn

def submit(c, fn, local):
    sched = 'local' if local else c.get('scheduler', 'type', fallback='local')
    if sched == 'slurm':
        o = sh('sbatch %s' % sq(fn)); m = re.search(r'(\d+)', o)
        if not m: die('sbatch failed: ' + o.strip(), 'check partition/account/time in culita.cfg')
        return 'slurm:' + m.group(1)
    if sched == 'pbs':
        o = sh('qsub %s' % sq(fn)).strip()
        if not o or 'error' in o.lower(): die('qsub failed: ' + o, 'check queue/account/time in culita.cfg')
        return 'pbs:' + o
    log = os.path.join(ROOT, 'logs', os.path.basename(fn)[:-3] + '.local.out')
    p = subprocess.Popen('nohup bash %s > %s 2>&1' % (sq(fn), sq(log)), shell=True)
    return 'local:%d' % p.pid

def cmd_run(args):
    c = load_cfg()
    if not os.path.isfile(CFG): die('not configured yet', 'run ./culita setup first')
    write_config_sh(c)
    if not preflight(c): die('pre-flight checks failed; nothing was submitted', 'fix the FAIL items above (./culita setup)')
    local = '--local' in args or c.get('scheduler', 'type', fallback='local') == 'local'
    do_bulk = '--core' not in args; do_core = '--bulk' not in args
    state = json.load(open(STATE)) if os.path.isfile(STATE) else {}
    if do_bulk:
        fn = job_script(c, 'culita_bulk', 'bash scripts/run_bulk.sh', c.get('scheduler', 'time_bulk', fallback='06:00:00'))
        state['bulk'] = submit(c, fn, local); say('bulk chemical potentials submitted (%s)' % state['bulk'], 'ok')
    if do_core:
        for X in core_jobs(c):
            fn = job_script(c, 'culita_core_' + X, 'bash scripts/run_core.sh %s' % X, c.get('scheduler', 'time_core', fallback='48:00:00'))
            state['core_' + X] = submit(c, fn, local); say('core substitution %s submitted (%s)' % (X, state['core_' + X]), 'ok')
    json.dump(state, open(STATE, 'w'), indent=1)
    print('\nFollow progress with ./culita status ; when everything is done: ./culita pack')

def cmd_status(args):
    c = load_cfg(); state = json.load(open(STATE)) if os.path.isfile(STATE) else {}
    if state:
        print(C['b'] + 'Jobs' + C['end'])
        for k, v in state.items():
            kind, jid = v.split(':', 1)
            q = ''
            if kind == 'slurm': q = sh('squeue -h -j %s -o %%T' % jid).strip() or 'finished'
            elif kind == 'pbs': q = (re.findall(r'job_state = (\S+)', sh('qstat -f %s' % jid)) or ['finished'])[0]
            elif kind == 'local': q = 'running' if os.path.exists('/proc/%s' % jid) else 'finished'
            print('  %-10s %-22s %s' % (k, v, q))
    print(C['b'] + 'Calculations' + C['end'])
    rows = [('bulk/%s/%s' % (el, st), os.path.join(ROOT, 'bulk', el, st)) for el in BULK for st in ('relax1', 'relax2', 'static', 'static_dense')]
    rows += [('core/%s' % X, os.path.join(ROOT, 'core', X)) for X in core_jobs(c)]
    n = {'done': 0}
    for lab, d in rows:
        st, msg = diagnose(d); n[st] = n.get(st, 0) + 1
        kind = {'done': 'ok', 'failed': 'err', 'incomplete': 'warn'}.get(st)
        if st != 'not started' or '--all' in args:
            fx = os.path.join(d, 'FIXES.txt')
            extra = ('  [%d automatic attempt(s), see FIXES.txt]' % len(open(fx).read().splitlines())) if os.path.isfile(fx) else ''
            say('%-24s %-11s %s%s' % (lab, st, msg, extra), kind)
    print('\n%d of %d steps done' % (n['done'], len(rows)) + ('' if n['done'] == len(rows) else '  (./culita status --all lists the rest)'))

def cmd_results(args):
    c = load_cfg(); r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'collect.py')], env=env_from_cfg(c),
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    print(r.stdout)

def cmd_pack(args):
    c = load_cfg(); cmd_results([])
    name = 'culita_results_%s.tar.gz' % datetime.datetime.now().strftime('%Y%m%d_%H%M')
    with tarfile.open(os.path.join(ROOT, name), 'w:gz') as t:
        for f in ['results.txt', 'culita.cfg']:
            if os.path.isfile(os.path.join(ROOT, f)): t.add(os.path.join(ROOT, f), f)
        for d in glob.glob(os.path.join(ROOT, 'bulk', '*', '*')) + glob.glob(os.path.join(ROOT, 'core', '*')):
            if not os.path.isdir(d): continue
            for f in ('INCAR', 'KPOINTS', 'POSCAR', 'CONTCAR', 'OSZICAR', 'FAILED.txt', 'SOURCE.txt') + (('OUTCAR',) if '/bulk/' in d else ()):
                p = os.path.join(d, f)
                if os.path.isfile(p): t.add(p, os.path.relpath(p, ROOT))
            p = os.path.join(d, 'OUTCAR')
            if '/core/' in d and os.path.isfile(p):   # core OUTCARs can be large: keep the parts that matter
                tmp = tempfile.NamedTemporaryFile('w', delete=False, suffix='.txt')
                L = open(p, errors='replace').read().splitlines(); tmp.write('\n'.join(L[:400] + ['...'] + L[-3000:])); tmp.close()
                t.add(tmp.name, os.path.relpath(p, ROOT) + '.head_tail.txt'); os.unlink(tmp.name)
    say('wrote %s - please send this file back' % name, 'ok')

# ------------------------------------------------------------------ selftest
FAKE_VASP = r'''
import os, sys
L = open('POSCAR').read().splitlines(); sp = L[5].split(); n = [int(x) for x in L[6].split()]; N = sum(n)
e = {'Cu': -3.728, 'Li': -2.051, 'Ta': -11.865, 'Nb': -10.10, 'Zr': -8.55, 'V': -9.08, 'Mo': -10.85, 'W': -12.95, 'Ag': -2.83}
import re
if N < 10:
    E0 = sum(e[s] * k for s, k in zip(sp, n)); TOT = E0
else:
    m = re.search(r'(?:X=|-> )([A-Z][a-z]?)', L[0])
    if m is None:                                   # host
        E0 = -3000.0
    else:
        X = m.group(1); E0 = -3000.0 + {'Ta': 0.81, 'Nb': 0.40, 'Zr': 0.20}[X] + e[X] - e['Cu']
    TOT = E0
inc = open('INCAR').read(); assert 'ICHARG = 1' not in inc
open('CONTCAR', 'w').write(open('POSCAR').read())
open('OSZICAR', 'w').write('DAV:   1  -1.0\nDAV:   2  -1.0\n   1 F= %.5f E0= %.5f\n' % (TOT, E0))
open('OUTCAR', 'w').write('   NIONS = %d\n POSITION   TOTAL-FORCE (eV/Angst)\n ---\n' % N + '\n'.join(['   0 0 0  0.001 0.0 0.0'] * N) +
    '\n  external pressure =  0.20 kB\n  free  energy   TOTEN  =  %.5f eV\n  energy  without entropy=  %.5f  energy(sigma->0) =  %.5f\n reached required accuracy\n General timing and accounting informations for this job:\n' % (TOT, E0, E0))
'''

def cmd_selftest(args):
    print(C['b'] + 'Self-test with a fake VASP (no cluster time used)' + C['end'])
    tmp = tempfile.mkdtemp(prefix='culita_selftest_')
    shutil.copytree(ROOT, os.path.join(tmp, 'repo'), ignore=shutil.ignore_patterns('culita.cfg', 'config.sh', '.culita_state.json',
                    'jobs', 'logs', 'results.txt', '*.tar.gz', 'relax1', 'relax2', 'static', 'static_dense', '.git'))
    R = os.path.join(tmp, 'repo')
    for d in glob.glob(os.path.join(R, 'core', '*')): shutil.rmtree(d)
    lib = os.path.join(tmp, 'lib')
    for vs in VARIANTS.values():
        for v in vs:
            os.makedirs(os.path.join(lib, v), exist_ok=True)
            open(os.path.join(lib, v, 'POTCAR'), 'w').write('   TITEL  = PAW_PBE %s 01Jan2000\n' % v)
    fake = os.path.join(tmp, 'fakevasp.py'); open(fake, 'w').write(FAKE_VASP)
    c = configparser.ConfigParser(); c.optionxform = str
    c.read_dict({'vasp': {'std': '%s %s' % (sys.executable, fake), 'gam': '%s %s' % (sys.executable, fake), 'modules': ''},
                 'paths': {'potcar_lib': lib},
                 'scheduler': {'type': 'local'}, 'run': {'core_elements': 'Nb Zr', 'ntasks': '1'}})
    with open(os.path.join(R, 'culita.cfg'), 'w') as f: c.write(f)
    env = dict(os.environ)
    r = subprocess.run([sys.executable, os.path.join(R, 'culita.py'), 'run', '--local', '--wait'], env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    res = open(os.path.join(R, 'results.txt')).read() if os.path.isfile(os.path.join(R, 'results.txt')) else ''
    checks = [('run finished (bulk; core host, Ta, Nb, Zr)', r.returncode == 0),
              ('bulk controls Cu, Li, Ta PASS', res.count('control PASS') == 3),
              ('E_subs(Ta) = +0.810 eV', 'E_subs(Ta) = +0.810 eV in this cell' in res),
              ('E_subs(Nb) = +0.400 eV', 'E_subs(Nb) = +0.400' in res),
              ('E_subs(Zr) - E_subs(Ta) = -0.610 eV', 'E_subs(Zr) - E_subs(Ta) = -0.610' in res)]
    ok = all(v for _, v in checks)
    for lab, v in checks: say(lab, 'ok' if v else 'err')
    if not ok:
        print(r.stdout[-3000:]); print(res)
    say('SELFTEST %s  (work folder: %s)' % ('PASSED' if ok else 'FAILED', tmp), 'ok' if ok else 'err')
    if ok: shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(0 if ok else 1)

# ------------------------------------------------------------------ main
def main():
    a = sys.argv[1:]
    if not a or a[0] in ('-h', '--help', 'help'):
        print(__doc__); return
    cmd, rest = a[0], a[1:]
    if cmd == 'run' and '--wait' in rest:       # foreground (used by selftest): run the job scripts directly
        c = load_cfg(); write_config_sh(c)
        if not preflight(c): sys.exit(1)
        for name, body in [('bulk', 'bash scripts/run_bulk.sh')] + [('core_' + X, 'bash scripts/run_core.sh ' + X) for X in core_jobs(c)]:
            fn = job_script(c, 'culita_' + name, body, '01:00:00')
            r = subprocess.run(['bash', fn], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            if r.returncode: print(r.stdout[-3000:]); sys.exit(1)
        cmd_results([]); return
    fn = {'setup': cmd_setup, 'check': lambda x: sys.exit(0 if preflight(load_cfg()) else 1), 'run': cmd_run,
          'status': cmd_status, 'results': cmd_results, 'pack': cmd_pack, 'selftest': cmd_selftest,
          '_check': lambda x: sys.exit(cmd_check_step(x) or 0), '_fix': cmd_fix_step}.get(cmd)
    if not fn: die('unknown command "%s"' % cmd, './culita help')
    try:
        fn(rest)
    except KeyboardInterrupt:
        print('\ninterrupted')

if __name__ == '__main__':
    main()
