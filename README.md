# culita-vasp

Runs the two remaining VASP calculations for the Cu–Ta–Li bilayer-complexion letter (S. Nagar, K. Solanki), with the
same settings as the original interface runs, and returns the numbers with convergence checks:

1. **Bulk chemical potentials** μ of Zr, Nb, V, Mo, W and Ag (replacing the Materials Project values in the element screen),
   plus Cu, Li and Ta as **controls** that must reproduce the paper's μCu = −3.728, μLi = −2.051, μTa = −11.865 eV.
2. **Nb (and Zr) substituted for one Cu atom in the core of the embedded Cu₃Li precipitate**,
   E_subs = E(core X) − E(host) − μ_X + μ_Cu, compared with Ta (+0.81 eV in the paper). Two modes:
   * **packaged (default, needs nothing from the old runs):** `core_inputs/` contains ready POSCARs of the 459-atom
     Cu | Ta-bilayer | Cu₃Li cell of Fig. S1 (relaxed DFT geometry, rebuilt from the Bader ACF.dat of the original Ta run):
     the host, and the host with Cu atom 196 (a full L1₂ Cu site: 8 Cu + 4 Li neighbours, 6.8 Å from the nearest Ta,
     6.2 Å from the surfaces) replaced by **Ta, Nb or Zr**. All four are relaxed with the same INCAR, so the set is
     self-contained; Ta is recomputed on the same site, so E_subs(Nb) − E_subs(Ta) is a like-for-like comparison.
   * **original (optional):** if the original Ta-core run and its host run are still on the cluster, give their folders
     in `setup`; Nb/Zr then replace the core Ta of that exact run, reusing its INCAR, KPOINTS and geometry.

No VASP binaries or POTCAR files are included; the tool uses the licensed installation on your cluster.
Requirements: Python ≥ 3.6 (standard library only), VASP 5.4.4, a PBE PAW library, Slurm, PBS or an interactive node.

> **For the two calculations requested for the letter, see [RUN_GUIDE.md](RUN_GUIDE.md)** (what to run and the answers for setup).

## Quick start (about 10 minutes of your time)

```bash
git clone https://github.com/<user>/culita-vasp.git
cd culita-vasp
./culita selftest     # 1 min, no VASP needed: proves the tool works on this machine
./culita setup        # answers a few questions; detects VASP, the PAW library and the scheduler
./culita run          # pre-flight checks, then submits the jobs (bulk; host, Ta, Nb, Zr core cells)
./culita status       # any time: what is done, running or failed (and why)
./culita pack         # when finished: creates culita_results_<date>.tar.gz - send this file back
```

`setup` asks for:

| Question | Example |
|---|---|
| command(s) that load VASP | `module load vasp/5.4.4` |
| MPI ranks per job and the vasp_std / vasp_gam commands | `32`, `mpirun -np 32 vasp_std` |
| PBE PAW library (folder with `Cu/POTCAR`, `Ta/POTCAR`, …) | `/packages/vasp/potpaw_PBE` |
| folder of the original Ta-in-core run (optional; Enter = packaged cell) | needs INCAR, KPOINTS, POTCAR, CONTCAR, OUTCAR |
| folder of the same cell without the core Ta (only if the above is given) | same files |
| optional: the original 459-atom element-screen runs | so μ uses exactly the same PAW variants (e.g. Nb_pv vs Nb_sv) |
| scheduler, partition, account, nodes, wall times | `slurm`, `general`, `grp_…`, `2`, `48:00:00` |

Everything is saved in `culita.cfg` (plain text; edit and run `./culita check` to re-validate).

## Safety checks before anything is submitted (`./culita check`)

* VASP command and modules resolve; the PAW library has every needed potential.
* The two original runs exist, have final energies, the same number of atoms, and differ by exactly one Ta ↔ Cu.
* The original Ta value is recomputed from them (should be +0.81 eV).

## Automatic error handling while running

Each VASP step is checked when it ends (normal termination, ionic convergence, electronic convergence, known VASP
error messages). Common failures are fixed and retried automatically, up to 4 times; every fix is logged in
`FIXES.txt` in that folder and shown by `./culita status`:

| Problem | Automatic action |
|---|---|
| wall-time stop / relaxation not converged within NSW / ZBRENT | restart from CONTCAR |
| electronic minimiser failure (EDDDAV, ZHEGV, …) or SCF not converged | ALGO = Normal, then ALGO = All; NELM = 200 |
| charge-mixing failure (BRMIX) | AMIX = 0.1, BMIX = 0.01 |
| k-point / symmetry error | ISYM = 0 |
| out of memory, segfault, VASP not found, POTCAR mismatch | stops that step with a plain-language explanation (needs a person) |

Re-running `./culita run` is always safe: finished steps are skipped and interrupted ones resume.

## How convergence is checked (reported in results.txt as PASS / CHECK)

| Check | Criterion |
|---|---|
| bulk k-points | μ from a 1.35× denser mesh differs by < 1 meV/atom |
| bulk cell | residual pressure < 2 kB; the second ISIF = 3 relaxation reached EDIFFG (removes Pulay stress) |
| settings / potentials | Cu, Li, Ta controls within 10 meV/atom of the paper's μ |
| energy field | the original Ta E_subs is recomputed with both `energy(sigma->0)` and TOTEN; the field that gives +0.81 eV is used for Nb/Zr |
| core relaxation | 'reached required accuracy' with the original EDIFFG, and the largest remaining force |
| cutoff | fixed at 500 eV as in the original runs (not varied) |

`results.txt` also reports E_subs(X) − E_subs(Ta), which depends only on the two core runs and the bulk μ of X and Ta.

## Calculation details

* Bulk: Zr hcp; Nb, Mo, W, V, Ta, Li bcc; Ag, Cu fcc (primitive cells), Γ-centred meshes of ~60 Å density
  (22×22×12 to 29×29×29). ENCUT = 500 eV, PREC = Accurate, PBE, non-spin-polarised, EDIFF = 1E-6;
  two ISIF = 3 relaxations (ISMEAR = 1, SIGMA = 0.1, EDIFFG = −0.005) and a static run with ISMEAR = −5.
* Core runs, packaged mode: `core_inputs/INCAR` = the paper's settings (PBE/PAW, ENCUT = 500 eV, EDIFF = 1E-3 eV,
  Γ point, quasi-Newton IBRION = 1, ISIF = 2) plus EDIFFG = −0.03 eV/Å, ISMEAR = 1 / SIGMA = 0.1, PREC = Normal,
  LREAL = Auto; standard Cu, Li, Ta PAW potentials. host, Ta, Nb and Zr run as four separate jobs.
* Core runs, original mode: the original INCAR is copied unchanged except ISTART = 0 / ICHARG = 2 (no WAVECAR/CHGCAR is copied) and,
  if present, MAGMOM re-mapped to the new atom order. The core Ta is the Ta atom with the fewest Ta neighbours within
  3.3 Å (bilayer Ta atoms have several); its neighbours are printed, and an ambiguous choice stops the tool
  (choose with `python3 scripts/prep_core.py X core/X --index N`).
* PAW variants: read from the original runs where available; otherwise Cu, Li, Ta standard (as stated in the
  supplement) and Zr_sv, Nb_pv, V_pv, Mo_pv, W_pv, Ag, flagged "CHECK" in the log.

## Files

| Path | Purpose |
|---|---|
| `culita`, `culita.py` | the command-line tool |
| `bulk/<El>/` | POSCAR, KPOINTS, INCAR.relax, INCAR.static for each element |
| `core_inputs/` | POSCARs of the 459-atom cell: `host/`, `Ta/`, `Nb/`, `Zr/` (one Cu → X at the same core site), shared INCAR and KPOINTS |
| `core/` | empty in the repository: the runs are set up and executed here (copied from `core_inputs/` or built from the original run) |
| `scripts/` | run_bulk.sh, run_core.sh, vasp_step.sh (run → check → fix → retry), prep_core.py, make_potcar.py, collect.py |
| `jobs/`, `logs/` | generated job scripts and scheduler logs |
| `results.txt` | the numbers and the convergence checks |
