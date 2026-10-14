# culita-vasp

Runs the two remaining VASP calculations for the Cu–Ta–Li bilayer-complexion letter (S. Nagar, K. Solanki), with the
same VASP settings as the paper's interface calculations, and returns the numbers with convergence checks:

1. **Bulk chemical potentials** μ of Zr, Nb, V, Mo, W and Ag (replacing the Materials Project values in the element screen),
   plus Cu, Li and Ta as **controls** that must reproduce the paper's μCu = −3.728, μLi = −2.051, μTa = −11.865 eV.
2. **Nb (and Zr) substituted for one Cu atom in the core of the embedded Cu₃Li precipitate**,
   E_subs = E(core X) − E(host) − μ_X + μ_Cu, compared with Ta (+0.81 eV in the paper).
   `core_inputs/` contains ready POSCARs of the 459-atom Cu | Ta-bilayer | Cu₃Li cell of Fig. S1 (relaxed DFT geometry):
   the host, and the host with Cu atom 196 (a full L1₂ Cu site: 8 Cu + 4 Li neighbours, 6.8 Å from the nearest Ta,
   6.2 Å from the surfaces) replaced by **Ta, Nb or Zr**. All four are relaxed here with the same INCAR, so the set is
   self-contained and needs no earlier calculations; Ta is recomputed on the same site, so E_subs(Nb) − E_subs(Ta) is a
   like-for-like comparison.

No VASP binaries or POTCAR files are included; the tool uses the licensed installation on your cluster.
Requirements: Python ≥ 3.6 (standard library only), VASP 5.4.4, a PBE PAW library, Slurm, PBS or an interactive node.

> **For the two calculations requested for the letter, see [RUN_GUIDE.md](RUN_GUIDE.md)** (what to run and the answers for setup).

## Quick start (about 10 minutes of your time)

```bash
git clone https://github.com/sainyamn/culita-vasp.git
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
| core substitutions to run | `Nb Zr` |
| scheduler, partition, account, nodes, wall times | `slurm`, `general`, `grp_…`, `2`, `48:00:00` |

Everything is saved in `culita.cfg` (plain text; edit and run `./culita check` to re-validate).

## Safety checks before anything is submitted (`./culita check`)

* VASP command and modules resolve; the PAW library has every needed potential.
* The 459-atom POSCARs for host, Ta and each requested element are present in `core_inputs/`.

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
| core reference | E_subs(Ta) recomputed in the same cell and compared with the published +0.81 eV (flagged if it differs by > 0.1 eV) |
| core relaxation | 'reached required accuracy' (EDIFFG = −0.03 eV/Å), and the largest remaining force |
| cutoff | fixed at 500 eV as in the paper (not varied) |

`results.txt` also reports E_subs(X) − E_subs(Ta), which depends only on the two core runs and the bulk μ of X and Ta.
All energies are `energy(sigma->0)`.

## Calculation details

* Bulk: Zr hcp; Nb, Mo, W, V, Ta, Li bcc; Ag, Cu fcc (primitive cells), Γ-centred meshes of ~60 Å density
  (22×22×12 to 29×29×29). ENCUT = 500 eV, PREC = Accurate, PBE, non-spin-polarised, EDIFF = 1E-6;
  two ISIF = 3 relaxations (ISMEAR = 1, SIGMA = 0.1, EDIFFG = −0.005) and a static run with ISMEAR = −5.
* Core runs: `core_inputs/INCAR` = the paper's settings (PBE/PAW, ENCUT = 500 eV, EDIFF = 1E-3 eV,
  Γ point, quasi-Newton IBRION = 1, ISIF = 2) plus EDIFFG = −0.03 eV/Å, ISMEAR = 1 / SIGMA = 0.1, PREC = Normal,
  LREAL = Auto; standard Cu, Li, Ta PAW potentials. host, Ta, Nb and Zr run as four separate jobs.
* PAW variants (fixed): Cu, Li, Ta standard (as stated in the supplement: 11, 1, 5 valence electrons); Zr_sv, Nb_pv,
  V_pv, Mo_pv, W_pv, Ag for the element screen (the usual Materials Project choices).

## Files

| Path | Purpose |
|---|---|
| `culita`, `culita.py` | the command-line tool |
| `bulk/<El>/` | POSCAR, KPOINTS, INCAR.relax, INCAR.static for each element |
| `core_inputs/` | POSCARs of the 459-atom cell: `host/`, `Ta/`, `Nb/`, `Zr/` (one Cu → X at the same core site), shared INCAR and KPOINTS |
| `core/` | empty in the repository: the runs are set up and executed here (copied from `core_inputs/`) |
| `scripts/` | run_bulk.sh, run_core.sh, vasp_step.sh (run → check → fix → retry), make_potcar.py, collect.py |
| `jobs/`, `logs/` | generated job scripts and scheduler logs |
| `results.txt` | the numbers and the convergence checks |
