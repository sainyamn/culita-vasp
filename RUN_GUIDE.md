# What to run for the two requested calculations

Requested by K. Solanki (9 Oct 2026), with the same VASP settings as the paper's interface calculations:

| Request | What it gives the paper | Folder in this repo | Jobs |
|---|---|---|---|
| 1. Bulk chemical potentials μ of **Zr, Nb, V, Mo, W, Ag** | replaces the Materials Project μ in Table S1 (γint of the element screen) | `bulk/Zr`, `bulk/Nb`, `bulk/V`, `bulk/Mo`, `bulk/W`, `bulk/Ag` | 1 job, all elements in sequence |
| (controls, run automatically) | must reproduce the paper's μCu = −3.728, μLi = −2.051, μTa = −11.865 eV; confirms the settings and potentials match | `bulk/Cu`, `bulk/Li`, `bulk/Ta` | same job |
| 2. **Nb** (and **Zr**) for one Cu in the core of the embedded precipitate | E_subs(Nb), E_subs(Zr) next to E_subs(Ta) = +0.81 eV (Fig. 2, Table S2) | `core_inputs/Nb`, `core_inputs/Zr` | 1 job each |
| (needed for request 2) | the same cell without the substitution, and Ta on the same site so all three compare like-for-like | `core_inputs/host`, `core_inputs/Ta` | 1 job each |

Nothing else is needed: no earlier calculation folders are used, and `./culita run` submits exactly these jobs.

## Steps on the cluster (about 10 minutes of hands-on time)

```bash
git clone https://github.com/sainyamn/culita-vasp.git
cd culita-vasp
./culita selftest      # 1 min, no VASP: confirms the tool works on this machine
./culita setup         # answer the questions below
./culita run --bulk    # request 1 first: small, finishes in about an hour
./culita run --core    # request 2: host, Ta, Nb, Zr as four separate jobs
./culita status        # any time; shows done / running / failed with the reason
./culita pack          # when everything is done -> culita_results_<date>.tar.gz
```
Send `culita_results_<date>.tar.gz` to S. Nagar (it contains results.txt, the final structures, OSZICARs and the
bulk OUTCARs; core OUTCARs are trimmed to the parts that matter).

### Answers for `./culita setup`

| Question | What to enter |
|---|---|
| command(s) to load VASP | the module that provides VASP 5.4.4, e.g. `module load vasp/5.4.4` (add `; ulimit -s unlimited` if your group usually does) |
| MPI ranks per job / vasp_std / vasp_gam commands | as your group normally launches VASP; the 459-atom core cells are Γ-only, so `vasp_gam` is used for them if available |
| PBE PAW library | the folder that holds `Cu/POTCAR`, `Li/POTCAR`, `Ta/POTCAR`, `Zr_sv/POTCAR`, `Nb_pv/POTCAR`, `V_pv/POTCAR`, `Mo_pv/POTCAR`, `W_pv/POTCAR`, `Ag/POTCAR` (these variants are fixed in the tool) |
| core substitutions | `Nb Zr` (or just `Nb` if computer time is short) |
| scheduler, partition, account, nodes, wall times | your usual values; suggested wall time 6 h for the bulk job and 48 h per core job (jobs resume automatically if re-run) |

## Cost (rough)

* Bulk: 9 elements × 4 small runs (1–2 atoms, dense k-mesh); typically under an hour each element on one node.
* Core: four relaxations of the 459-atom cell at the Γ point, each comparable to one interface relaxation in the paper
  (typically 1–3 days on 2–4 nodes). The four jobs are independent and can run at the same time.

## What comes back (results.txt)

* μ of every element next to the paper's values, with PASS/CHECK for k-point convergence, cell relaxation and the
  Cu/Li/Ta controls.
* E_subs(Ta), E_subs(Nb), E_subs(Zr) in the same cell, and E_subs(Nb) − E_subs(Ta), E_subs(Zr) − E_subs(Ta), with the
  relaxation convergence and the largest remaining force of each run.
* If the recomputed E_subs(Ta) differs from the published +0.81 eV by more than 0.1 eV, results.txt says so; the
  differences X − Ta are unaffected because Ta, Nb and Zr sit on the same site of the same cell.
