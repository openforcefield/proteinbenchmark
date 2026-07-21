# devtools/idpcg/

Pixi environment for generating IDP starting conformers with
[IDPConformerGenerator](https://github.com/julie-forman-kay-lab/IDPConformerGenerator)
(Teixeira et al. 2022, J. Phys. Chem. A, 126:5985).

## Contents

- `pixi.toml` — environment spec (python <3.12, numpy, pybind11, dssp, mdanalysis)
  with IDPConformerGenerator pinned to commit `e25a7d6`
- `cullpdb_pc90_res2.0_R0.25_d201015_chains24003` — PISCES culled PDB list
  (Dunbrack lab, 90% seq identity, ≤2.0 Å, R≤0.25, Oct 2020). Same list used
  in the IDPConformerGenerator paper; extracted from
  `IDPConformerGenerator/example/drksh3_example/cull.tar`
- `results/` — diagnostic plots (Rg scatter, RMSD heatmaps) from generation runs

## Gitignored (large, regenerable)

- `.pixi/` — installed environment (~500 MB)
- `idpcg_database/` — torsion-angle database built from the PISCES list:
  - `pdbs.tar` — downloaded PDB structures (~3.7 GB)
  - `idpconfgen_database.json` — extracted torsion angles (~318 MB)
- `idpcg_workdir/` — intermediate conformer PDBs from generation runs

## How the database was built

```bash
cd devtools/idpcg
pixi install && pixi run install-idpcg
```

All three stages were run on macOS (Apple Silicon) using `forkserver`
multiprocessing with 8 cores, via direct Python API calls:

1. `idpconfgen pdbdl` — downloaded 24,002 chains from RCSB (8 failed/obsolete)
2. `idpconfgen sscalc` — assigned secondary structure with DSSP (mkdssp v3)
3. `idpconfgen torsions` — extracted backbone/sidechain torsion angles

The generation script (`devtools/scripts/generate_idp_conformers.py`) uses this
database to sample conformers. Selected PDBs are written to
`proteinbenchmark/data/pdbs/`.

## macOS note

IDPConformerGenerator's multiprocessing deadlocks on macOS with both `spawn`
(numba JIT pickling) and `fork` (stale locks on successive Pool creation).
The `forkserver` start method avoids both failure modes. This is set at the
top of `generate_idp_conformers.py` before any imports.
