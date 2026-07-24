#!/usr/bin/env python
"""
Generate IDP starting conformers using IDPConformerGenerator.

Generates N conformers for each IDP target using backbone torsion angles
sampled from PDB fragment statistics (Teixeira et al. 2022, J. Phys. Chem. A,
126:5985), with FASPR sidechain packing.  Selects one conformer per target
whose Rg is closest to the experimental value.

Pipeline:
  1. (Optional) Build a torsion-angle database from PISCES-culled PDB structures.
  2. Generate N conformers via ``idpconfgen build`` with ``--dany --dloop-off``.
  3. Compute radius of gyration (Rg) for every conformer.
  4. Select the conformer closest to the experimental Rg.
  5. Copy the selected PDB to the output directory.

Requires the ``idpcg`` pixi environment (see ``devtools/idpcg/pixi.toml``).

On macOS, multiprocessing uses 'fork' to avoid deadlocks with
IDPConformerGenerator's internal Pool usage.  The set_start_method call
must happen before importing IDPConformerGenerator or numba.

Usage examples
--------------
Generate 100 ab40 conformers, select one closest to experimental Rg::

    python generate_idp_conformers.py --target ab40 --database db.json

Run both targets::

    python generate_idp_conformers.py --all --database db.json

Build the torsion database and generate conformers::

    python generate_idp_conformers.py --all --build-database --database-dir ./idpcg_database
"""

from __future__ import annotations

import multiprocessing
import platform
import sys

if platform.system() == "Darwin":
    multiprocessing.set_start_method("fork")

import argparse
import shutil
import subprocess
from pathlib import Path

import numpy as np
import MDAnalysis as mda


# ---------------------------------------------------------------------------
# Target IDP sequences
# These must match the aa_sequence values in
# proteinbenchmark/benchmark_targets.py.  They are duplicated here so that
# this script can run from the idpcg pixi environment without installing
# the proteinbenchmark package.
# ---------------------------------------------------------------------------

TARGETS: dict[str, dict] = {
    "ab40": {
        "sequence": "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV",
        "experimental_rg": 12.0,
    },
    "asyn": {
        "sequence": (
            "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVH"
            "GVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQL"
            "GKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
        ),
        "experimental_rg": 31.0,
    },
}


# ---------------------------------------------------------------------------
# 1. Database building (optional prerequisite)
# ---------------------------------------------------------------------------


def build_database(database_dir: Path, n_workers: int = -1) -> Path:
    """
    Build an IDPConformerGenerator torsion-angle database from a PISCES-culled
    PDB set.  This is a one-time step whose output can be reused across many
    ``idpconfgen build`` runs.

    The three stages are:
      a) ``idpconfgen pdbdl``   -- download culled PDB structures
      b) ``idpconfgen sscalc``  -- assign secondary structure with DSSP
      c) ``idpconfgen torsions`` -- extract backbone/sidechain torsions

    This basically follows the example/small_peptide workflow

    Parameters
    ----------
    database_dir
        Working directory for intermediate files and the final database JSON.
    n_workers
        Number of CPU workers (``-1`` means all-but-one, mirroring the
        ``-n`` flag of IDPConformerGenerator).

    Returns
    -------
    Path to the generated ``idpconfgen_database.json``.
    """
    database_dir.mkdir(parents=True, exist_ok=True)
    db_json = database_dir / "idpconfgen_database.json"
    if db_json.exists():
        print(f"[build_database] Database already exists at {db_json}, skipping.")
        return db_json

    worker_flag = ["-n"] if n_workers == -1 else ["-n", str(n_workers)]

    # Stage a: download PDB structures from a PISCES culled list.
    # The culled list (Dunbrack lab, 90% seq identity, ≤2.0 Å, R≤0.25)
    # is the same one used in the IDPConformerGenerator paper (Teixeira 2022).
    # Expected at database_dir.parent (i.e. alongside the database directory).
    cull_list = database_dir.parent / "cullpdb_pc90_res2.0_R0.25_d201015_chains24003"
    if not cull_list.exists():
        sys.exit(
            f"[build_database] Culled PDB list not found: {cull_list}\n"
            "  Copy it from IDPConformerGenerator/example/drksh3_example/cull.tar"
        )

    pdbs_tar = database_dir / "pdbs.tar"
    if not pdbs_tar.exists():
        print("[build_database] Stage 1/3: downloading culled PDB structures ...")
        _run_idpconfgen(
            ["pdbdl", str(cull_list), "-u", "-d", str(pdbs_tar)] + worker_flag,
            cwd=database_dir,
        )

    # Stage b: compute secondary-structure assignments with DSSP.
    # sscalc writes ``sscalc.json`` and ``sscalc_splitted.tar`` into cwd
    # by default (see IDPConformerGenerator cli_sscalc.py).
    sscalc_json = database_dir / "sscalc.json"
    if not sscalc_json.exists():
        print("[build_database] Stage 2/3: computing secondary structure (DSSP) ...")
        _run_idpconfgen(
            ["sscalc", str(pdbs_tar), "-rd"] + worker_flag,
            cwd=database_dir,
        )

    # Stage c: extract torsion angles from the DSSP-annotated structures.
    sscalc_tar = database_dir / "sscalc_splitted.tar"
    print("[build_database] Stage 3/3: extracting torsion angles ...")
    _run_idpconfgen(
        [
            "torsions",
            str(sscalc_tar),
            "-sc",
            str(sscalc_json),
            "-o",
            str(db_json),
        ]
        + worker_flag,
        cwd=database_dir,
    )

    if not db_json.exists():
        sys.exit(
            f"[build_database] ERROR: expected database at {db_json} but not found."
        )
    print(f"[build_database] Database built successfully: {db_json}")
    return db_json


# ---------------------------------------------------------------------------
# 2. Conformer generation
# ---------------------------------------------------------------------------


def generate_conformers(
    target: str,
    sequence: str,
    n_conformers: int,
    database_path: Path,
    work_dir: Path,
    random_seed: int = 0,
) -> list[Path]:
    """
    Call ``idpconfgen build`` to sample backbone conformers for *sequence*.

    Parameters
    ----------
    target
        Short name used for directory naming (e.g. ``"ab40"``).
    sequence
        One-letter amino-acid sequence.
    n_conformers
        Number of conformers to request.
    database_path
        Path to the pre-built torsion-angle database JSON.
    work_dir
        Parent working directory; conformers are written to
        ``work_dir/<target>_conformers/``.

    Returns
    -------
    Sorted list of PDB file paths that were actually produced.
    """
    if not database_path.exists():
        sys.exit(
            f"[generate] Database not found: {database_path}\n"
            "  Run with --build-database first, or supply a valid --database path."
        )

    conformer_dir = work_dir / f"{target}_conformers"
    conformer_dir.mkdir(parents=True, exist_ok=True)

    # Check whether conformers already exist (allow resuming).
    existing = sorted(conformer_dir.glob("*.pdb"))
    if len(existing) >= n_conformers:
        print(
            f"[generate] {len(existing)} conformers already present in "
            f"{conformer_dir}, skipping generation."
        )
        return existing[:n_conformers]

    print(
        f"[generate] Generating {n_conformers} conformers for {target} "
        f"({len(sequence)} residues) ..."
    )

    # Call idpconfgen.cli_build.main() directly instead of via subprocess so
    # that the fork start method (set at module level) is inherited.
    # --dany samples from all DSSP classes (Teixeira 2022, §3.2).
    # --dloop-off is required when --dany is active (mutual exclusion group).
    from idpconfgen.cli_build import main as idpcg_build

    idpcg_build(
        input_seq=sequence,
        database=str(database_path),
        custom_sampling=None,
        dany=True,
        dloop_off=True,
        nconfs=n_conformers,
        ncores=1,
        output_folder=str(conformer_dir),
        # Must explicitly pass these — conformer_generator() defaults differ
        # from the CLI defaults when main() is called directly.
        disable_sidechains=False,
        sidechain_method="faspr",
        energy_threshold_backbone=100,
        energy_threshold_sidechains=250,
        random_seed=random_seed,
    )

    produced = sorted(conformer_dir.glob("*.pdb"))
    if not produced:
        sys.exit(f"[generate] ERROR: no PDB files found in {conformer_dir}.")
    if len(produced) < n_conformers:
        print(
            f"[generate] WARNING: requested {n_conformers} conformers but only "
            f"{len(produced)} were produced."
        )
    return produced


# ---------------------------------------------------------------------------
# 3. Structural analysis
# ---------------------------------------------------------------------------


def compute_rg(pdb_path: Path) -> float:
    """
    Compute the radius of gyration (Rg) of a single PDB structure.

    Uses MDAnalysis ``AtomGroup.radius_of_gyration()`` which returns Rg in
    angstroms, mass-weighted by default.

    Parameters
    ----------
    pdb_path
        Path to a PDB file.

    Returns
    -------
    Rg in angstroms.
    """
    u = mda.Universe(str(pdb_path))
    return u.atoms.radius_of_gyration()


# ---------------------------------------------------------------------------
# 4. Helpers
# ---------------------------------------------------------------------------


def _run_idpconfgen(args: list[str], cwd: Path | None = None) -> None:
    """
    Run an ``idpconfgen`` CLI subcommand, streaming output to stdout.

    Exits on non-zero return code.
    """
    cmd = ["idpconfgen"] + args
    print(f"[cmd] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        sys.exit(f"[cmd] idpconfgen exited with code {result.returncode}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate and select diverse IDP conformers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Target selection.
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument(
        "--target",
        choices=list(TARGETS.keys()),
        help="Run a single target (ab40 or asyn).",
    )
    target_group.add_argument(
        "--all",
        action="store_true",
        dest="run_all",
        help="Run all targets.",
    )

    # Generation parameters.
    parser.add_argument(
        "--n-conformers",
        type=int,
        default=100,
        help="Number of conformers to generate per target (default: 100).",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=0,
        help="Random seed for IDPConformerGenerator (default: 0, deterministic).",
    )

    # Database.
    parser.add_argument(
        "--database",
        type=Path,
        default=None,
        help="Path to a pre-built IDPConformerGenerator torsion-angle database JSON.",
    )
    parser.add_argument(
        "--build-database",
        action="store_true",
        help=(
            "Build the torsion-angle database from scratch before generating "
            "conformers.  Requires DSSP to be installed.  The database is "
            "written to --database-dir."
        ),
    )
    parser.add_argument(
        "--database-dir",
        type=Path,
        default=Path("idpcg_database"),
        help=(
            "Working directory for database building " "(default: ./idpcg_database)."
        ),
    )

    # Output.
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "proteinbenchmark"
        / "data"
        / "pdbs",
        help=(
            "Directory for final selected PDB files "
            "(default: proteinbenchmark/data/pdbs/)."
        ),
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("idpcg_workdir"),
        help="Working directory for intermediate conformer files (default: ./idpcg_workdir).",
    )

    args = parser.parse_args()

    if args.n_conformers < 1:
        parser.error("--n-conformers must be at least 1.")

    # Resolve the database path.
    if args.build_database:
        db_path = build_database(args.database_dir.resolve())
    elif args.database is not None:
        db_path = args.database.resolve()
    else:
        parser.error("Provide --database or use --build-database.")

    # Determine which targets to process.
    if args.run_all:
        targets_to_run = list(TARGETS.keys())
    else:
        targets_to_run = [args.target]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    for target in targets_to_run:
        target_info = TARGETS[target]
        sequence = target_info["sequence"]
        exp_rg = target_info["experimental_rg"]
        print(f"\n{'='*60}")
        print(f"Processing target: {target} ({len(sequence)} residues)")
        print(f"Experimental Rg: {exp_rg:.1f} angstrom")
        print(f"{'='*60}")

        # --- Step 1: generate conformers ---
        pdb_paths = generate_conformers(
            target=target,
            sequence=sequence,
            n_conformers=args.n_conformers,
            database_path=db_path,
            work_dir=args.work_dir.resolve(),
            random_seed=args.random_seed,
        )

        # --- Step 2: compute Rg for every conformer ---
        print(f"[analysis] Computing Rg for {len(pdb_paths)} conformers ...")
        rg_values = np.array([compute_rg(p) for p in pdb_paths])
        print(
            f"[analysis] Rg range: {rg_values.min():.2f} -- "
            f"{rg_values.max():.2f} angstrom"
        )

        # --- Step 3: select the conformer closest to experimental Rg ---
        best_idx = int(np.argmin(np.abs(rg_values - exp_rg)))
        print(
            f"[select] Conformer {pdb_paths[best_idx].name}: "
            f"Rg = {rg_values[best_idx]:.2f} angstrom "
            f"(target {exp_rg:.1f}, delta = {abs(rg_values[best_idx] - exp_rg):.2f})"
        )

        # --- Step 4: copy selected PDB to output directory ---
        output_dir = args.output_dir.resolve()
        dest = output_dir / f"{target}.pdb"
        shutil.copy2(pdb_paths[best_idx], dest)
        print(f"[output] {pdb_paths[best_idx].name} -> {dest}")



if __name__ == "__main__":
    main()
