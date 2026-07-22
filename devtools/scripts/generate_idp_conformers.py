#!/usr/bin/env python
"""
Generate diverse IDP conformer ensembles using IDPConformerGenerator.

Pipeline:
  1. (Optional) Build a torsion-angle database from culled PDB structures.
  2. Generate N conformers for a target IDP sequence via ``idpconfgen build``.
  3. Compute radius of gyration (Rg) and pairwise CA-RMSD for every conformer.
  4. Select K diverse representatives via greedy farthest-point sampling.
  5. Copy selected PDBs to the output directory and write a summary.

Requires the ``idpcg`` pixi environment (see ``devtools/idpcg/pixi.toml``).

On macOS, multiprocessing must use 'forkserver' (not the default 'spawn') to
avoid deadlocks.  'spawn' deadlocks on numba JIT pickling; 'fork' deadlocks
on stale locks when creating successive Pools.  'forkserver' starts a clean
server process before any threads exist, avoiding both failure modes.
The set_start_method call below must happen before any other import.

Usage examples
--------------
Generate 100 ab40 conformers, select 5 diverse ones::

    python generate_idp_conformers.py --target ab40 --database db.json

Run both targets::

    python generate_idp_conformers.py --all --database db.json

Build a database first (requires DSSP)::

    python generate_idp_conformers.py --build-database --database-dir ./db_workdir
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd as mda_rmsd


# ---------------------------------------------------------------------------
# Target IDP sequences
# These must match the aa_sequence values in
# proteinbenchmark/benchmark_targets.py.  They are duplicated here so that
# this script can run from the idpcg pixi environment without installing
# the proteinbenchmark package.
# ---------------------------------------------------------------------------

TARGETS: dict[str, str] = {
    "ab40": "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV",
    "asyn": (
        "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVH"
        "GVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQL"
        "GKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
    ),
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
    # It ships with the drksh3 example and is copied into devtools/idpcg/.
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
    # that the forkserver start method (set at module level) is inherited.
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


def compute_pairwise_rmsd(pdb_paths: list[Path]) -> np.ndarray:
    """
    Compute the symmetric pairwise CA-RMSD matrix for a set of conformers.

    Each pair is optimally superimposed before computing RMSD (Kabsch
    alignment via MDAnalysis ``rms.rmsd`` with ``superposition=True``).

    Parameters
    ----------
    pdb_paths
        List of PDB file paths.

    Returns
    -------
    (N, N) symmetric distance matrix in angstroms.
    """
    n = len(pdb_paths)
    if n == 0:
        return np.empty((0, 0), dtype=np.float64)

    # Pre-extract CA coordinates from every conformer.
    ca_coords: list[np.ndarray] = []
    for path in pdb_paths:
        u = mda.Universe(str(path))
        ca = u.select_atoms("name CA")
        if ca.n_atoms == 0:
            sys.exit(f"[rmsd] No CA atoms found in {path.name}.")
        # Copy coordinates so the Universe can be garbage-collected.
        ca_coords.append(ca.positions.copy())

    # Verify all conformers have the same number of CA atoms.
    n_ca = ca_coords[0].shape[0]
    for i, coords in enumerate(ca_coords):
        if coords.shape[0] != n_ca:
            sys.exit(
                f"[rmsd] Conformer {pdb_paths[i].name} has {coords.shape[0]} "
                f"CA atoms, expected {n_ca}."
            )

    # Build the symmetric RMSD matrix.
    rmsd_matrix = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            # mda_rmsd returns the RMSD in angstroms after optimal
            # superposition (center=True, superposition=True).
            r = mda_rmsd(
                ca_coords[i],
                ca_coords[j],
                center=True,
                superposition=True,
            )
            rmsd_matrix[i, j] = r
            rmsd_matrix[j, i] = r

    return rmsd_matrix


# ---------------------------------------------------------------------------
# 4. Diversity selection
# ---------------------------------------------------------------------------


def select_diverse_conformers(
    pdb_paths: list[Path],
    n_select: int,
    rmsd_matrix: np.ndarray,
) -> list[int]:
    """
    Select *n_select* maximally diverse conformers using a greedy
    farthest-point (max-min) algorithm on the pairwise RMSD matrix.

    This is a greedy heuristic that approximately maximizes the minimum
    pairwise RMSD among the selected set.  The selected conformers are
    diverse initialization structures for MD, not a statistically
    representative IDP ensemble.

    Procedure:
      1. Seed with the conformer that has the largest mean RMSD to all others
         (the most "globally outlying" structure).
      2. Iteratively add the conformer whose minimum distance to the
         already-selected set is largest.

    If fewer conformers are available than requested, all indices are returned.

    Parameters
    ----------
    pdb_paths
        List of conformer PDB paths (only used for length).
    n_select
        Number of representatives to pick.
    rmsd_matrix
        (N, N) symmetric RMSD distance matrix.

    Returns
    -------
    List of integer indices into *pdb_paths* for the selected conformers.
    """
    n = len(pdb_paths)
    if n <= n_select:
        print(
            f"[select] Only {n} conformers available, returning all "
            f"(requested {n_select})."
        )
        return list(range(n))

    # Seed: pick the conformer with the highest mean RMSD to all others.
    selected: list[int] = [int(np.argmax(rmsd_matrix.mean(axis=1)))]

    for _ in range(n_select - 1):
        # For each candidate, find its minimum distance to the selected set.
        min_dists = rmsd_matrix[:, selected].min(axis=1)
        # Exclude already-selected conformers from consideration.
        min_dists[selected] = -1.0
        # Pick the candidate farthest from its nearest selected neighbor.
        selected.append(int(np.argmax(min_dists)))

    return selected


# ---------------------------------------------------------------------------
# 5. Plotting
# ---------------------------------------------------------------------------


def plot_summary(
    rg_values: np.ndarray,
    rmsd_matrix: np.ndarray,
    selected_indices: list[int],
    target: str,
    plot_dir: Path,
) -> None:
    """
    Two summary plots:
      1. Rg vs. conformer index with selected representatives highlighted.
      2. Pairwise CA-RMSD heatmap with selected conformers marked.

    Parameters
    ----------
    rg_values
        Array of Rg values (angstroms), one per conformer.
    rmsd_matrix
        (N, N) symmetric pairwise CA-RMSD matrix in angstroms.
    selected_indices
        Indices of the selected diverse conformers.
    target
        Target name (used in titles and filenames).
    plot_dir
        Directory to save PNG figures (300 dpi).
    """
    sel = np.array(selected_indices)

    # --- Plot 1: Rg scatter ---
    fig, ax = plt.subplots(figsize=(8, 4))
    indices = np.arange(len(rg_values))

    ax.scatter(
        indices,
        rg_values,
        s=20,
        color="0.6",
        edgecolors="none",
        label="all conformers",
    )
    ax.scatter(
        sel,
        rg_values[sel],
        s=80,
        color="tab:red",
        edgecolors="black",
        linewidths=0.5,
        zorder=5,
        label="selected representatives",
    )
    ax.set_xlabel("Conformer index")
    ax.set_ylabel("Radius of gyration (angstrom)")
    ax.set_title(f"{target}: Rg across generated conformers")
    ax.legend(frameon=False)
    fig.tight_layout()

    rg_path = plot_dir / f"{target}_rg_summary.png"
    fig.savefig(str(rg_path), dpi=300)
    plt.close(fig)
    print(f"[plot] Saved Rg plot to {rg_path}")

    # --- Plot 2: Pairwise RMSD heatmap of the 5 selected conformers ---
    sel_rmsd = rmsd_matrix[np.ix_(sel, sel)]
    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(sel_rmsd, cmap="viridis", origin="lower")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("CA-RMSD (angstrom)")

    # Label ticks as conf1..conf5
    tick_labels = [f"conf{i+1}" for i in range(len(sel))]
    ax.set_xticks(range(len(sel)))
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels(tick_labels)

    # Annotate cells with RMSD values
    for i in range(len(sel)):
        for j in range(len(sel)):
            ax.text(j, i, f"{sel_rmsd[i, j]:.1f}", ha="center", va="center",
                    fontsize=8, color="white" if sel_rmsd[i, j] > sel_rmsd.max() * 0.5 else "black")

    ax.set_title(f"{target}: pairwise CA-RMSD of selected conformers")
    fig.tight_layout()

    rmsd_path = plot_dir / f"{target}_rmsd_matrix.png"
    fig.savefig(str(rmsd_path), dpi=300)
    plt.close(fig)
    print(f"[plot] Saved RMSD matrix to {rmsd_path}")


# ---------------------------------------------------------------------------
# 6. Summary table
# ---------------------------------------------------------------------------


def print_summary_table(
    pdb_paths: list[Path],
    rg_values: np.ndarray,
    rmsd_matrix: np.ndarray,
    selected_indices: list[int],
) -> None:
    """
    Print a formatted table of selected conformers with Rg and mean RMSD to
    all other conformers.
    """
    # Mean RMSD of each conformer to every other conformer in the full set.
    n = rmsd_matrix.shape[0]
    mean_rmsd_all = rmsd_matrix.sum(axis=1) / (n - 1) if n > 1 else np.zeros(n)

    header = f"{'Conf':>12s}  {'Rg (A)':>8s}  {'Mean RMSD (A)':>14s}  {'File':>s}"
    print()
    print(header)
    print("-" * len(header))
    for rank, idx in enumerate(selected_indices, start=1):
        print(
            f"  {rank:>8d}    {rg_values[idx]:8.2f}  {mean_rmsd_all[idx]:14.2f}  "
            f"{pdb_paths[idx].name}"
        )
    print()


# ---------------------------------------------------------------------------
# Helpers
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
        "--n-select",
        type=int,
        default=5,
        help="Number of diverse representatives to select (default: 5).",
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
        "--plot-dir",
        type=Path,
        default=None,
        help="Directory for summary plots (default: same as --output-dir).",
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
    if args.n_select < 1:
        parser.error("--n-select must be at least 1.")

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
        sequence = TARGETS[target]
        print(f"\n{'='*60}")
        print(f"Processing target: {target} ({len(sequence)} residues)")
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

        # --- Step 3: compute pairwise CA-RMSD matrix ---
        print("[analysis] Computing pairwise CA-RMSD matrix ...")
        rmsd_matrix = compute_pairwise_rmsd(pdb_paths)
        triu_vals = rmsd_matrix[np.triu_indices_from(rmsd_matrix, k=1)]
        mean_rmsd = triu_vals.mean() if triu_vals.size > 0 else 0.0
        print(f"[analysis] Mean pairwise RMSD: {mean_rmsd:.2f} angstrom")

        # --- Step 4: select diverse representatives ---
        n_select = min(args.n_select, len(pdb_paths))
        selected = select_diverse_conformers(pdb_paths, n_select, rmsd_matrix)

        # --- Step 5: copy selected PDBs to output directory ---
        output_dir = args.output_dir.resolve()
        for rank, idx in enumerate(selected, start=1):
            dest = output_dir / f"{target}_conf{rank}.pdb"
            shutil.copy2(pdb_paths[idx], dest)
            print(f"[output] {pdb_paths[idx].name} -> {dest}")

        # --- Step 6: print summary table ---
        print_summary_table(pdb_paths, rg_values, rmsd_matrix, selected)

        # --- Step 7: summary plots (Rg scatter + RMSD heatmap) ---
        plot_dir = (args.plot_dir or args.output_dir).resolve()
        plot_dir.mkdir(parents=True, exist_ok=True)
        plot_summary(rg_values, rmsd_matrix, selected, target, plot_dir)


if __name__ == "__main__":
    main()
