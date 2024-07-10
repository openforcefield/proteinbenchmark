from pathlib import Path
import shutil
import subprocess
from typing import List

import loos
import numpy
import pandas
from loos.pyloos import Trajectory
from openff.toolkit import Molecule
from openff.units import unit
from pymbar import timeseries

from proteinbenchmark.analysis_parameters import *
from proteinbenchmark.benchmark_targets import benchmark_targets, experimental_datasets
from proteinbenchmark.utilities import list_of_dicts_to_csv


def get_timeseries_mean(correlated_timeseries: numpy.ndarray):
    # Get burn-in time, statistical inefficiency, and maximum number of
    # uncorrelated samples from pymbar.timeseries
    t0, g, Neff_max = timeseries.detect_equilibration(correlated_timeseries, nskip=10)

    # Get an uncorrelated sample from the correlated timeseries
    truncated_timeseries = correlated_timeseries[t0:]
    uncorrelated_sample_indices = timeseries.subsample_correlated_data(
        truncated_timeseries, g=g
    )
    uncorrelated_timeseries = truncated_timeseries[uncorrelated_sample_indices]

    # Get mean and standard error of the mean for the correlated, truncated, and
    # uncorrelated timeseries
    N_correlated = correlated_timeseries.size
    correlated_mean = correlated_timeseries.mean()
    correlated_sem = correlated_timeseries.std(ddof=1) / numpy.sqrt(N_correlated)
    N_truncated = truncated_timeseries.size
    truncated_mean = truncated_timeseries.mean()
    truncated_sem = truncated_timeseries.std(ddof=1) / numpy.sqrt(N_truncated)
    N_uncorrelated = uncorrelated_timeseries.size
    uncorrelated_mean = uncorrelated_timeseries.mean()
    uncorrelated_sem = uncorrelated_timeseries.std(ddof=1) / numpy.sqrt(N_uncorrelated)

    return {
        "Statistical Inefficiency": g,
        "Correlated Mean": correlated_mean,
        "Correlated SEM": correlated_sem,
        "Correlated N": N_correlated,
        "Truncated Mean": truncated_mean,
        "Truncated SEM": truncated_sem,
        "Truncated N": N_truncated,
        "Uncorrelated Mean": uncorrelated_mean,
        "Uncorrelated SEM": uncorrelated_sem,
        "Uncorrelated N": N_uncorrelated,
    }


def align_trajectory(
    topology_path: str,
    trajectory_path: str,
    output_prefix: str,
    output_selection: str = 'chainid == "A"',
    align_selection: str = None,
    reference_path: str = None,
    reference_selection: str = None,
):
    """
    Center and align a subset of atoms in a trajectory.

    Parameters
    ---------
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    output_prefix
        The prefix for the path to write the aligned topology and trajectory.
    output_selection
        LOOS selection string for atoms to write to the output. Default is all
        atoms.
    align_selection
        LOOS selection string for atoms to use for alignment. Default is output
        selection.
    reference_path
        The path to the structure used as a reference for alignment. Default is
        the system topology.
    reference_selection
        LOOS selection string for atoms to use for alignment in the reference
        structure. Default is align selection.
    """

    # Load topology and trajectory
    topology = loos.createSystem(topology_path)
    trajectory = Trajectory(trajectory_path, topology)

    # Set up reference structure for alignment
    if reference_path is None:
        reference = topology.copy()
    else:
        reference = loos.createSystem(reference_path)

    # Atom selections
    output_atoms = loos.selectAtoms(topology, output_selection)

    if align_selection is None:
        align_atoms = output_atoms
    else:
        align_atoms = loos.selectAtoms(topology, align_selection)

    if reference_selection is None:
        if align_selection is None:
            reference_atoms = loos.selectAtoms(reference, output_selection)
        else:
            reference_atoms = loos.selectAtoms(reference, align_selection)

    else:
        reference_atoms = loos.selectAtoms(reference, reference_selection)

    # Set up writer for output trajectory
    output_trajectory = loos.DCDWriter(f"{output_prefix}.dcd")

    first_frame = True

    for frame in trajectory:
        # Align frame onto reference
        transform_matrix = align_atoms.superposition(reference_atoms)

        # Apply transformation to output atoms
        output_atoms.applyTransform(loos.loos.XForm(transform_matrix))

        # Recenter output atoms
        output_atoms.centerAtOrigin()

        # Write current frame for output atoms
        output_trajectory.writeFrame(output_atoms)

        # Write a PDB file of the first frame for use as a topology file
        if first_frame:
            first_frame = False
            pdb = loos.PDB.fromAtomicGroup(output_atoms)
            pdb.clearBonds()

            with open(f"{output_prefix}.pdb", "w") as pdb_file:
                pdb_file.write(str(pdb))


def measure_dihedrals(
    topology_path: str,
    trajectory_path: str,
    frame_length: unit.Quantity,
    output_path: str,
):
    """
    Measure the tau angle, backbone dihedrals, and sidechain dihedrals over a
    trajectory.

    Parameters
    ---------
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    frame_length
       The amount of simulation time between frames in the trajectory.
    output_path
        The path to write the time series of dihedrals.
    """

    # Load topology
    topology = loos.createSystem(topology_path)
    min_resid = topology.minResid()
    max_resid = topology.maxResid()
    cterm_resname = topology[len(topology) - 1].resname()

    # Select atoms for dihedrals
    dihedrals_by_residue = list()
    atom_selections = dict()

    for residue in topology.splitByResidue():
        resid = residue[0].resid()
        resname = residue[0].resname()

        if resname in {"NME", "NH2"}:
            continue

        residue_dihedrals = dict()

        for dihedral, dihedral_atom_dict in DIHEDRAL_ATOMS[resname].items():
            if resid == min_resid and dihedral == "phi":
                continue

            if resid == max_resid and dihedral in {"psi", "omega"}:
                continue

            if (
                cterm_resname == "NH2"
                and resid == (max_resid - 1)
                and dihedral == "omega"
            ):
                continue

            # Get atoms in dihedral
            atom_names = dihedral_atom_dict["atom_names"]

            if "resid_offsets" in dihedral_atom_dict:
                atom_resids = [
                    resid + offset for offset in dihedral_atom_dict["resid_offsets"]
                ]

            elif dihedral == "tau":
                atom_resids = [resid] * 3

            else:
                atom_resids = [resid] * 4

            atom_list = list()

            for atom_name, atom_resid in zip(atom_names, atom_resids):
                atom_full_name = f"{atom_name}-{atom_resid}"

                # Atom selection is slow, so only do this once for each atom
                if atom_full_name not in atom_selections:
                    atom_selection = loos.selectAtoms(
                        topology, f'resid == {atom_resid} && name == "{atom_name}"'
                    )

                    if len(atom_selection) == 0:
                        raise ValueError(
                            f"Unable to select atom {atom_name} with resid "
                            f"{atom_resid} for dihedral {dihedral}."
                        )

                    atom_selections[atom_full_name] = atom_selection[0]

                atom_list.append(atom_selections[atom_full_name])

            residue_dihedrals[dihedral] = atom_list

        dihedrals_by_residue.append(
            {
                "resid": resid,
                "resname": resname,
                "dihedrals": residue_dihedrals,
            }
        )

    # Set up trajectory
    trajectory = Trajectory(trajectory_path, topology)
    frame_time = 0.0 * unit.picosecond
    fragment_index = 0
    dihedrals = list()

    # Load one frame into memory at a time
    for frame in trajectory:
        frame_time += frame_length

        frame_index = trajectory.index()
        frame_time_ns = frame_time.m_as(unit.nanosecond)

        # Write dihedrals to file every 10 000 frames to avoid pandas
        # out-of-memory
        if frame_index % 10000 == 9999 and frame_index > 0:
            list_of_dicts_to_csv(dihedrals, f"{output_path}-{fragment_index}")
            fragment_index += 1
            dihedrals = list()

        # Measure dihedrals
        for residue_dict in dihedrals_by_residue:
            for dihedral_name, atoms in residue_dict["dihedrals"].items():
                if len(atoms) == 3:
                    dihedral = loos.angle(atoms[0], atoms[1], atoms[2])

                else:
                    dihedral = loos.torsion(atoms[0], atoms[1], atoms[2], atoms[3])

                dihedrals.append(
                    {
                        "Frame": frame_index,
                        "Time (ns)": frame_time_ns,
                        "Resid": residue_dict["resid"],
                        "Resname": residue_dict["resname"],
                        "Dihedral Name": dihedral_name,
                        "Dihedral (deg)": dihedral,
                    }
                )

    if fragment_index == 0:
        list_of_dicts_to_csv(dihedrals, output_path)

    else:
        list_of_dicts_to_csv(dihedrals, f"{output_path}-{fragment_index}")

    return fragment_index


def measure_h_bond_geometries(
    topology_path: str,
    trajectory_path: str,
    frame_length: unit.Quantity,
    output_path: str,
    h_bond_distance_cutoff: unit.Quantity = 2.5 * unit.angstrom,
    h_bond_angle_cutoff: float = 30,
    occupancy_threshold: float = 0.01,
):
    """
    Measure the donor-acceptor distance, hydrogen-acceptor distance, and donor-
    hydrogen-acceptor angle for hydrogen bonds over a trajectory. For backbone
    N-H---O=C hydrogen bonds, also measure the H-O-C angle and H-O-C-N dihedral.

    Parameters
    ---------
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    frame_length
       The amount of simulation time between frames in the trajectory.
    output_path
        The path to write the hydrogen bond geometries.
    h_bond_distance_cutoff
        Distance between hydrogen donor and acceptor below which a hydrogen bond
        is considered to be occupied.
    h_bond_angle_cutoff
        Deviation in donor-hydrogen-acceptor angle from linear in degrees below
        which a hydrogen bond is considered to be occupied.
    occupancy_threshold
        Fraction of frames in which a putative hydrogen bond must be occupied to
        be considered observed and be measured.
    """

    # Load topology
    topology = loos.createSystem(topology_path)
    min_resid = topology.minResid()
    max_resid = topology.maxResid()

    # Select OFF atoms that can participate in hydrogen bonds
    offtop = Topology.from_pdb(topology_path).molecule(0)
    putative_donors = [
        match.topology_atom_indices
        for match in offtop.chemical_environment_matches("[#7,#8,#16:1]-[#1:2]")
    ]
    putative_acceptors = [
        match.topology_atom_indices
        for match in offtop.chemical_environment_matches("[#7,#8,#16:1]")
    ]

    # Construct a list of [donor, hydrogen, acceptor] LOOS atom selections
    putative_h_bonds = list()
    atom_selections = dict()

    for donor_index, hydrogen_index in putative_donors:
        # Non-hydrogen atoms bonded to donor
        donor_bonded_atoms = [
            atom
            for atom in offtop.atom(donor_index).bonded_atoms
            if atom.atomic_number == 1
        ]

        # Get LOOS atom selections for donor and hydrogen
        if donor_index not in atom_selections:
            donor_resid = offtop.atom(donor_index).metadata["residue_number"]
            donor_name = offtop.atom(donor_index).name

            atom_selection = loos.selectAtoms(
                topology, f'resid == {donor_resid} && name == "{donor_name}"'
            )

            if len(atom_selection) == 0:
                raise ValueError(
                    f"Unable to select donor {donor_name} with resid " f"{donor_resid}."
                )

            atom_selections[donor_index] = atom_selection[0]

        if hydrogen_index not in atom_selections:
            hydrogen_resid = offtop.atom(hydrogen_index).metadata["residue_number"]
            hydrogen_name = offtop.atom(hydrogen_index).name

            atom_selection = loos.selectAtoms(
                topology, f'resid == {hydrogen_resid} && name == "{hydrogen_name}"'
            )

            if len(atom_selection) == 0:
                raise ValueError(
                    f"Unable to select hydrogen {hydrogen_name} with resid "
                    f"{hydrogen_resid}."
                )

            atom_selections[hydrogen_index] = atom_selection[0]

        donor_atom = atom_selections[donor_index]
        hydrogen_atom = atom_selections[hydrogen_index]

        # Find acceptors for this donor
        for (acceptor_index,) in putative_acceptors:
            # Don't consider acceptors within two covalent bonds of the donor
            if acceptor_index == donor_index or any(
                [
                    acceptor_index == bonded_atom.molecule_atom_index
                    or offtop.atom(acceptor_index).is_bonded_to(bonded_atom)
                    for bonded_atom in donor_bonded_atoms
                ]
            ):
                continue

            # Get LOOS atom selection for acceptor
            if acceptor_index not in atom_selections:
                acceptor_resid = offtop.atom(acceptor_index).metadata["residue_number"]
                acceptor_name = offtop.atom(acceptor_index).name

                atom_selection = loos.selectAtoms(
                    topology, f'resid == {acceptor_resid} && name == "{acceptor_name}"'
                )

                if len(atom_selection) == 0:
                    raise ValueError(
                        f"Unable to select acceptor {acceptor_name} with resid "
                        f"{acceptor_resid}."
                    )

                atom_selections[acceptor_index] = atom_selection[0]

            acceptor_atom = atom_selections[acceptor_index]

            putative_h_bonds.append([donor_atom, hydrogen_atom, acceptor_atom])

    # Set up trajectory
    trajectory = Trajectory(trajectory_path, topology)

    # Count observations of putative hydrogen bonds
    h_bond_observations = numpy.zeros(len(putative_h_bonds))

    for frame in trajectory:
        # Set up HBondDetector to determine H bond occupancy from a distance and
        # angle cutoff. Do this for each frame to grab the periodic box vectors.
        h_bond_detector = loos.HBondDetector(
            h_bond_distance_cutoff.m_as(unit.angstrom),
            h_bond_angle_cutoff,
            frame,
        )

        # Check occupancy of remaining putative hydrogen bonds in this frame
        for hb_index, h_bond in enumerate(putative_h_bonds):
            if h_bond_detector.hBonded(h_bond[0], h_bond[1], h_bond[2]):
                h_bond_observations[hb_index] += 1

    # Find hydrogen bonds with occupancy above threshold
    observation_threshold = occupancy_threshold * len(trajectory)
    observed_h_bonds = [
        {
            "donor": h_bond[0],
            "hydrogen": h_bond[1],
            "acceptor": h_bond[2],
        }
        for hb_index, h_bond in enumerate(putative_h_bonds)
        if h_bond_observations[hb_index] >= observation_threshold
    ]

    # Get amide atoms for observed H bonds in backbone
    acceptor_amide_atoms = dict()

    for h_bond in observed_h_bonds:
        donor_resid = h_bond["donor"].resid()
        donor_name = h_bond["donor"].name()
        hydrogen_name = h_bond["hydrogen"].name()
        acceptor_resid = h_bond["acceptor"].resid()
        acceptor_name = h_bond["acceptor"].name()

        # Only consider non-terminal backbone hydrogen bonds
        if (
            donor_name != "N"
            or hydrogen_name != "H"
            or acceptor_name != "O"
            or donor_resid == min_resid
            or acceptor_resid == max_resid
        ):
            continue

        acceptor_index = h_bond["acceptor"].index()

        # Keep track of this in a dictionary because a given acceptor may
        # participate in multiple observed hydrogen bonds
        if acceptor_index not in acceptor_amide_atoms:
            amide_c_selection = loos.selectAtoms(
                topology, f'resid == {acceptor_resid} && name == "C"'
            )

            if len(amide_c_selection) == 0:
                raise ValueError(
                    f"Unable to select amide C with resid {acceptor_resid}."
                )

            amide_n_selection = loos.selectAtoms(
                topology, f'resid == {acceptor_resid + 1} && name == "N"'
            )

            if len(amide_n_selection) == 0:
                raise ValueError(
                    f"Unable to select amide N with resid {acceptor_resid + 1}."
                )

            acceptor_amide_atoms[acceptor_index] = {
                "C": amide_c_selection[0],
                "N": amide_n_selection[0],
            }

        h_bond["acceptor_amide_atoms"] = acceptor_amide_atoms[acceptor_index]

    # Reset Trajectory iterator
    trajectory.reset()

    # Measure hydrogen bond geometries
    frame_time = 0.0 * unit.picosecond
    fragment_index = 0
    h_bond_geometries = list()

    # Load one frame into memory at a time
    for frame in trajectory:
        frame_time += frame_length

        frame_index = trajectory.index()
        frame_time_ns = frame_time.m_as(unit.nanosecond)

        # Write hydrogen bond geomtries to file every 10 000 frames to avoid
        # pandas out-of-memory
        if frame_index % 10000 == 9999 and frame_index > 0:
            list_of_dicts_to_csv(h_bond_geometries, f"{output_path}-{fragment_index}")
            fragment_index += 1
            h_bond_geometries = list()

        # Measure h_bond_geometries
        for h_bond in observed_h_bonds:
            donor = h_bond["donor"]
            hydrogen = h_bond["hydrogen"]
            acceptor = h_bond["acceptor"]

            DA_distance = donor.coords().distance(acceptor.coords())
            HA_distance = hydrogen.coords().distance(acceptor.coords())
            DHA_angle = loos.angle(donor, hydrogen, acceptor)

            h_bond_dict = {
                "Frame": frame_index,
                "Time (ns)": frame_time_ns,
                "Donor Resid": donor.resid(),
                "Donor Resname": donor.resname(),
                "Donor Name": donor.name(),
                "Hydrogen Name": hydrogen.name(),
                "Acceptor Resid": acceptor.resid(),
                "Acceptor Resname": acceptor.resname(),
                "Acceptor Name": acceptor.name(),
                "DA Distance (Angstrom)": DA_distance,
                "HA Distance (Angstrom)": HA_distance,
                "DHA Angle (deg)": DHA_angle,
            }

            if "acceptor_amide_atoms" in h_bond:
                amide_c = h_bond["acceptor_amide_atoms"]["C"]
                amide_n = h_bond["acceptor_amide_atoms"]["N"]

                h_bond_dict["HOC Angle (deg)"] = loos.angle(hydrogen, acceptor, amide_c)
                h_bond_dict["HOCN Dihedral (deg)"] = loos.torsion(
                    hydrogen, acceptor, amide_c, amide_n
                )

            h_bond_geometries.append(h_bond_dict)

    if fragment_index == 0:
        list_of_dicts_to_csv(h_bond_geometries, output_path)

    else:
        list_of_dicts_to_csv(h_bond_geometries, f"{output_path}-{fragment_index}")

    return fragment_index


def compute_chemical_shifts_shiftx2(
    topology_path: str,
    trajectory_path: str,
    frame_length: unit.Quantity,
    output_path: str,
    ph: float,
    temperature: unit.Quantity,
    shiftx2_output_dir: str = None,
    shiftx2_install_dir: str = None,
    python2_path: str = None,
):
    """
    Compute the chemical shifts of protein atoms using ShiftX2.

    Parameters
    ---------
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    frame_length
       The amount of simulation time between frames in the trajectory.
    output_path
        The path to write the computed chemical shifts.
    ph
        The pH of the benchmark target.
    temperature
        The temperature of the benchmark target.
    shiftx2_output_dir
        The directory to write output from ShiftX2.
    shiftx2_install_dir
        The root directory for the installation of ShiftX2.
    python2_path
        The path to the python2 executable to pass to subprocess.
    """

    # ShiftX2 v1.13 has a few quirks that make it unwieldy for analyzing MD
    # trajectories. The ShiftX2 estimate is a combination of the sequence-based
    # ShiftY+ estimator and the structure-based ShiftX+ estimator. ShiftX2 has
    # an NMR mode that will extract models from a PDB file into single-model PDB
    # files, run ShiftX+ on each model PDB, average them, then run ShiftY+ once
    # and combine it with the average ShiftX+ prediction. We want the combined
    # ShiftX+ and ShiftY+ estimates for each frame. Additionally, the main
    # ShiftX2 script calls the shell `python` through `os.system()` but expects
    # python2, making it fail in environments that default to python3. To solve
    # this, we are going to call ShiftY+ once on the topology, extract frames to
    # temporary PDBs, call ShiftX+ on each frame PDB in batch mode, then combine
    # the ShiftX+ and ShiftY+ estimates for each frame. We will use `python2`
    # explicitly in calls to `subprocess.run()`.

    # Get python2 executable
    if python2_path is None:
        python2_path = shutil.which("python2")

    # Get ShiftX2 installation directory
    if shiftx2_install_dir is None:
        shiftx2_install_dir = Path(shutil.which("shiftx2.py")).parent

    # Set up directory to store shiftx2 output
    if shiftx2_output_dir is None:
        shiftx2_output_dir = Path(Path(output_path).parent, "shiftx2")
    else:
        shiftx2_output_dir = Path(shiftx2_output_dir)

    shiftx2_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load topology
    topology = loos.createSystem(topology_path)

    resname_by_resid = {
        residue[0].resid(): residue[0].resname()
        for residue in topology.splitByResidue()
    }

    # Set up trajectory
    trajectory = Trajectory(trajectory_path, topology)

    # Write each frame to temporary PDB
    tmp_pdb_prefix = Path(shiftx2_output_dir, Path(trajectory_path).stem)
    for frame_index, frame in enumerate(trajectory):
        frame_pdb_path = f"{tmp_pdb_prefix}-{frame_index}.pdb"
        with open(frame_pdb_path, "w") as pdb_file:
            pdb_file.write(str(loos.PDB.fromAtomicGroup(frame)))

    # Compute sequence-based estimate of chemical shifts by runnnig ShiftY+ on
    # the topology PDB
    shifty_output_path = Path(
        shiftx2_output_dir,
        f"{Path(topology_path).stem}.shifty",
    )
    shifty_sequence_similarity_cutoff = "40"
    shifty_return_value = subprocess.run(
        [
            python2_path,
            "shifty3.py",
            "-i",
            Path(topology_path).resolve(),
            "-o",
            shifty_output_path.resolve(),
            "-c",
            shifty_sequence_similarity_cutoff,
            "-t",
            shiftx2_output_dir.resolve(),
        ],
        cwd=Path(shiftx2_install_dir, "shifty3"),
    )

    # Compute structure-based estimate of chemical shifts by running ShiftX+
    # on the temporary frame PDBs
    shiftxp_output_path = f"{str(frame_pdb_path)}.sxp"
    shiftxp_bin_dir = Path(shiftx2_install_dir, "bin")
    shiftxp_lib_jar = Path(shiftx2_install_dir, "lib", "weka.jar")
    shiftxp_return_value = subprocess.run(
        [
            "java",
            "-Xmx1900m",
            "-cp",
            f"{shiftxp_bin_dir}:{shiftxp_lib_jar}",
            "ShiftXp",
            "-b",
            f"{tmp_pdb_prefix}-*.pdb",
            "-atoms",
            "ALL",
            "-ph",
            str(ph),
            "-temp",
            str(temperature.m_as(unit.kelvin)),
            "-dir",
            shiftx2_install_dir,
        ]
    )

    # Set up ShiftX2 combine options
    shiftx2_combine_path = Path(shiftx2_install_dir, "script", "combine_cs.py")
    shiftx2_combine_csv_format = "1"
    shiftx2_combine_all_atoms = "1"

    # Set up list of dicts to store chemical shifts
    frame_time = 0.0 * unit.picosecond
    fragment_index = 0
    chemical_shifts = list()

    N_frames = frame_index + 1
    for frame_index in range(N_frames):
        frame_time += frame_length
        frame_time_ns = frame_time.m_as(unit.nanosecond)

        # Write chemical shifts to file every 10 000 frames to avoid pandas
        # out-of-memory
        if frame_index % 10000 == 9999 and frame_index > 0:
            list_of_dicts_to_csv(
                chemical_shifts,
                f"{output_path}-{fragment_index}",
            )
            fragment_index += 1
            chemical_shifts = list()

        # Compute the ShiftX2 estimate of chemical shifts by combining the
        # SHIFTX+ and SHIFTY+ estimates
        shiftxp_output_path = f"{tmp_pdb_prefix}-{frame_index}.pdb.sxp"
        shiftx2_combine_output_path = f"{tmp_pdb_prefix}-{frame_index}.cs"
        shiftx2_combine_return_value = subprocess.run(
            [
                python2_path,
                shiftx2_combine_path,
                shiftxp_output_path,
                shifty_output_path,
                shiftx2_combine_output_path,
                shiftx2_combine_csv_format,
                shiftx2_combine_all_atoms,
            ]
        )

        # Read ShiftX2 estimate
        for index, row in pandas.read_csv(shiftx2_combine_output_path).iterrows():
            resid = row["NUM"]
            chemical_shifts.append(
                {
                    "Frame": frame_index,
                    "Time (ns)": frame_time_ns,
                    "Resid": resid,
                    "Resname": resname_by_resid[resid],
                    "Atom": row["ATOMNAME"],
                    "Chemical Shift": row["SHIFT"],
                }
            )

    if fragment_index == 0:
        list_of_dicts_to_csv(chemical_shifts, output_path)

    else:
        list_of_dicts_to_csv(chemical_shifts, f"{output_path}-{fragment_index}")

    # Remove temporary files
    shutil.rmtree(shiftx2_output_dir)

    return fragment_index


def compute_chemical_shifts_sparta_plus(
    topology_path: str,
    trajectory_path: str,
    frame_length: unit.Quantity,
    output_path: str,
    spartap_output_dir: str = None,
):
    """
    Compute the chemical shifts of protein atoms using SPARTA+.

    Parameters
    ---------
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    frame_length
       The amount of simulation time between frames in the trajectory.
    output_path
        The path to write the computed chemical shifts.
    shiftx2_output_dir
        The directory to write output from SPARTA+.
    """

    # Set up directory to store SPARTA+ output
    if spartap_output_dir is None:
        spartap_output_dir = Path(Path(output_path).parent, "spartap")
    else:
        spartap_output_dir = Path(spartap_output_dir)

    spartap_output_dir.mkdir(parents=True, exist_ok=True)

    # Load topology
    topology = loos.createSystem(topology_path)

    resname_by_resid = {
        residue[0].resid(): residue[0].resname()
        for residue in topology.splitByResidue()
    }

    # Set up trajectory
    trajectory = Trajectory(trajectory_path, topology)

    # Write each frame to temporary PDB
    pdb_filename_prefix = Path(trajectory_path).stem
    tmp_pdb_prefix = Path(spartap_output_dir, pdb_filename_prefix)
    for frame_index, frame in enumerate(trajectory):
        pass
        frame_pdb_path = f"{tmp_pdb_prefix}-{frame_index}.pdb"
        with open(frame_pdb_path, "w") as pdb_file:
            pdb_file.write(str(loos.PDB.fromAtomicGroup(frame)))

        # Estimate chemical shifts by running SPARTA+ on the temporary frame
        # PDBs every 2000 frames to avoid shell error "argument list too long"
        if frame_index % 2000 == 1999 and frame_index > 0:
            spartap_return_value = subprocess.run(
                [
                    "sparta+",
                    "-in",
                    f"{pdb_filename_prefix}-*.pdb",
                ],
                cwd=spartap_output_dir,
            )
            for tmp_file in Path(spartap_output_dir).glob(f"{pdb_filename_prefix}-*.pdb"):
                tmp_file.unlink()

    spartap_return_value = subprocess.run(
        [
            "sparta+",
            "-in",
            f"{pdb_filename_prefix}-*.pdb",
        ],
        cwd=spartap_output_dir,
    )

    # Set up list of dicts to store chemical shifts
    frame_time = 0.0 * unit.picosecond
    fragment_index = 0
    chemical_shifts = list()

    N_frames = frame_index + 1
    for frame_index in range(N_frames):
        frame_time += frame_length
        frame_time_ns = frame_time.m_as(unit.nanosecond)

        # Write chemical shifts to file every 10 000 frames to avoid pandas
        # out-of-memory
        if frame_index % 10000 == 9999 and frame_index > 0:
            list_of_dicts_to_csv(
                chemical_shifts,
                f"{output_path}-{fragment_index}",
            )
            fragment_index += 1
            chemical_shifts = list()

        # Read SPARTA+ estimate
        first_data_line = numpy.infty
        with open(f"{tmp_pdb_prefix}-{frame_index}_pred.tab", "r") as spartap_file:
            for line_index, line in enumerate(spartap_file):
                fields = line.split()
                if line_index >= first_data_line and fields[4] != "9999.000":
                    resid = int(fields[0])
                    chemical_shifts.append(
                        {
                            "Frame": frame_index,
                            "Time (ns)": frame_time_ns,
                            "Resid": resid,
                            "Resname": resname_by_resid[resid],
                            "Atom": fields[2],
                            "Chemical Shift": float(fields[4]),
                            "Secondary Shift": float(fields[3]),
                            "Random Coil Shift": float(fields[5]),
                            "Ring Current Shift": float(fields[6]),
                            "Electric Field Shift": float(fields[7]),
                            "Shift Uncertainty": float(fields[8]),
                        }
                    )
                elif fields and fields[0] == "FORMAT":
                    first_data_line = line_index + 2

    if fragment_index == 0:
        list_of_dicts_to_csv(chemical_shifts, output_path)

    else:
        list_of_dicts_to_csv(chemical_shifts, f"{output_path}-{fragment_index}")

    # Remove temporary files
    shutil.rmtree(spartap_output_dir)

    return fragment_index


def assign_dihedral_clusters(
    dihedrals_path: str,
    output_path: str,
    ramachandran: str = "hollingsworth",
    rotamer: str = "hintze",
):
    """
    Assign frames to Ramachandran clusters based on backbone dihedrals and to
    sidechain rotamers based on sidechain dihedrals.

    Parameters
    ---------
    dihedrals_path
        Path to the time series of dihedrals.
    output_path
        Path to write Ramachandran cluster assignments.
    ramachandran
        The name of the definitions of Ramachandran clusters.
    rotamer
        The name of the rotamer library.
    """

    # Check Ramachandran cluster definition
    ramachandran = ramachandran.lower()

    if ramachandran == "hollingsworth":
        ramachandran_clusters = HOLLINGSWORTH_RAMACHANDRAN_CLUSTERS

    else:
        raise ValueError("Argument `ramachandran` must be one of\n    hollinsworth")

    # Check rotamer library
    rotamer = rotamer.lower()

    if rotamer == "hintze":
        rotamer_library = HINTZE_ROTAMER_LIBRARY

    else:
        raise ValueError("Argument `rotamer` must be one of\n    hintze")

    # Read time series of dihedrals
    dihedral_df = pandas.read_csv(dihedrals_path, index_col=0)

    dihedral_df = (
        dihedral_df.pivot(
            index=["Frame", "Time (ns)", "Resid", "Resname"],
            columns="Dihedral Name",
            values="Dihedral (deg)",
        )
        .add_suffix(" (deg)")
        .reset_index()
        .rename_axis(columns=None)
    )

    # Wrap phi into [0, 360) and psi into [-120, 240) so that no clusters are
    # split across a periodic boundary
    # wrapped_value = (value - lower) % (upper - lower) + lower
    phi = dihedral_df["phi (deg)"] % 360
    psi = (dihedral_df["psi (deg)"] + 120) % 360 - 120

    # Assign residues with defined phi and psi values to the outlier cluster
    dihedral_df["Ramachandran Cluster"] = None
    dihedral_df.loc[~((phi.isna()) | (psi.isna())), "Ramachandran Cluster"] = "Outlier"

    # Assign Ramachandran clusters
    for cluster in ramachandran_clusters:
        dihedral_df.loc[
            (phi > cluster["phi"][0])
            & (phi < cluster["phi"][1])
            & (psi > cluster["psi"][0])
            & (psi < cluster["psi"][1]),
            "Ramachandran Cluster",
        ] = cluster["cluster"]

    # Get list of sidechain dihedrals in target
    target_sidechain_dihedrals = [
        f"chi{i}" for i in range(1, 5) if f"chi{i} (deg)" in dihedral_df.columns
    ]

    # Assign sidechain rotamers
    dihedral_df["Sidechain Rotamer"] = None

    for resname in dihedral_df["Resname"].unique():
        if resname not in rotamer_library:
            continue

        # Assign residues with non-trivial sidechains to the outlier rotamer
        resname_rows = dihedral_df["Resname"] == resname
        dihedral_df.loc[resname_rows, "Sidechain Rotamer"] = "Outlier"

        # Assign rotamers for this residue type
        for rotamer_index, rotamer in rotamer_library[resname].iterrows():
            selected_rows = resname_rows

            for chi in ["chi1", "chi2", "chi3", "chi4"]:
                if f"{chi}_mean" in residue_df.columns:
                    diff = dihedral_df[f"{chi} (deg)"] - rotamer[f"{chi}_mean"]
                    abs_wrapped_diff = ((diff + 180) % 360 - 180).abs()
                    selected_rows = selected_rows & (
                        abs_wrapped_diff < rotamer[f"{chi}_esd"]
                    )

            dihedral_df.loc[selected_rows, "Sidechain Rotamer"] = rotamer["rotamer"]

    dihedral_df.to_csv(output_path)


def compute_scalar_couplings(
    observable_path: str,
    dihedrals_path: str,
    output_path: str,
    karplus: str = "vogeli",
    time_series_output_path: str = None,
):
    """
    Compute NMR scalar couplings using a Karplus model.

    Parameters
    ---------
    observable_path
        The path to the data for experimental observables.
    dihedrals_path
        The path to the time series of dihedrals.
    output_path
        The path to write the computed mean scalar couplings.
    karplus
        The name of the set of Karplus parameters to use.
    time_series_output_path
        The path to write the time series of computed scalar couplings.
    """

    # Get Karplus parameters for scalar couplings associated with phi
    karplus = karplus.lower()

    if karplus == "vogeli":
        karplus_parameters = VOGELI_KARPLUS_PARAMETERS
    elif karplus == "schmidt":
        karplus_parameters = SCHMIDT_KARPLUS_PARAMETERS
    elif karplus == "hu":
        karplus_parameters = HU_KARPLUS_PARAMETERS
    elif karplus == "case_dft1":
        karplus_parameters = CASE_DFT1_KARPLUS_PARAMETERS
    elif karplus == "case_dft2":
        karplus_parameters = CASE_DFT2_KARPLUS_PARAMETERS

    else:
        raise ValueError(
            "Argument `karplus` must be one of\n    vogeli\n    schmidt\n    hu"
            "\n    case_dft1\n    case_dft2"
        )

    # Get Karplus parameters for scalar couplings associated with psi
    karplus_parameters.update(WIRMER_KARPLUS_PARAMETERS)
    karplus_parameters.update(DING_KARPLUS_PARAMETERS)
    karplus_parameters.update(HENNIG_KARPLUS_PARAMETERS)

    # Get Karplus parameters for scalar couplings associated with chi1
    karplus_parameters.update(PEREZ_KARPLUS_PARAMETERS)

    if karplus != "schmidt":
        karplus_parameters.update(CHOU_KARPLUS_PARAMETERS)

    # Load data for experimental observables
    # Take uncertainty from Karplus parameters, not from experiment
    observable_df = pandas.read_csv(
        observable_path,
        sep="\s+",
        usecols=["Observable", "Resid", "Resname", "Experiment"],
    )

    # Read time series of dihedrals
    dihedral_df = pandas.read_csv(dihedrals_path, index_col=0)

    # Skip observables that depend on dihedrals that are undefined
    max_resid = dihedral_df["Resid"].max()
    indices_to_drop = list()
    for index, row in observable_df.iterrows():
        observable = row["Observable"]
        observable_resid = row["Resid"]
        observable_dihedrals = karplus_parameters[observable]["dihedral"].split(",")

        if (observable_resid == 1 and "phi" in observable_dihedrals) or (
            observable_resid == max_resid and "psi" in observable_dihedrals
        ):
            indices_to_drop.append(index)

    observable_df.drop(indices_to_drop, inplace=True)
    observable_df.reset_index(drop=True, inplace=True)

    # Compute observables
    if time_series_output_path is not None:
        observable_timeseries = list()
    computed_observables = list()

    for index, row in observable_df.iterrows():
        observable = row["Observable"]
        observable_resid = row["Resid"]
        observable_parameters = karplus_parameters[observable]

        if observable == "3j_hn_ca":
            # Compute sine and cosine of phi and psi
            # Reset indices for products, e.g. cos phi * cos psi
            karplus_df = dihedral_df[
                (dihedral_df["Dihedral Name"] == "phi")
                & (dihedral_df["Resid"] == observable_resid)
            ]
            phi = numpy.deg2rad(karplus_df["Dihedral (deg)"].values)

            karplus_df = dihedral_df[
                (dihedral_df["Dihedral Name"] == "psi")
                & (dihedral_df["Resid"] == observable_resid - 1)
            ]
            prev_psi = numpy.deg2rad(karplus_df["Dihedral (deg)"].values)

            cos_phi = numpy.cos(phi)
            sin_phi = numpy.sin(phi)
            cos_psi = numpy.cos(prev_psi)
            sin_psi = numpy.sin(prev_psi)

            # Compute estimate for scalar coupling
            computed_coupling = (
                observable_parameters["cos_phi"] * cos_phi
                + observable_parameters["sin_phi"] * sin_phi
                + observable_parameters["cos_psi"] * cos_psi
                + observable_parameters["sin_psi"] * sin_psi
                + observable_parameters["cos_phi_cos_psi"] * cos_phi * cos_psi
                + observable_parameters["cos_phi_sin_psi"] * cos_phi * sin_psi
                + observable_parameters["sin_phi_cos_psi"] * sin_phi * cos_psi
                + observable_parameters["sin_phi_sin_psi"] * sin_phi * sin_psi
                + observable_parameters["C"]
            ).m_as(unit.second**-1)

        else:
            # Get relevant dihedral angle
            if observable_parameters["dihedral"] == "prev_psi":
                dihedral_name = "psi"
                dihedral_resid = observable_resid - 1

            else:
                dihedral_name = observable_parameters["dihedral"]
                dihedral_resid = observable_resid

            karplus_df = dihedral_df[
                (dihedral_df["Dihedral Name"] == dihedral_name)
                & (dihedral_df["Resid"] == dihedral_resid)
            ]

            # Get residue specific parameters for sidechains
            if observable in {"3j_n_cg1", "3j_n_cg2", "3j_co_cg1", "3j_co_cg2"}:
                dihedral_resname = row["Resname"]
                observable_parameters = observable_parameters[dihedral_resname]

            elif observable in {"3j_ha_hb2", "3j_ha_hb3"}:
                dihedral_resname = PEREZ_KARPLUS_RESIDUE_MAP[row["Resname"]]
                observable_parameters = observable_parameters[dihedral_resname]

            # Compute cos(theta + delta)
            karplus_angle = numpy.deg2rad(
                observable_parameters["delta"] + karplus_df["Dihedral (deg)"].values
            )
            cos_angle = numpy.cos(karplus_angle)

            # Compute estimate for scalar coupling
            # <J> = A <cos^2(theta)> + B <cos(theta)> + C
            computed_coupling = (
                observable_parameters["A"] * numpy.square(cos_angle)
                + observable_parameters["B"] * cos_angle
                + observable_parameters["C"]
            ).m_as(unit.second**-1)

        # Get experimental uncertainty from Karplus model
        experiment_uncertainty = observable_parameters["sigma"].m_as(unit.second**-1)

        # Truncate experimental coupling to Karplus extrema
        experimental_coupling = row["Experiment"] / unit.second
        truncated_experimental_coupling = min(
            max(experimental_coupling, observable_parameters["minimum"]),
            observable_parameters["maximum"],
        ).m_as(unit.second**-1)

        if time_series_output_path is not None:
            # Get mean and SEM of correlated, truncated, and uncorrelated timeseries
            # for computed scalar coupling
            computed_coupling_mean = get_timeseries_mean(computed_coupling)

            # Write time series of observable
            observable_timeseries.append(
                {
                    "Frame": karplus_df["Frame"],
                    "Time (ns)": karplus_df["Time (ns)"],
                    "Observable": observable,
                    "Resid": observable_resid,
                    "Resname": row["Resname"],
                    "Experiment": row["Experiment"],
                    "Experiment Uncertainty": experiment_uncertainty,
                    "Truncated Experiment": truncated_experimental_coupling,
                    "Computed": computed_coupling,
                }
            )

        else:
            computed_coupling_mean = {
                "Correlated Mean": computed_coupling.mean(),
                "Correlated SEM": (
                    computed_coupling.std(ddof=1) / numpy.sqrt(computed_coupling.size)
                ),
            }

        # Write computed means of observable
        computed_observables.append(
            {
                "Experiment Uncertainty": experiment_uncertainty,
                "Truncated Experiment": truncated_experimental_coupling,
                **computed_coupling_mean,
            }
        )

    if time_series_output_path is not None:
        observable_timeseries_df = pandas.concat(
            [pandas.DataFrame(df) for df in observable_timeseries]
        ).reset_index(drop=True)
        observable_timeseries_df.to_csv(time_series_output_path)

    scalar_coupling_df = pandas.concat(
        [observable_df, pandas.DataFrame(computed_observables)], axis=1
    )
    scalar_coupling_df.to_csv(output_path)


def compute_h_bond_scalar_couplings(
    observable_path: str,
    h_bond_geometries_path: str,
    output_path: str,
    karplus: str = "barfield",
    time_series_output_path: str = None,
):
    """
    Compute NMR 3J_N_CO scalar couplings and chi^2 with respect to experimental
    values.

    Parameters
    ---------
    observable_path
        The path to the data for experimental observables.
    h_bond_geometries_path
        The path to the time series of hydrogen bond geometries.
    output_path
        The path to write the computed scalar couplings and chi^2 values.
    karplus
        The name of the set of Karplus parameters to use.
    time_series_output_path
        The path to write the time series of computed scalar couplings.
    """

    # Check Karplus parameters
    karplus = karplus.lower()

    if karplus == "barfield":
        karplus_parameters = BARFIELD_KARPLUS_PARAMETERS["3j_n_co"]
    else:
        raise ValueError("Argument `karplus` must be one of\n    barfield")

    # 3J(R, theta, phi) = exp(-a * (R - R_0))
    #     * ((A cos^2 phi + B cos phi + C) sin^2 theta + D cos^2 theta)
    # 3J(R, theta, phi) = (
    #     (exp(a R_0) A cos^2 phi + exp(a R_0) B cos phi + exp(a R_0) C)
    #     * sin^2 theta + exp(a R_0) D cos^2 theta) * exp(-a R)
    karplus_exp = karplus_parameters["exponent"]
    karplus_constant = numpy.exp(karplus_exp * karplus_parameters["min_dist"])
    karplus_cos_sq = karplus_constant * karplus_parameters["D"]
    karplus_sin_sq_cos_sq = karplus_constant * karplus_parameters["A"]
    karplus_sin_sq_cos = karplus_constant * karplus_parameters["B"]
    karplus_sin_sq = karplus_constant * karplus_parameters["C"]

    # Load data for experimental observables
    observable_df = pandas.read_csv(
        observable_path,
        sep="\s+",
        skiprows=1,
        names=[
            "Observable",
            "Resid N",
            "Resname N",
            "Resid CO",
            "Resname CO",
            "Experiment",
            "Uncertainty",
        ],
        usecols=[
            "Observable",
            "Resid N",
            "Resname N",
            "Resid CO",
            "Resname CO",
            "Experiment",
        ],
    )

    # Skip observables that don't have a good Karplus model
    indices_to_drop = list()
    for index, row in observable_df.iterrows():
        if row["Observable"] in {"3j_n_cg", "3j_n_cd"}:
            indices_to_drop.append(index)

    observable_df.drop(indices_to_drop, inplace=True)
    observable_df.reset_index(drop=True, inplace=True)

    # Read time series of hydrogen bond geometries for backbone amide H bonds
    chunk_size = 1E6
    h_bond_df = pandas.concat(
        [
            chunk[
                (chunk["Donor Name"] == "N")
                & (chunk["Hydrogen Name"] == "H")
                & (chunk["Acceptor Name"] == "O")
            ]
            for chunk in pandas.read_csv(
                h_bond_geometries_path, index_col=0, chunksize=chunk_size,
            )
        ]
    )

    # Compute observables
    if time_series_output_path is not None:
        observable_timeseries = list()
    computed_observables = list()

    for index, row in observable_df.iterrows():
        observable = row["Observable"]

        observable_resid_n = row["Resid N"]
        observable_resid_co = row["Resid CO"]

        geometry_df = h_bond_df[
            (h_bond_df["Donor Resid"] == observable_resid_n)
            & (h_bond_df["Acceptor Resid"] == observable_resid_co)
        ]

        HO_distance = geometry_df["HA Distance (Angstrom)"].values * unit.angstrom
        HOC_angle = numpy.deg2rad(geometry_df["HOC Angle (deg)"].values)
        HOCN_dihedral = numpy.deg2rad(geometry_df["HOCN Dihedral (deg)"].values)

        exp_HO_distance = numpy.exp(-karplus_exp * HO_distance)
        cos_sq_HOC_angle = numpy.square(numpy.cos(HOC_angle))
        sin_sq_HOC_angle = numpy.square(numpy.sin(HOC_angle))
        cos_HOCN_dihedral = numpy.cos(HOCN_dihedral)
        cos_sq_HOCN_dihedral = numpy.square(cos_HOCN_dihedral)

        # 3J(R, theta, phi) = (
        #     (exp(a R_0) A cos^2 phi + exp(a R_0) B cos phi + exp(a R_0) C)
        #     * sin^2 theta + exp(a R_0) D cos^2 theta) * exp(-a R)
        computed_coupling = (
            (
                (
                    karplus_sin_sq_cos_sq * cos_sq_HOCN_dihedral
                    + karplus_sin_sq_cos * cos_HOCN_dihedral
                    + karplus_sin_sq
                )
                * sin_sq_HOC_angle
                + karplus_cos_sq * cos_sq_HOC_angle
            )
            * exp_HO_distance
        ).m_as(unit.second**-1)

        # Get experimental uncertainty from Karplus model
        experiment_uncertainty = karplus_parameters["sigma"].m_as(unit.second**-1)

        if time_series_output_path is not None:
            # Get mean and SEM of correlated, truncated, and uncorrelated timeseries
            # for computed scalar coupling
            computed_coupling_mean = get_timeseries_mean(computed_coupling)

            # Write time series of observable
            observable_timeseries.append(
                {
                    "Frame": geometry_df["Frame"],
                    "Time (ns)": geometry_df["Time (ns)"],
                    "Observable": observable,
                    "Resid N": observable_resid_n,
                    "Resname N": row["Resname N"],
                    "Resid CO": observable_resid_co,
                    "Resname CO": row["Resname CO"],
                    "Experiment": row["Experiment"],
                    "Experiment Uncertainty": experiment_uncertainty,
                    "Computed": computed_coupling,
                }
            )

        else:
            computed_coupling_mean = {
                "Correlated Mean": computed_coupling.mean(),
                "Correlated SEM": (
                    computed_coupling.std(ddof=1) / numpy.sqrt(computed_coupling.size)
                ),
            }

        # Write computed means of observable
        computed_observables.append(
            {
                "Experiment Uncertainty": experiment_uncertainty,
                **computed_coupling_mean,
            }
        )

    if time_series_output_path is not None:
        observable_timeseries_df = pandas.concat(
            [pandas.DataFrame(df) for df in observable_timeseries]
        ).reset_index(drop=True)
        observable_timeseries_df.to_csv(time_series_output_path)

    scalar_coupling_df = pandas.concat(
        [observable_df, pandas.DataFrame(computed_observables)], axis=1
    )

    scalar_coupling_df.to_csv(output_path)


def compute_fraction_helix(
    observable_path: str,
    dihedral_clusters_path: str,
    h_bond_geometries_path: str,
    chemical_shifts_path: str,
    output_path: str,
    h_bond_distance_cutoff: unit.Quantity = 3.5 * unit.angstrom,
    h_bond_angle_cutoff: float = 30,
):
    """
    Compute fraction of helix by residue based on hydrogen bond occupancy.

    Parameters
    ---------
    observable_path
        The path to the data for experimental observables.
    dihedral_clusters_path
        The path to the time series of dihedral cluster assignments.
    h_bond_geometries_path
        The path to the time series of hydrogen bond geometries.
    chemical_shifts_path
        The path to the time series of chemical shifts.
    output_path
        The path to write the computed helical fractions.
    h_bond_distance_cutoff
        Distance between non-hydrogen donor and acceptor below which a hydrogen
        bond is considered to be occupied.
    h_bond_angle_cutoff
        Deviation in donor-hydrogen-acceptor angle from linear in degrees below
        which a hydrogen bond is considered to be occupied.
    """

    h_bond_distance_threshold = h_bond_distance_cutoff.m_as(unit.angstrom)
    h_bond_angle_threshold = 180 - h_bond_angle_cutoff

    # Load data for experimental observables
    observable_df = pandas.read_csv(
        observable_path,
        sep="\s+",
        skiprows=1,
        names=[
            "Resid", "Resname", "Experiment", "Shift Coil", "Delta Shift",
            "Delta Shift Residue", "SPARTA+ Alpha", "SPARTA+ PII",
            "ShiftX2 Alpha", "ShiftX2 PII",
        ],
    )

    chunk_size = 1E6

    # Read time series of dihedral cluster assignments
    dihedral_df = pandas.read_csv(dihedral_clusters_path, index_col=0)

    # Read time series of hydrogen bond geometries for backbone amide H bonds
    h_bond_df = pandas.concat(
        [
            chunk[
                (chunk["Donor Name"] == "N")
                & (chunk["Hydrogen Name"] == "H")
                & (chunk["Acceptor Name"] == "O")
            ]
            for chunk in pandas.read_csv(
                h_bond_geometries_path, index_col=0, chunksize=chunk_size,
            )
        ]
    )

    # Read time series of chemical shifts for backbone carbonyl carbons
    chemical_shift_df = pandas.concat(
        [
            chunk[chunk["Atom"] == "C"]
            for chunk in pandas.read_csv(
                chemical_shifts_path, index_col=0, chunksize=chunk_size,
            )
        ]
    )

    # Compute observables
    computed_observables = list()

    for index, row in observable_df.iterrows():
        observable_resid = row["Resid"]

        # observable_resid is 0-based, but some measured properties use 1-based
        # residue indices due to the capped initial structure
        capped_resid = observable_resid + 1
        residue_dihedral_df = dihedral_df[dihedral_df["Resid"] == capped_resid]

        # Get frames where (phi, psi) is closer than 30 deg to ideal alpha helix
        # at (-63, -43)
        helical_dihedrals = (
            (residue_dihedral_df["phi (deg)"] + 63) ** 2
            + (residue_dihedral_df["psi (deg)"] + 43) ** 2
        ) < 900

        computed_dihedral_fraction_helix = numpy.mean(helical_dihedrals)
        helical_dihedral_frames = residue_dihedral_df[helical_dihedrals]["Frame"]

        # Get hydrogen bonds to i-4 and i+4 residues
        residue_h_bond_df = h_bond_df[
            (
                (h_bond_df["Donor Resid"] == observable_resid)
                & (h_bond_df["Acceptor Resid"] == observable_resid - 4)
            )
            | (
                (h_bond_df["Donor Resid"] == observable_resid + 4)
                & (h_bond_df["Acceptor Resid"] == observable_resid)
            )
        ]

        if len(residue_h_bond_df) > 0:
            computed_h_bond_fraction_helix = numpy.mean(
                (
                    residue_h_bond_df["DA Distance (Angstrom)"]
                    < h_bond_distance_threshold
                )
                & (residue_h_bond_df["DHA Angle (deg)"] > h_bond_angle_threshold)
                & (residue_h_bond_df["Frame"].isin(helical_dihedral_frames))
            )

        else:
            computed_h_bond_fraction_helix = 0.0

        # Get mean chemical shift of backbone carbonyl carbon
        residue_chemical_shift_df = chemical_shift_df[
            chemical_shift_df["Resid"] == capped_resid
        ]
        carbonyl_chemical_shift = numpy.mean(
            residue_chemical_shift_df["Chemical Shift"]
        )

        # Compute fraction helix from carbonyl chemical shift using experimental
        # reference values for helix and coil
        coil_chemical_shift = row["Shift Coil"]
        delta_helix_coil_shift = row["Delta Shift Residue"]
        computed_chemical_shift_fraction_helix = (
            carbonyl_chemical_shift - coil_chemical_shift
        ) / delta_helix_coil_shift
        computed_chemical_shift_fraction_helix = numpy.clip(
            computed_chemical_shift_fraction_helix, 0.0, 1.0,
        )

        # Compute fraction helix from carbonyl chemical shift using ShiftX2
        # reference values for helix and coil
        #coil_chemical_shift = row["ShiftX2 PII"]
        #delta_helix_coil_shift = row["ShiftX2 Alpha"] - row["ShiftX2 PII"]
        coil_chemical_shift = row["SPARTA+ PII"]
        delta_helix_coil_shift = row["SPARTA+ Alpha"] - row["SPARTA+ PII"]
        computed_spartap_fraction_helix = (
            carbonyl_chemical_shift - coil_chemical_shift
        ) / delta_helix_coil_shift
        computed_spartap_fraction_helix = numpy.clip(
            computed_spartap_fraction_helix, 0.0, 1.0,
        )

        computed_observables.append(
            {
                "Computed Dihedrals": computed_dihedral_fraction_helix,
                "Computed H Bonds": computed_h_bond_fraction_helix,
                "Computed Chemical Shifts": computed_chemical_shift_fraction_helix,
                "Computed SPARTA+": computed_spartap_fraction_helix,
            }
        )

    helix_fraction_df = pandas.concat(
        [observable_df, pandas.DataFrame(computed_observables)], axis=1
    )

    helix_fraction_df.to_csv(output_path)
