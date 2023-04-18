import loos
from loos.pyloos import Trajectory
import numpy
from openff.toolkit import Molecule
from openmm import unit
import pandas
from pathlib import Path
from proteinbenchmark.analysis_parameters import *
from proteinbenchmark.benchmark_targets import (
    benchmark_targets,
    experimental_datasets,
)
from proteinbenchmark.utilities import list_of_dicts_to_csv
from typing import List


DEG_TO_RAD = numpy.pi / 180.0


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
    output_trajectory = loos.DCDWriter(f'{output_prefix}.dcd')

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

            with open(f'{output_prefix}.pdb', 'w') as pdb_file:
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

    # Select atoms for dihedrals
    dihedrals_by_residue = list()
    atom_selections = dict()

    for residue in topology.splitByResidue():

        resid = residue[0].resid()
        resname = residue[0].resname()

        residue_dihedrals = dict()

        for dihedral, dihedral_atom_dict in DIHEDRAL_ATOMS[resname].items():

            if resid == min_resid and dihedral == 'phi':
                continue

            if resid == max_resid and dihedral in ['psi', 'omega']:
                continue

            # Get atoms in dihedral
            atom_names = dihedral_atom_dict['atom_names']

            if 'resid_offsets' in dihedral_atom_dict:

                atom_resids = [
                    resid + offset
                    for offset in dihedral_atom_dict['resid_offsets']
                ]

            elif dihedral == 'tau':
                atom_resids = [resid] * 3

            else:
                atom_resids = [resid] * 4

            atom_list = list()

            for (atom_name, atom_resid) in zip(atom_names, atom_resids):

                atom_full_name = f'{atom_name}-{atom_resid}'

                # Atom selection is slow, so only do this once for each atom
                if atom_full_name not in atom_selections:

                    atom_selection = loos.selectAtoms(
                        topology,
                        f'resid == {atom_resid} && name == "{atom_name}"'
                    )

                    if len(atom_selection) == 0:

                        raise ValueError(
                            f'Unable to select atom {atom_name} with resid '
                            f'{atom_resid} for dihedral {dihedral}.'
                        )

                    atom_selections[atom_full_name] = atom_selection[0]

                atom_list.append(atom_selections[atom_full_name])

            residue_dihedrals[dihedral] = atom_list

        dihedrals_by_residue.append({
            'resid': resid,
            'resname': resname,
            'dihedrals': residue_dihedrals,
        })

    # Set up trajectory
    trajectory = Trajectory(trajectory_path, topology)
    frame_time = 0.0 * unit.picosecond
    fragment_index = 0
    dihedrals = list()

    # Load one frame into memory at a time
    for frame in trajectory:

        frame_time += frame_length

        frame_index = trajectory.index()
        frame_time_ns = frame_time.value_in_unit(unit.nanosecond)

        # Write dihedrals to file every 10 000 frames to avoid pandas
        # out-of-memory
        if frame_index % 10000 == 0 and frame_index > 0:

            list_of_dicts_to_csv(dihedrals, f'{output_path}-{fragment_index}')
            fragment_index += 1
            dihedrals = list()

        # Measure dihedrals
        for residue_dict in dihedrals_by_residue:
            for dihedral_name, atoms in residue_dict['dihedrals'].items():

                if len(atoms) == 3:
                    dihedral = loos.angle(atoms[0], atoms[1], atoms[2])

                else:

                    dihedral = loos.torsion(
                        atoms[0], atoms[1], atoms[2], atoms[3]
                    )

                dihedrals.append({
                    'Frame': frame_index,
                    'Time (ns)': frame_time_ns,
                    'Resid': residue_dict['resid'],
                    'Resname': residue_dict['resname'],
                    'Dihedral Name': dihedral_name,
                    'Dihedral (deg)': dihedral,
                })

    if fragment_index == 0:
        list_of_dicts_to_csv(dihedrals, output_path)

    else:
        list_of_dicts_to_csv(dihedrals, f'{output_path}-{fragment_index}')

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
        Distance between non-hydrogen donor and acceptor below which a hydrogen
        bond is considered to be occupied.
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
    offmol = Molecule.from_polymer_pdb(topology_path)
    putative_donors = offmol.chemical_environment_matches(
        '[#7,#8,#16:1]-[#1:2]'
    )
    putative_acceptors = offmol.chemical_environment_matches('[#7,#8,#16:1]')

    # Construct a list of [donor, hydrogen, acceptor] LOOS atom selections
    putative_h_bonds = list()
    atom_selections = dict()

    for (donor_index, hydrogen_index) in putative_donors:

        # Non-hydrogen atoms bonded to donor
        donor_bonded_atoms = [
            atom for atom in offmol.atoms[donor_index].bonded_atoms
            if atom.atomic_number == 1
        ]

        # Get LOOS atom selections for donor and hydrogen
        if donor_index not in atom_selections:

            donor_resid = offmol.atoms[donor_index].metadata['residue_number']
            donor_name = offmol.atoms[donor_index].name

            atom_selection = loos.selectAtoms(
                topology, f'resid == {donor_resid} && name == "{donor_name}"'
            )

            if len(atom_selection) == 0:

                raise ValueError(
                    f'Unable to select donor {donor_name} with resid '
                    f'{donor_resid}.'
                )

            atom_selections[donor_index] = atom_selection[0]

        if hydrogen_index not in atom_selections:

            hydrogen_resid = (
                offmol.atoms[hydrogen_index].metadata['residue_number']
            )
            hydrogen_name = offmol.atoms[hydrogen_index].name

            atom_selection = loos.selectAtoms(
                topology,
                f'resid == {hydrogen_resid} && name == "{hydrogen_name}"'
            )

            if len(atom_selection) == 0:

                raise ValueError(
                    f'Unable to select hydrogen {hydrogen_name} with resid '
                    f'{hydrogen_resid}.'
                )

            atom_selections[hydrogen_index] = atom_selection[0]

        donor_atom = atom_selections[donor_index]
        hydrogen_atom = atom_selections[hydrogen_index]

        # Find acceptors for this donor
        for (acceptor_index,) in putative_acceptors:

            # Don't consider acceptors within two covalent bonds of the donor
            if acceptor_index == donor_index or any([
                acceptor_index == bonded_atom.molecule_atom_index
                or offmol.atoms[acceptor_index].is_bonded_to(bonded_atom)
                for bonded_atom in donor_bonded_atoms
            ]):

                continue

            # Get LOOS atom selection for acceptor
            if acceptor_index not in atom_selections:

                acceptor_resid = (
                    offmol.atoms[acceptor_index].metadata['residue_number']
                )
                acceptor_name = offmol.atoms[acceptor_index].name

                atom_selection = loos.selectAtoms(
                    topology,
                    f'resid == {acceptor_resid} && name == "{acceptor_name}"'
                )

                if len(atom_selection) == 0:

                    raise ValueError(
                        f'Unable to select acceptor {acceptor_name} with resid '
                        f'{acceptor_resid}.'
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
            h_bond_distance_cutoff.value_in_unit(unit.angstrom),
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
            'donor': h_bond[0],
            'hydrogen': h_bond[1],
            'acceptor': h_bond[2],
        }
        for hb_index, h_bond in enumerate(putative_h_bonds)
        if h_bond_observations[hb_index] >= observation_threshold
    ]

    # Get atoms within two covalent bonds of acceptors in observed H bonds
    acceptor_bonded_atoms = dict()

    for h_bond in observed_h_bonds:

        acceptor_index = h_bond['acceptor'].index()

        # Keep track of this in a dictionary because a given acceptor may
        # participate in multiple observed hydrogen bonds
        if acceptor_index not in acceptor_bonded_atoms:

            acceptor_bonded_atoms[acceptor_index] = list()

            # Loop over non-hydrogen atoms covalently bonded to acceptor
            acceptor_bonded_indices = [
                atom.molecule_atom_index
                for atom in offmol.atoms[acceptor_index].bonded_atoms
                if atom.atomic_number != 1
            ]

            for bonded_index in acceptor_bonded_indices:

                # Get LOOS atom selection for bonded atom
                if bonded_index not in atom_selections:

                    atom_selection = loos.selectAtoms(
                        topology, f'index == {bonded_index}'
                    )

                    if len(atom_selection) == 0:
                        raise ValueError(f'No atom with index {bonded_index}')

                    atom_selections[bonded_index] = atom_selection[0]

                # Loop over non-hydrogen atoms two covalent bonds from the
                # acceptor
                acceptor_bonded_2_indices = [
                    atom.molecule_atom_index
                    for atom in offmol.atoms[bonded_index].bonded_atoms
                    if atom.atomic_number != 1
                    and atom.molecule_atom_index != acceptor_index
                ]

                bonded_2_atom_list = list()

                for bonded_2_index in acceptor_bonded_2_indices:

                    # Get LOOS atom selection for bonded atom
                    if bonded_2_index not in atom_selections:

                        atom_selection = loos.selectAtoms(
                            topology, f'index == {bonded_2_index}'
                        )

                        if len(atom_selection) == 0:

                            raise ValueError(
                                f'No atom with index {bonded_2_index}'
                            )

                        atom_selections[bonded_2_index] = atom_selection[0]

                    bonded_2_atom_list.append(atom_selections[bonded_2_index])

                acceptor_bonded_atoms[acceptor_index].append({
                    'bonded_atom': atom_selections[bonded_index],
                    'bonded_2_atoms': bonded_2_atom_list,
                })

        h_bond['acceptor_bonded_atoms'] = acceptor_bonded_atoms[acceptor_index]

    # Reset Trajectory iterator
    trajectory.reset()

    # Measure hydrogen bond geometries
    frame_time = 0.0 * unit.picosecond
    fragment_index = 0
    h_bond_geometries = list()

    for frame in trajectory:

        frame_time += frame_length

        frame_index = trajectory.index()
        frame_time_ns = frame_time.value_in_unit(unit.nanosecond)

        # Write hydrogen bond geomtries to file every 10 000 frames to avoid
        # pandas out-of-memory
        if frame_index % 10000 == 0 and frame_index > 0:

            list_of_dicts_to_csv(
                h_bond_geometries, f'{output_path}-{fragment_index}'
            )
            fragment_index += 1
            h_bond_geometries = list()

        # Measure h_bond_geometries
        for h_bond in observed_h_bonds:

            donor = h_bond['donor']
            hydrogen = h_bond['hydrogen']
            acceptor = h_bond['acceptor']

            donor_acceptor_distance = donor.coords().distance(acceptor.coords())
            hydrogen_acceptor_distance = (
                hydrogen.coords().distance(acceptor.coords())
            )

            donor_hydrogen_acceptor_angle = loos.angle(
                donor, hydrogen, acceptor
            )

            for acceptor_bonded_atoms in h_bond['acceptor_bonded_atoms']:

                bonded_atom = acceptor_bonded_atoms['bonded_atom']

                hydrogen_acceptor_bonded_angle = loos.angle(
                    hydrogen, acceptor, bonded_atom
                )

                for bonded_2_atom in acceptor_bonded_atoms['bonded_2_atoms']:

                    dihedral = loos.torsion(
                        hydrogen, acceptor, bonded_atom, bonded_2_atom
                    )

                    h_bond_geometries.append({
                        'Frame': frame_index,
                        'Time (ns)': frame_time_ns,
                        'Donor Resid': donor.resid(),
                        'Donor Resname': donor.resname(),
                        'Donor Name': donor.name(),
                        'Hydrogen Name': hydrogen.name(),
                        'Acceptor Resid': acceptor.resid(),
                        'Acceptor Resname': acceptor.resname(),
                        'Acceptor Name': acceptor.name(),
                        'Acceptor Bonded Name': bonded_atom.name(),
                        'Acceptor Bonded 2 Name': bonded_2_atom.name(),
                        'Donor Acceptor Distance (Angstrom)': (
                            donor_acceptor_distance
                        ),
                        'Hydrogen Acceptor Distance (Angstrom)': (
                            hydrogen_acceptor_distance
                        ),
                        'Donor Hydrogen Acceptor Angle (deg)': (
                            donor_hydrogen_acceptor_angle
                        ),
                        'Hydrogen Acceptor Bonded Angle (deg)': (
                            hydrogen_acceptor_bonded_angle
                        ),
                        'Dihedral (deg)': dihedral,
                    })

    # df['Occupied'] = (
    #     (df['Hydrogen Acceptor Distance (Angstrom)'] < 2.5)
    #     & (df['Donor Hydrogen Acceptor Angle (deg)'] > 150)
    # )
    # df.groupby(by = [
    #     'Donor Resid', 'Donor Name', 'Hydrogen Name', 'Acceptor Resid',
    #     'Acceptor Name'
    # ])['Occupied'].mean().reset_index()

    if fragment_index == 0:
        list_of_dicts_to_csv(h_bond_geometries, output_path)

    else:

        list_of_dicts_to_csv(
            h_bond_geometries, f'{output_path}-{fragment_index}'
        )

    return fragment_index


def assign_dihedral_clusters(
    dihedrals_path: str,
    output_path: str,
    ramachandran: str = 'hollingsworth',
    rotamer: str = 'hintze',
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

    if ramachandran == 'hollingsworth':
        ramachandran_clusters = HOLLINGSWORTH_RAMACHANDRAN_CLUSTERS

    else:

        raise ValueError(
            'Argument `ramachandran` must be one of\n    hollinsworth'
        )

    # Check rotamer library
    rotamer = rotamer.lower()

    if rotamer == 'hintze':
        rotamer_library = HINTZE_ROTAMER_LIBRARY

    else:
        raise ValueError('Argument `rotamer` must be one of\n    hintze')

    # Read time series of dihedrals
    dihedral_df = pandas.read_csv(dihedrals_path, index_col = 0)

    dihedral_df = dihedral_df.pivot(
        index = ['Frame', 'Time (ns)', 'Resid', 'Resname'],
        columns = 'Dihedral Name',
        values = 'Dihedral (deg)',
    ).add_suffix(' (deg)').reset_index().rename_axis(columns = None)

    # Wrap phi into [0, 360) and psi into [-120, 240) so that no clusters are
    # split across a periodic boundary
    # wrapped_value = (value - lower) % (upper - lower) + lower
    phi = dihedral_df['phi (deg)'] % 360
    psi = (dihedral_df['psi (deg)'] + 120) % 360 - 120

    # Assign residues with defined phi and psi values to the outlier cluster
    dihedral_df['Ramachandran Cluster'] = None
    dihedral_df.loc[
        ~((phi.isna()) | (psi.isna())), 'Ramachandran Cluster'
    ] = 'Outlier'

    # Assign Ramachandran clusters
    for cluster in ramachandran_clusters:

        dihedral_df.loc[
            (phi > cluster['phi'][0]) & (phi < cluster['phi'][1])
                & (psi > cluster['psi'][0]) & (psi < cluster['psi'][1]),
            'Ramachandran Cluster'
        ] = cluster['cluster']

    # Get list of sidechain dihedrals in target
    target_sidechain_dihedrals = [
        f'chi{i}' for i in range(1, 5) if f'chi{i} (deg)' in dihedral_df.columns
    ]

    # Assign sidechain rotamers
    dihedral_df['Sidechain Rotamer'] = None

    for resname in dihedral_df['Resname'].unique():

        if resname not in rotamer_library:
            continue

        # Assign residues with non-trivial sidechains to the outlier rotamer
        resname_rows = dihedral_df['Resname'] == resname
        dihedral_df.loc[resname_rows, 'Sidechain Rotamer'] = 'Outlier'

        # Assign rotamers for this residue type
        for rotamer_index, rotamer in rotamer_library[resname].iterrows():

            selected_rows = resname_rows

            for chi in ['chi1', 'chi2', 'chi3', 'chi4']:
                if f'{chi}_mean' in residue_df.columns:

                    diff = dihedral_df[f'{chi} (deg)'] - rotamer[f'{chi}_mean']
                    abs_wrapped_diff = ((diff + 180) % 360 - 180).abs()
                    selected_rows = selected_rows & (
                        abs_wrapped_diff < rotamer[f'{chi}_esd']
                    )

            dihedral_df.loc[
                selected_rows, 'Sidechain Rotamer'
            ] = rotamer['rotamer']

    dihedral_df.to_csv(output_path)


def compute_scalar_couplings(
    observable_path: str,
    dihedrals_path: str,
    output_path: str,
    karplus: str = 'best',
):
    """
    Compute NMR scalar couplings and chi^2 with respect to experimental values.

    Parameters
    ---------
    observable_path
        The path to the data for experimental observables.
    dihedrals_path
        The path to the time series of dihedrals.
    output_path
        The path to write the computed scalar couplings and chi^2 values.
    karplus
        The name of the set of Karplus parameters to use.
    """

    # Check Karplus parameters
    karplus = karplus.lower()

    if karplus == 'best':
        karplus_parameters = BEST_KARPLUS_PARAMETERS
    else:
        raise ValueError('Argument `karplus` must be one of\n    best')

    # Load data for experimental observables
    observable_df = pandas.read_csv(
        observable_path,
        sep = '\s+',
        skiprows = 1,
        names = ['Observable', 'Resid', 'Experiment', 'Uncertainty']
    )

    # Read time series of dihedrals
    dihedral_df = pandas.read_csv(dihedrals_path, index_col = 0)

    # Skip observables that depend on dihedrals that are undefined
    max_resid = dihedral_df['Resid'].max()
    indices_to_drop = list()
    for index, row in observable_df.iterrows():

        observable = row['Observable']
        observable_resid = row['Resid']
        observable_parameters = karplus_parameters[observable]
        observable_dihedrals = observable_parameters['dihedral'].split(',')

        if observable_resid == 1 and 'phi' in observable_dihedrals:
            indices_to_drop.append(index)

        if observable_resid == max_resid and 'psi' in observable_dihedrals:
            indices_to_drop.append(index)

    observable_df.drop(indices_to_drop, inplace = True)
    observable_df.reset_index(drop = True, inplace = True)

    # Compute observables
    trig_values = dict()
    trajectory_averages = dict()
    computed_observables = list()

    for index, row in observable_df.iterrows():

        observable = row['Observable']
        observable_resid = row['Resid']
        observable_parameters = karplus_parameters[observable]
        observable_dihedrals = observable_parameters['dihedral'].split(',')

        dihedral_names = list()

        for dihedral_name in observable_dihedrals:

            if dihedral_name == 'prev_psi':

                dihedral_name = 'psi'
                dihedral_resid = observable_resid - 1

            else:
                dihedral_resid = observable_resid

            dihedral_names.append((dihedral_name, dihedral_resid))

        if observable == '3j_hn_ca':

            # Get averages of trig functions of phi and psi
            karplus_cos_phi = observable_parameters['C_cos_phi']
            karplus_sin_phi = observable_parameters['C_sin_phi']
            karplus_cos_psi = observable_parameters['C_cos_psi']
            karplus_sin_psi = observable_parameters['C_sin_psi']
            karplus_cos_cos = observable_parameters['C_cos_phi_cos_psi']
            karplus_cos_sin = observable_parameters['C_cos_phi_sin_psi']
            karplus_sin_cos = observable_parameters['C_sin_phi_cos_psi']
            karplus_sin_sin = observable_parameters['C_sin_phi_sin_psi']
            karplus_C0 = observable_parameters['C_0']

            angle_names = list()

            for dihedral_name, dihedral_resid in dihedral_names:

                angle_name = f'{dihedral_name}-{dihedral_resid}'
                angle_names.append(angle_name)

                karplus_angle = DEG_TO_RAD * dihedral_df[
                    (dihedral_df['Dihedral Name'] == dihedral_name)
                    & (dihedral_df['Resid'] == dihedral_resid)
                ]['Dihedral (deg)']

                if f'Cos-{angle_name}' not in trajectory_averages:

                    cos_angle = karplus_angle.apply(numpy.cos)
                    cos_angle.reset_index(drop = True, inplace = True)
                    trig_values[f'Cos-{angle_name}'] = cos_angle
                    trajectory_averages[f'Cos-{angle_name}'] = cos_angle.mean()

                if f'Sin-{angle_name}' not in trajectory_averages:

                    sin_angle = karplus_angle.apply(numpy.sin)
                    sin_angle.reset_index(drop = True, inplace = True)
                    trig_values[f'Sin-{angle_name}'] = sin_angle
                    trajectory_averages[f'Sin-{angle_name}'] = sin_angle.mean()

            phi, psi = tuple(angle_names)

            for phi_name in [f'Cos-{phi}', f'Sin-{phi}']:
                for psi_name in [f'Cos-{psi}', f'Sin-{psi}']:

                    product_name = f'{phi_name}-{psi_name}'
                    product = trig_values[phi_name] * trig_values[psi_name]
                    trajectory_averages[product_name] = product.mean()

            # Compute estimate for scalar coupling
            computed_coupling = (
                karplus_cos_phi * trajectory_averages[f'Cos-{phi}']
                + karplus_sin_phi * trajectory_averages[f'Sin-{phi}']
                + karplus_cos_psi * trajectory_averages[f'Cos-{psi}']
                + karplus_sin_psi * trajectory_averages[f'Sin-{psi}']
                + karplus_cos_cos * trajectory_averages[f'Cos-{phi}-Cos-{psi}']
                + karplus_cos_sin * trajectory_averages[f'Cos-{phi}-Sin-{psi}']
                + karplus_sin_cos * trajectory_averages[f'Sin-{phi}-Cos-{psi}']
                + karplus_sin_sin * trajectory_averages[f'Sin-{phi}-Sin-{psi}']
                + karplus_C0
            )

        else:

            # Get averages of cos(theta + delta) and cos^2(theta + delta)
            dihedral_name, dihedral_resid = dihedral_names[0]
            karplus_A = observable_parameters['A']
            karplus_B = observable_parameters['B']
            karplus_C = observable_parameters['C']
            delta = observable_parameters['delta']

            angle_name = f'{dihedral_name}-{dihedral_resid}'
            if delta != 0.0:
                angle_name = f'{angle_name}+{delta:.1f}'

            if f'Cos^2-{angle_name}' not in trajectory_averages:

                karplus_angle = DEG_TO_RAD * (
                    dihedral_df[
                        (dihedral_df['Dihedral Name'] == dihedral_name)
                        & (dihedral_df['Resid'] == dihedral_resid)
                    ]['Dihedral (deg)'] + delta
                )

                cos_angle = karplus_angle.apply(numpy.cos)
                cos_angle.reset_index(drop = True, inplace = True)
                trig_values[f'Cos-{angle_name}'] = cos_angle
                trajectory_averages[f'Cos-{angle_name}'] = cos_angle.mean()

                cos_sq_angle = cos_angle.apply(numpy.square)
                trajectory_averages[f'Cos^2-{angle_name}'] = cos_sq_angle.mean()

            # Compute estimate for scalar coupling
            # <J> = A <cos^2(theta)> + B <cos(theta)> + C
            computed_coupling = (
                karplus_A * trajectory_averages[f'Cos^2-{angle_name}']
                + karplus_B * trajectory_averages[f'Cos-{angle_name}']
                + karplus_C
            )

        # Compute contribution to chi^2
        experimental_coupling = row['Experiment'] / unit.second
        uncertainty = observable_parameters['sigma']
        chi_sq = numpy.square(
            (computed_coupling - experimental_coupling) / uncertainty
        )

        computed_observables.append({
            'Computed': computed_coupling.value_in_unit(unit.second**-1),
            'Chi^2': chi_sq,
        })

    scalar_coupling_df = pandas.concat(
        [observable_df, pandas.DataFrame(computed_observables)],
        axis = 1
    )

    scalar_coupling_df.to_csv(output_path)

