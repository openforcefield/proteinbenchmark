import loos
from loos.pyloos import Trajectory
import numpy
from openmm import unit
import pandas
from pathlib import Path
from proteinbenchmark.analysis_parameters import *
from proteinbenchmark.benchmark_targets import (
    benchmark_targets,
    experimental_datasets,
)
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
    Measure backbone and sidechain dihedrals over the trajectory.

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

    dihedrals = list()
    frame_time = 0.0 * unit.picosecond

    # Load one frame into memory at a time
    for frame in trajectory:

        frame_time += frame_length

        frame_index = trajectory.index()
        frame_time_ns = frame_time.value_in_unit(unit.nanosecond)

        # Measure dihedrals
        for residue_dict in dihedrals_by_residue:
            for dihedral_name, atoms in residue_dict['dihedrals'].items():

                dihedral = loos.torsion(atoms[0], atoms[1], atoms[2], atoms[3])

                dihedrals.append({
                    'Frame': frame_index,
                    'Time (ns)': frame_time_ns,
                    'Resid': residue_dict['resid'],
                    'Resname': residue_dict['resname'],
                    'Dihedral Name': dihedral_name,
                    'Dihedral (deg)': dihedral,
                })

    dihedral_df = pandas.DataFrame(dihedrals)
    dihedral_df.to_csv(output_path)


def assign_ramachandran_cluster(
    dihedrals_path: str,
    output_path: str,
):
    """
    Assign frames to Ramachandran clusters based on backbone dihedrals.

    Parameters
    ---------
    dihedrals_path
        Path to the time series of dihedrals.
    output_path
        Path to write Ramachandran cluster assignments.
    """

    # Read time series of dihedrals
    dihedrals_df = pandas.read_csv(dihedrals_path, index_col = 0)

    # Ramachandran clusters
    


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

    if karplus.lower() == 'best':
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

    scalar_coupling_df.to_csv(output_path, sep = ' ')

