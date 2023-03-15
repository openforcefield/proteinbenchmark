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


def reimage_trajectory(
    topology_path: str,
    trajectory_path: str,
    output_prefix: str,
    output_selection: str = None,
    align_selection: str = None,
    reference_path: str = None,
    reference_selection: str = None,
):
    """
    Center, reimage, and align a subset of atoms in a trajectory.

    Parameters
    ---------
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    output_prefix
        The prefix for the path to write the reimaged topology and trajectory.
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
        #output_atoms.applyTransform(transform_matrix)

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


def compute_scalar_couplings(
    observable_path: str,
    topology_path: str,
    trajectory_path: str,
    frame_length: unit.Quantity,
    output_prefix: str,
    karplus: str = 'best',
):
    """
    Compute NMR scalar couplings and chi^2 with respect to experimental values.

    Parameters
    ---------
    observable_path
        The path to the data for experimental observables.
    topology_path
        The path to the system topology, e.g. a PDB file.
    trajectory_path
        The path to the trajectory.
    frame_length
       The amount of simulation time between frames in the trajectory. 
    output_prefix
        The path to write the time series of dihedrals and computed scalar
        couplings and chi^2.
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

    # Load topology
    topology = loos.createSystem(topology_path)
    max_resid = topology.maxResid()

    # Skip 1j_n_ca scalar couplings for C terminal residues because psi is
    # undefined
    indices_to_drop = list()
    for index, row in observable_df.iterrows():

        observable = row['Observable']
        observable_resid = row['Resid']
        observable_name = f'{observable}-{observable_resid}'

        if observable_name == f'1j_n_ca-{max_resid}':
            observable_df.drop(index, inplace = True)

    observable_df.drop(indices_to_drop, inplace = True)
    observable_df.reset_index(drop = True, inplace = True)

    # Select atoms for observables
    observable_dihedral_names = dict()
    dihedral_atoms = dict()
    atom_selections = dict()

    for index, row in observable_df.iterrows():

        observable = row['Observable']
        observable_resid = row['Resid']
        observable_name = f'{observable}-{observable_resid}'
        observable_parameters = karplus_parameters[observable]
        observable_dihedrals = observable_parameters['dihedral'].split(',')

        observable_dihedral_names[observable_name] = list()

        for dihedral in observable_dihedrals:

            if dihedral == 'prev_psi':

                dihedral = 'psi'
                dihedral_resid = observable_resid - 1

            else:
                dihedral_resid = observable_resid

            dihedral_name = f'{dihedral}-{dihedral_resid}'
            observable_dihedral_names[observable_name].append(dihedral_name)

            # Get atoms in dihedral
            if dihedral_name in dihedral_atoms:
                continue

            atom_names = DIHEDRAL_ATOMS[dihedral]['atom_names']
            atom_resids = [
                dihedral_resid + offset
                for offset in DIHEDRAL_ATOMS[dihedral]['resid_offsets']
            ]

            atom_list = list()

            for (name, resid) in zip(atom_names, atom_resids):

                # Atom selection is slow, so only do this once for each atom
                if f'{name}-{resid}' not in atom_selections:

                    atom_selection = loos.selectAtoms(
                        topology, f'resid == {resid} && name == "{name}"'
                    )

                    if len(atom_selection) == 0:

                        raise ValueError(
                            f'Unable to select atom {name} with resid {resid} '
                            f'for dihedral {dihedral_name} and observable '
                            f'{observable_name}.'
                        )

                    atom_selections[f'{name}-{resid}'] = atom_selection[0]

                atom_list.append(atom_selections[f'{name}-{resid}'])

            dihedral_atoms[dihedral_name] = atom_list


    # Set up trajectory
    trajectory = Trajectory(trajectory_path, topology)

    dihedrals = list()
    frame_time = 0.0 * unit.picosecond

    # Load one frame into memory at a time
    for frame in trajectory:

        frame_time += frame_length

        frame_dihedrals = {
            'Frame': trajectory.index(),
            'Time(ns)': frame_time.value_in_unit(unit.nanosecond),
        }

        # Measure dihedrals
        for dihedral_name, atoms in dihedral_atoms.items():

            dihedral = loos.torsion(atoms[0], atoms[1], atoms[2], atoms[3])
            frame_dihedrals[f'{dihedral_name}(deg)'] = dihedral

        dihedrals.append(frame_dihedrals)

    dihedral_df = pandas.DataFrame(dihedrals)

    dihedral_df.to_csv(f'{output_prefix}-dihedrals.dat', sep = ' ')

    # Compute observables
    trajectory_averages = dict()
    computed_observables = list()

    for index, row in observable_df.iterrows():

        observable = row['Observable']
        observable_resid = row['Resid']
        observable_name = f'{observable}-{observable_resid}'
        dihedral_names = observable_dihedral_names[observable_name]

        observable_parameters = karplus_parameters[observable]

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
            karplus_C = observable_parameters['C_0']

            for angle_name in dihedral_names:

                karplus_angle = DEG_TO_RAD * dihedral_df[f'{angle_name}(deg)']

                if f'Cos-{angle_name}' not in trajectory_averages:

                    cos_angle = karplus_angle.apply(numpy.cos)
                    dihedral_df[f'Cos-{angle_name}'] = cos_angle
                    trajectory_averages[f'Cos-{angle_name}'] = cos_angle.mean()

                if f'Sin-{angle_name}' not in trajectory_averages:

                    sin_angle = karplus_angle.apply(numpy.sin)
                    dihedral_df[f'Sin-{angle_name}'] = sin_angle
                    trajectory_averages[f'Sin-{angle_name}'] = sin_angle.mean()

            phi, psi = tuple(dihedral_names)

            for phi_name in [f'Cos-{phi}', f'Sin-{phi}']:
                for psi_name in [f'Cos-{psi}', f'Sin-{psi}']:

                    product_name = f'{phi_name}-{psi_name}'
                    product = dihedral_df[phi_name] * dihedral_df[psi_name]
                    dihedral_df[product_name] = product
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
                + karplus_C
            )

        else:

            # Get averages of cos(theta + delta) and cos^2(theta + delta)
            dihedral_name = dihedral_names[0]
            karplus_A = observable_parameters['A']
            karplus_B = observable_parameters['B']
            karplus_C = observable_parameters['C']
            delta = observable_parameters['delta']

            angle_name = dihedral_name
            if delta != 0.0:
                angle_name = f'{dihedral_name}+{delta:.1f}'

            if f'Cos^2-{angle_name}' not in trajectory_averages:

                karplus_angle = DEG_TO_RAD * (
                    dihedral_df[f'{dihedral_name}(deg)'] + delta
                )

                cos_angle = karplus_angle.apply(numpy.cos)
                dihedral_df[f'Cos-{angle_name}'] = cos_angle
                trajectory_averages[f'Cos-{angle_name}'] = cos_angle.mean()

                cos_sq_angle = cos_angle.apply(numpy.square)
                dihedral_df[f'Cos^2-{angle_name}'] = cos_sq_angle
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

    scalar_coupling_df.to_csv(f'{output_prefix}-observables.dat', sep = ' ')

