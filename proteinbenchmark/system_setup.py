import numpy
from openff.toolkit import (
    ForceField as OFFForceField,
    Molecule as OFFMolecule,
    Topology as OFFTopology,
)
from openmm import app, unit
import openmm
from pathlib import Path
from proteinbenchmark.utilities import (
    read_xml,
    remove_model_lines,
    write_pdb,
    write_xml
)


class _OFFForceField(OFFForceField):
    """
    Dummy class that defines a `createSystem()` method so that this force field
    can be passed to the `addSolvent() and `_addIons()` methods of
    `openmm.app.Modeller`.
    """

    def __init__(self, unique_molecules, *args, **kwargs):

        self._from_openmm_unique_molecules = unique_molecules
        super().__init__(*args, **kwargs)

    def createSystem(self, openmm_topology: app.topology.Topology):
        """Return an OpenMM system from an OpenMM topology."""

        off_topology = OFFTopology.from_openmm(
            openmm_topology,
            unique_molecules = self._from_openmm_unique_molecules,
        )

        return self.create_openmm_system(off_topology)


def build_initial_coordinates(
    build_method: str,
    ph: float,
    initial_pdb: str,
    protonated_pdb: str,
    aa_sequence: str = None,
    nterm_cap: str = None,
    cterm_cap: str = None,
):
    """
    Build initial coordinates and set protonation state.

    Parameters
    ---------
    build_method
        The method for building the initial coordinates. "extended" will build
        coordinates from a given amino acid sequence using pmx with phi and psi
        set to -180 deg for every residue. "pdb" will build coordinates from a
        given PDB file.
    ph
        The pH used to assign protonation state.
    initial_pdb
        The path to the initial PDB file. Input for the "pdb" build method,
        output for the "extended" build method.
    protonated_pdb
        The path to write the protonated PDB with initial coordinates.
    aa_sequence
        The primary amino acid sequence for the "extended" build method.
    nterm_cap
        The capping group to add on the N terminus. Must be "ace" or None.
    cterm_cap
        The capping group to add on the C terminus. Must be "nh2", "nme", or
        None.
    """

    from pdb2pqr.main import build_main_parser as pdb2pqr_build_main_parser
    from pdb2pqr.main import main_driver as pdb2pqr_main_driver
    import pmx

    # Get initial PDB from file or from pmx
    if build_method == 'pdb':

        print(
            f'Building initial coordinates at pH {ph:.2f} using build_method '
            f'"pdb" with initial_pdb\n    {initial_pdb}'
        )

    elif build_method in {'extended', 'helical'}:

        if aa_sequence is None:

            raise ValueError(
                'aa_sequence must be set to use build_method "extended" or '
                "helical"
            )

        print(
            f'Building initial coordinates at pH {ph:.3f} using build_method '
            f'"{build_method}" with aa_sequence\n    {aa_sequence}'
        )

        # Build sequence in pmx
        chain = pmx.Chain().create(aa_sequence)

        # Set backbone dihedrals
        if build_method == 'helical':

            # Alpha helix values from
            # Hollingsworth SA, Karplus PA. (2010). BioMol Concepts 1, 271-283.
            build_phi = -63
            build_psi = -43

        else:

            # Extended conformation with all backbone dihedrals at -180 deg
            build_phi = -180
            build_psi = -180

        for residue in chain.residues:

            residue.set_phi(build_phi, propagate = True)
            residue.set_psi(build_psi, propagate = True)

        # Add N terminal cap
        if nterm_cap is not None:

            nterm_cap = nterm_cap.lower()

            if nterm_cap == 'ace':
                chain.add_nterm_cap()

            else:
                raise ValueError('Argument `nterm_cap` must be one of\n    ace')

        # Add C terminal cap
        if cterm_cap is not None:

            cterm_cap = cterm_cap.lower()

            if cterm_cap == 'nme':
                chain.add_cterm_cap()

            elif cterm_cap == 'nh2':

                # pmx only supports Nme caps, so add the Nhe cap manually in a
                # similar way to pmx.Chain.add_cterm_cap()
                chain.cbuild('GLY')
                cterm = chain.cterminus()

                for atom in ['O', 'C', 'HA1', 'HA2']:
                    del cterm[atom]

                n, h, ca = cterm.fetchm(['N', 'H', 'CA'])
                h.name = 'HN1'
                ca.name = 'HN2'

                n_to_h = numpy.array(n.x) - numpy.array(h.x)
                n_to_ca = numpy.array(n.x) - numpy.array(ca.x)
                n_h_dist = numpy.linalg.norm(n_to_h)
                n_ca_dist = numpy.linalg.norm(n_to_ca)

                ca.x = numpy.array(n.x) - n_to_ca * n_h_dist / n_ca_dist

                cterm.set_resname('NH2')

            else:

                raise ValueError(
                    'Argument `cterm_cap` must be one of\n    nh2\n    nme'
                )

        # Write pmx system to PDB file
        chain.write(initial_pdb)

        # PMX writes MODEL lines that cannot be read by OpenMM
        remove_model_lines(initial_pdb)

    else:

        raise ValueError(
            'Argument `build_method` must be one of\n    "extended"'
            '\n    "helical"\n    "pdb"'
        )

    # Clean up initial structure and assign protonation states using pdb2pqr

    # Use PARSE force field because it supports more protonation states than
    # AMBER
    cterm_pka = 3.2
    nterm_pka = 8.0

    if ph < cterm_pka and cterm_cap is not None:
        pdb2pqr_ff = 'PARSE --neutralc'
    elif ph > nterm_pka and nterm_cap is not None:
        pdb2pqr_ff = 'PARSE --neutraln'
    else:
        pdb2pqr_ff = 'PARSE'

    protonated_pqr = f'{protonated_pdb[:-4]}.pqr'
    pdb2pqr_args = (
        f'--ff {pdb2pqr_ff} --nodebump --titration-state-method propka '
        f'--with-ph {ph:.1f} -o {ph:.1f} --pdb-output {protonated_pdb} '
        f'{initial_pdb} {protonated_pqr}'
    )
    pdb2pqr_parser = pdb2pqr_build_main_parser()
    pdb2pqr_main_driver(pdb2pqr_parser.parse_args(pdb2pqr_args.split()))

    # Special handling for protonated C terminus
    if ph < cterm_pka and cterm_cap is not None:

        protonated_model = pmx.Model(protonated_pdb)
        c_term_residue = protonated_model.residues[-1]

        # pdb2pqr adds the C-terminal carboxyl hydrogen to the oxygen
        # farther from the carbon. The OpenMM residue template expects the
        # carboxyl hydrogen to be bonded to the "OXT" oxygen. If "O" is
        # farther from "C" than "OXT", swap "O" and "OXT" coordinates
        c = c_term_residue['C']
        o = c_term_residue['O']
        oxt = c_term_residue['OXT']

        if o - c > oxt - c:

            print(
                'Swapping coordinates of "O" and "OXT" in residue '
                f'{o.resname} {o.resnr:d}'
            )

            oxt_to_o = numpy.array(o.x) - numpy.array(oxt.x)
            oxt.translate(oxt_to_o)
            o.translate(-oxt_to_o)

            protonated_model.write(protonated_pdb)

            # Remove MODEL lines written by PMX that cannot be read by OpenMM
            # Also remove END so that we can write CONECT records
            remove_model_lines(protonated_pdb, remove_end = True)

        # Add CONECT record to enforce bond between "HO" and "OXT"
        ho = c_term_residue['HO']
        with open(protonated_pdb, 'a') as pdb_file:

            pdb_file.write(f'CONECT{oxt.id:5d}{ho.id:5d}\n')
            pdb_file.write(f'CONECT{ho.id:5d}{oxt.id:5d}\n')
            pdb_file.write('END')


def solvate(
    solvent_padding: unit.Quantity,
    ionic_strength: unit.Quantity,
    nonbonded_cutoff: unit.Quantity,
    vdw_switch_width: unit.Quantity,
    protonated_pdb_file: str,
    solvated_pdb_file: str,
    openmm_system_xml: str,
    water_model: str,
    force_field_file: str,
    water_model_file: str = None,
):
    """
    Add water and salt ions and write OpenMM System to XML.

    Parameters
    ----------
    solvent_padding
        The padding distance used to setup the solvent box.
    ionic_strength
        The bulk ionic strength of the desired thermodynamic state.
    nonbonded_cutoff
        The cutoff for the Lennard-Jones potential and PME direct space
        summation.
    vdw_switch_width
        The distance from the nonbonded cutoff at which to apply the
        switching function.
    protonated_pdb_file
        The path to the protonated PDB with initial coordinates.
    solvated_pdb_file
        The path to write the solvated PDB with ions.
    openmm_system_xml
        The path to write the parametrized OpenMM system as a serialized XML.
    water_model
        The name of the water model used to parametrize the water.
    force_field_file
        The path to the force field to parametrize the system.
    water_model_file
        The path to the force field containing the water model.
    """

    # Check units of arguments
    if not solvent_padding.unit.is_compatible(unit.nanometer):
        raise ValueError('solvent_padding does not have units of Length')

    if not ionic_strength.unit.is_compatible(unit.molar):

        raise ValueError(
            'ionic_strength does not have units of Amount Length^-3'
        )

    if not nonbonded_cutoff.unit.is_compatible(unit.nanometer):
        raise ValueError('nonbonded_cutoff does not have units of Length')
    if not vdw_switch_width.unit.is_compatible(unit.nanometer):
        raise ValueError('vdw_switch_width does not have units of Length')

    print(
        f'Solvating system with water model {water_model} and'
        '\n    solvent_padding '
        f'{solvent_padding.value_in_unit(unit.nanometer):.3f} nm'
        '\n    ionic_strength '
        f'{ionic_strength.value_in_unit(unit.molar):.3f} M'
        '\n    nonbonded_cutoff '
        f'{nonbonded_cutoff.value_in_unit(unit.nanometer):.3f} nm'
        '\n    vdw_switch_width '
        f'{vdw_switch_width.value_in_unit(unit.nanometer):.3f} nm'
    )

    # Set up force field
    smirnoff = Path(force_field_file).suffix == '.offxml'

    if smirnoff:

        # SMIRNOFF force field

        # Get unique molecules (solute, water, sodium ion, chloride ion) needed
        # for openff.toolkit.topology.Topology.from_openmm()
        solute_offmol = OFFMolecule.from_polymer_pdb(protonated_pdb_file)
        unique_molecules = [
            solute_offmol,
            OFFMolecule.from_smiles('O'),
            OFFMolecule.from_smiles('[Na+1]'),
            OFFMolecule.from_smiles('[Cl-1]'),
        ]

        # Use the dummy class _OFFForceField so we can pass this to Modeller
        force_field = _OFFForceField(unique_molecules, force_field_file)
        print(f'Force field read from\n    {force_field_file}')

        # Set up solute topology and positions
        solute_off_topology = solute_offmol.to_topology()
        solute_interchange = force_field.create_interchange(solute_off_topology)
        solute_topology = solute_interchange.topology.to_openmm()
        solute_positions = solute_interchange.positions.to_openmm()

    else:

        if water_model_file is None:

            # OpenMM force field with no separate water model
            force_field = app.ForceField(force_field_file)
            print(f'Force field read from\n    {force_field_file}')

        else:

            # OpenMM force field with separate water model
            force_field = app.ForceField(force_field_file, water_model_file)
            print(
                f'Force field read from\n    {force_field_file}'
                f'\n    and {water_model_file}'
            )

        # Set up solute topology and positions
        solute_pdb = app.PDBFile(protonated_pdb_file)
        solute_topology = solute_pdb.topology
        solute_positions = solute_pdb.positions

    # Add water in a rhombic dodecahedral box with no ions
    modeller = app.Modeller(solute_topology, solute_positions)

    modeller.addSolvent(
        forcefield = force_field,
        model = water_model,
        padding = solvent_padding,
        boxShape = 'dodecahedron',
        ionicStrength = 0 * unit.molar,
        neutralize = False,
    )

    # Add salt ions using the SLTCAP method

    # Get total charge of system without ions
    system = force_field.createSystem(modeller.topology)

    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), openmm.NonbondedForce):

            nonbonded_force = system.getForce(i)
            break

    else:
        raise ValueError('The ForceField does not specify a NonbondedForce')

    total_charge = unit.Quantity(0, unit.elementary_charge)
    for i in range(nonbonded_force.getNumParticles()):
        total_charge += nonbonded_force.getParticleParameters(i)[0]

    # Round to nearest integer
    total_charge = int(
        numpy.round(total_charge.value_in_unit(unit.elementary_charge))
    )

    print(f'Total charge is {total_charge} e')

    # Get the number of water molecules and their positions
    water_positions = {}
    _oxygen = app.element.oxygen
    for chain in modeller.topology.chains():
        for residue in chain.residues():
            if residue.name == 'HOH':
                for atom in residue.atoms():
                    if atom.element == _oxygen:
                        water_positions[residue] = (
                            modeller.positions[atom.index]
                        )

    n_water = len(water_positions)

    print(f'Added {n_water} waters')

    # Calculate effective ionic strength from desired bulk ionic strength
    # using SLTCAP
    # N+ = V_w C_0 (sqrt(1 + Q^2 / (2 V_w C_0)^2) - Q / (2 V_w C_0))
    # N- = V_w C_0 (sqrt(1 + Q^2 / (2 V_w C_0)^2) + Q / (2 V_w C_0))
    # N+ = V_w C_eff + (|Q| - Q) / 2
    # N- = V_w C_eff + (|Q| + Q) / 2
    # N+ + N- = 2 V_w C_0 sqrt(1 + Q^2 / (2 V_w C_0)^2) = 2 V_w C_eff + |Q|
    # C_eff = C_0 sqrt(1 + Q^2 / (2 V_w C_0)^2) - |Q| / (2 V_w C_0))

    # If ionic strength is zero, then Q / (2 V_w C_0) is undefined
    if ionic_strength.value_in_unit(unit.molar) > 0:

        bulk_water_concentration = 55.4 * unit.molar

        charge_magnitude = numpy.abs(total_charge)
        solvent_volume = (n_water - charge_magnitude) / bulk_water_concentration
        solute_ion_ratio = charge_magnitude / (
            2 * solvent_volume * ionic_strength
        )
        sltcap_effective_ionic_strength = ionic_strength * (
            numpy.sqrt(1 + solute_ion_ratio * solute_ion_ratio)
            - solute_ion_ratio
        )

        # Compute number of ions expected for neutral solute and for SLTCAP
        n_cation_expected_neutral = int(numpy.round(
            solvent_volume * ionic_strength - total_charge / 2
        ))
        n_anion_expected_neutral = int(numpy.round(
            solvent_volume * ionic_strength + total_charge / 2
        ))
        n_cation_expected_sltcap = int(numpy.round(
            solvent_volume * ionic_strength
            * numpy.sqrt(1 + solute_ion_ratio * solute_ion_ratio)
            - total_charge / 2
        ))
        n_anion_expected_sltcap = int(numpy.round(
            solvent_volume * ionic_strength
            * numpy.sqrt(1 + solute_ion_ratio * solute_ion_ratio)
            + total_charge / 2
        ))

        print(
            f'Solute-to-ion ratio |Q| / (2 e V_w C_0): {solute_ion_ratio:f}'
            '\nEffective ionic strength: '
            f'{sltcap_effective_ionic_strength.value_in_unit(unit.molar):f} M'
            '\nExpected number of ions for neutral solute: '
            f'{int(numpy.round(n_cation_expected_neutral)):d} Na+ '
            f'{int(numpy.round(n_anion_expected_neutral)):d} Cl-'
            '\nExpected number of ions for SLTCAP: '
            f'{int(numpy.round(n_cation_expected_sltcap)):d} Na+ '
            f'{int(numpy.round(n_anion_expected_sltcap)):d} Cl-'
        )

        # Add ions using modeller with effective ionic strength to achieve
        # expected number of ion pairs from SLTCAP
        modeller._addIons(
            force_field, n_water, water_positions,
            ionicStrength = sltcap_effective_ionic_strength,
            neutralize = True,
        )

    elif total_charge != 0:

        # Neutralize if there are no bulk ions
        modeller._addIons(
            force_field, n_water, water_positions,
            ionicStrength = ionic_strength,
            neutralize = True,
        )

    # Report the number of ions added
    n_cation, n_anion = 0, 0
    _sodium_residue = app.element.sodium.symbol.upper()
    _chlorine_residue = app.element.chlorine.symbol.upper()
    for chain in modeller.topology.chains():
        for residue in chain.residues():
            if residue.name == _sodium_residue:
                n_cation += 1
            elif residue.name == _chlorine_residue:
                n_anion += 1

    print(f'Actual number of ions added: {n_cation} Na+ {n_anion} Cl-')

    # Write solvated system to PDB file
    write_pdb(solvated_pdb_file, modeller.topology, modeller.positions)

    # Create an OpenMM System from the solvated system
    if smirnoff:

        openmm_system = force_field.createSystem(modeller.topology)

    else:

        switch_distance = nonbonded_cutoff - vdw_switch_width
        openmm_system = force_field.createSystem(
            modeller.topology,
            nonbondedMethod = app.PME,
            nonbondedCutoff = nonbonded_cutoff,
            switchDistance = switch_distance,
            constraints = app.HBonds,
        )

    # Validate total charge of solvated system
    for i in range(openmm_system.getNumForces()):
        if isinstance(openmm_system.getForce(i), openmm.NonbondedForce):

            nonbonded_force = openmm_system.getForce(i)
            break

    else:
        raise ValueError('The ForceField does not specify a NonbondedForce')

    solvated_total_charge = unit.Quantity(0, unit.elementary_charge)
    for i in range(nonbonded_force.getNumParticles()):
        solvated_total_charge += nonbonded_force.getParticleParameters(i)[0]

    solvated_total_charge = int(
        numpy.round(solvated_total_charge.value_in_unit(unit.elementary_charge))
    )

    if solvated_total_charge != 0:

        raise ValueError(
            f'Total charge of solvated system is {solvated_total_charge:d}'
        )

    # Write OpenMM system to XML file
    write_xml(openmm_system_xml, openmm_system)


def minimize(
    restraint_energy_constant: unit.Quantity,
    openmm_system_xml: str,
    solvated_pdb_file: str,
    minimized_pdb_file: str,
):
    """
    Minimize energy of solvated system with Cartesian restraints on non-hydrogen
    solute atoms.

    Parameters
    ----------
    restraint_energy_constant
        Energy constant for Cartesian restraints (units Energy Length^-2).
    openmm_system_xml
        The path to the parametrized OpenMM system as a serialized XML.
    solvated_pdb_file
        The path to the solvated PDB with ions.
    minimized_pdb_file
        The path to write the minimized PDB.
    """

    # Check units of arguments
    if not restraint_energy_constant.unit.is_compatible(
        unit.kilojoules_per_mole / unit.nanometer**2
    ):

        raise ValueError(
            'restraint_energy_constant does not have units of Energy Length^-2'
        )

    k = restraint_energy_constant.value_in_unit(
        unit.kilocalories_per_mole / unit.angstrom**2
    )

    print(
        f'Minimizing energy with Cartesian restraint energy constant {k:.4f} '
        'kcal mol^-1 angstrom^-2'
    )

    # Load OpenMM system and solvated PDB
    openmm_system = read_xml(openmm_system_xml)
    solvated_pdb = app.PDBFile(solvated_pdb_file)

    # Create Cartesian restraints on non-hydrogen solute atoms
    cartesian_restraint = openmm.CustomExternalForce(
        'k * periodicdistance(x, y, z, x0, y0, z0)^2'
    )
    openmm_system.addForce(cartesian_restraint)
    cartesian_restraint.addGlobalParameter('k', restraint_energy_constant)
    cartesian_restraint.addPerParticleParameter('x0')
    cartesian_restraint.addPerParticleParameter('y0')
    cartesian_restraint.addPerParticleParameter('z0')

    _hydrogen = openmm.app.element.hydrogen
    _sodium_residue = app.element.sodium.symbol.upper()
    _chlorine_residue = app.element.chlorine.symbol.upper()
    solvent_residue_names = ['HOH', _sodium_residue, _chlorine_residue]
    for chain in solvated_pdb.topology.chains():
        for residue in chain.residues():
            if residue.name not in solvent_residue_names:
                for atom in residue.atoms():
                    if atom.element != _hydrogen:
                        cartesian_restraint.addParticle(
                            atom.index, solvated_pdb.positions[atom.index]
                        )

    # Set up minimization and print initial energy
    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    simulation = app.Simulation(
        solvated_pdb.topology,
        openmm_system,
        integrator,
        openmm.Platform.getPlatformByName('CUDA'),
    )
    simulation.context.setPositions(solvated_pdb.positions)
    initial_state = simulation.context.getState(getEnergy = True)
    initial_energy = initial_state.getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole
    )
    print(
        f'Initial energy of solvated system: {initial_energy:.4f} kcal mol^-1'
    )

    # Run minimization with Cartesian restraints and print final energy
    simulation.minimizeEnergy()
    final_state = simulation.context.getState(
        getEnergy = True, getPositions = True
    )
    final_energy = final_state.getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole
    )
    print(f'Final energy of minimized system: {final_energy:.4f} kcal mol^-1')

    # Write minimized coordinates to PDB
    write_pdb(
        minimized_pdb_file,
        simulation.topology,
        final_state.getPositions(),
    )

