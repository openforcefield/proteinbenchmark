from pathlib import Path

import numpy
import openmm
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from openmm import app
from openmm import unit as openmm_unit
from openmmforcefields.generators import EspalomaTemplateGenerator

from proteinbenchmark.force_fields import force_fields
from proteinbenchmark.utilities import (
    read_xml,
    remove_model_lines,
    write_pdb,
    write_xml,
)


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

    import pmx
    from pdb2pqr.main import build_main_parser as pdb2pqr_build_main_parser
    from pdb2pqr.main import main_driver as pdb2pqr_main_driver

    # Get initial PDB from file or from pmx
    if build_method == "pdb":
        print(
            f"Building initial coordinates at pH {ph:.2f} using build_method "
            f'"pdb" with initial_pdb\n    {initial_pdb}'
        )

    elif build_method in {"extended", "helical"}:
        if aa_sequence is None:
            raise ValueError(
                'aa_sequence must be set to use build_method "extended" or ' "helical"
            )

        print(
            f"Building initial coordinates at pH {ph:.3f} using build_method "
            f'"{build_method}" with aa_sequence\n    {aa_sequence}'
        )

        # Build sequence in pmx
        chain = pmx.Chain().create(aa_sequence)

        # Set backbone dihedrals
        if build_method == "helical":
            # Alpha helix values from
            # Hollingsworth SA, Karplus PA. (2010). BioMol Concepts 1, 271-283.
            build_phi = -63
            build_psi = -43

        else:
            # Extended conformation with all backbone dihedrals at -180 deg
            build_phi = -180
            build_psi = -180

        for residue in chain.residues:
            residue.set_phi(build_phi, propagate=True)
            residue.set_psi(build_psi, propagate=True)

        # Add N terminal cap
        if nterm_cap is not None:
            nterm_cap = nterm_cap.lower()

            if nterm_cap == "ace":
                chain.add_nterm_cap()

            else:
                raise ValueError("Argument `nterm_cap` must be one of\n    ace")

        # Add C terminal cap
        if cterm_cap is not None:
            cterm_cap = cterm_cap.lower()

            if cterm_cap == "nme":
                chain.add_cterm_cap()

            elif cterm_cap == "nh2":
                # pmx only supports Nme caps, so add the Nhe cap manually in a
                # similar way to pmx.Chain.add_cterm_cap()
                chain.cbuild("GLY")
                cterm = chain.cterminus()

                for atom in ["O", "C", "HA1", "HA2"]:
                    del cterm[atom]

                n, h, ca = cterm.fetchm(["N", "H", "CA"])
                h.name = "HN1"
                ca.name = "HN2"

                n_to_h = numpy.array(h.x) - numpy.array(n.x)
                n_to_ca = numpy.array(ca.x) - numpy.array(n.x)
                n_h_dist = numpy.linalg.norm(n_to_h)
                n_ca_dist = numpy.linalg.norm(n_to_ca)

                ca.x = numpy.array(n.x) + n_to_ca * n_h_dist / n_ca_dist

                cterm.set_resname("NH2")

            else:
                raise ValueError(
                    "Argument `cterm_cap` must be one of\n    nh2\n    nme"
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

    if ph < cterm_pka and cterm_cap is None:
        pdb2pqr_ff = "PARSE --neutralc"
    elif ph > nterm_pka and nterm_cap is None:
        pdb2pqr_ff = "PARSE --neutraln"
    else:
        pdb2pqr_ff = "PARSE"

    protonated_pqr = f"{protonated_pdb[:-4]}.pqr"
    pdb2pqr_args = (
        f"--ff {pdb2pqr_ff} --nodebump --titration-state-method propka "
        f"--with-ph {ph:.1f} -o {ph:.1f} --pdb-output {protonated_pdb} "
        f"{initial_pdb} {protonated_pqr}"
    )
    pdb2pqr_parser = pdb2pqr_build_main_parser()
    pdb2pqr_main_driver(pdb2pqr_parser.parse_args(pdb2pqr_args.split()))

    # Special handling for protonated C terminus
    if ph < cterm_pka and cterm_cap is None:
        protonated_model = pmx.Model(protonated_pdb)
        c_term_residue = protonated_model.residues[-1]

        # pdb2pqr adds the C-terminal carboxyl hydrogen to the oxygen
        # farther from the carbon. The OpenMM residue template expects the
        # carboxyl hydrogen to be bonded to the "OXT" oxygen. If "O" is
        # farther from "C" than "OXT", swap "O" and "OXT" coordinates
        c = c_term_residue["C"]
        o = c_term_residue["O"]
        oxt = c_term_residue["OXT"]

        if o - c > oxt - c:
            print(
                'Swapping coordinates of "O" and "OXT" in residue '
                f"{o.resname} {o.resnr:d}"
            )

            oxt_to_o = numpy.array(o.x) - numpy.array(oxt.x)
            oxt.translate(oxt_to_o)
            o.translate(-oxt_to_o)

            protonated_model.write(protonated_pdb)

            # Remove MODEL lines written by PMX that cannot be read by OpenMM
            # Also remove END so that we can write CONECT records
            remove_model_lines(protonated_pdb, remove_end=True)

        # Add CONECT record to enforce bond between "HO" and "OXT"
        ho = c_term_residue["HO"]
        with open(protonated_pdb, "a") as pdb_file:
            pdb_file.write(f"CONECT{oxt.id:5d}{ho.id:5d}\n")
            pdb_file.write(f"CONECT{ho.id:5d}{oxt.id:5d}\n")
            pdb_file.write("END")


def solvate(
    ionic_strength: unit.Quantity,
    protonated_pdb_file: str,
    solvated_pdb_file: str,
    water_model: str,
    solvent_padding: unit.Quantity = None,
    n_solvent: int = None,
):
    """
    Add water and salt ions. Exactly one of solvent_padding or n_solvent must be
    specified.

    Parameters
    ----------
    ionic_strength
        The bulk ionic strength of the desired thermodynamic state.
    protonated_pdb_file
        The path to the protonated PDB with initial coordinates.
    solvated_pdb_file
        The path to write the solvated PDB with ions.
    water_model
        The name of the water model used to parametrize the water.
    solvent_padding
        The padding distance used to setup the solvent box.
    n_solvent
        The number of solvent molecules used to setup the solvent box.
    """

    # Check arguments
    if (solvent_padding is None) == (n_solvent is None):
        raise ValueError(
            "Exactly one of solvent_padding or n_solvent must be specified."
        )

    if solvent_padding is not None and not solvent_padding.is_compatible_with(
        unit.nanometer
    ):
        raise ValueError("solvent_padding does not have units of Length")

    if not ionic_strength.is_compatible_with(unit.molar):
        raise ValueError("ionic_strength does not have units of Amount Length^-3")

    if solvent_padding is not None:
        solvent_arg_str = (
            f"\n    solvent_padding {solvent_padding.m_as(unit.nanometer):.3f} nm"
        )

    else:
        solvent_arg_str = f"\n    n_solvent {n_solvent:d}"

    print(
        f"Solvating system with water model {water_model} and"
        f"{solvent_arg_str}"
        "\n    ionic_strength "
        f"{ionic_strength.m_as(unit.molar):.3f} M"
    )

    if water_model in {"opc3", "tip3p", "tip3p-fb"}:
        water_has_virtual_sites = False
        modeller_water_model = "tip3p"
    elif water_model in {"opc", "tip4p-fb"}:
        water_has_virtual_sites = True
        modeller_water_model = "tip4pew"
    else:
        raise ValueError(
            "water_model must be one of\n    opc\n    opc3\n    tip3p"
            "\n    tip3p-fb\n    tip4p-fb"
        )

    # Use Amber ff14SB as a reference for building solvent coordinates
    force_field = app.ForceField(
        force_fields["ff14sb-tip3p"]["force_field_file"],
        force_fields["ff14sb-tip3p"]["water_model_file"],
    )

    # Set up solute topology and positions
    solute_pdb = app.PDBFile(protonated_pdb_file)
    modeller = app.Modeller(solute_pdb.topology, solute_pdb.positions)

    # Add water in a rhombic dodecahedral box with no ions
    if solvent_padding is not None:
        modeller.addSolvent(
            forcefield=force_field,
            model=modeller_water_model,
            padding=solvent_padding.to_openmm(),
            boxShape="dodecahedron",
            ionicStrength=0 * openmm_unit.molar,
            neutralize=False,
        )

    else:
        modeller.addSolvent(
            forcefield=force_field,
            numAdded=n_solvent,
            model=modeller_water_model,
            boxShape="dodecahedron",
            ionicStrength=0 * openmm_unit.molar,
            neutralize=False,
        )

        # Set dodecahedron box vectors correctly
        from openmm.vec3 import Vec3

        box_width = modeller.topology.getPeriodicBoxVectors()[0][0].value_in_unit(
            openmm_unit.nanometer
        )
        box_vectors = (
            Vec3(box_width, 0, 0),
            Vec3(1 / 3, 2 * numpy.sqrt(2) / 3, 0) * box_width,
            Vec3(-1 / 3, numpy.sqrt(2) / 3, numpy.sqrt(6) / 3) * box_width,
        )
        modeller.topology.setPeriodicBoxVectors(box_vectors * openmm_unit.nanometer)

    # Add salt ions using the SLTCAP method

    # Get total charge of system without ions
    system = force_field.createSystem(modeller.topology)

    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nonbonded_force = force
            break

    else:
        raise ValueError("The ForceField does not specify a NonbondedForce")

    total_charge = openmm_unit.Quantity(0, openmm_unit.elementary_charge)
    for i in range(nonbonded_force.getNumParticles()):
        total_charge += nonbonded_force.getParticleParameters(i)[0]

    # Round to nearest integer
    total_charge = int(
        numpy.round(total_charge.value_in_unit(openmm_unit.elementary_charge))
    )

    print(f"Total charge is {total_charge} e")

    # Get the number of water molecules and their positions
    water_positions = {}
    _oxygen = app.element.oxygen
    for chain in modeller.topology.chains():
        for residue in chain.residues():
            if residue.name == "HOH":
                for atom in residue.atoms():
                    if atom.element == _oxygen:
                        water_positions[residue] = modeller.positions[atom.index]

    n_water = len(water_positions)

    print(f"Added {n_water} waters")

    # Calculate effective ionic strength from desired bulk ionic strength
    # using SLTCAP
    # N+ = V_w C_0 (sqrt(1 + Q^2 / (2 V_w C_0)^2) - Q / (2 V_w C_0))
    # N- = V_w C_0 (sqrt(1 + Q^2 / (2 V_w C_0)^2) + Q / (2 V_w C_0))
    # N+ = V_w C_eff + (|Q| - Q) / 2
    # N- = V_w C_eff + (|Q| + Q) / 2
    # N+ + N- = 2 V_w C_0 sqrt(1 + Q^2 / (2 V_w C_0)^2) = 2 V_w C_eff + |Q|
    # C_eff = C_0 sqrt(1 + Q^2 / (2 V_w C_0)^2) - |Q| / (2 V_w C_0))

    # If ionic strength is zero, then Q / (2 V_w C_0) is undefined
    if ionic_strength.m_as(unit.molar) > 0:
        bulk_water_concentration = 55.4 * unit.molar

        charge_magnitude = numpy.abs(total_charge)
        solvent_volume = (n_water - charge_magnitude) / bulk_water_concentration
        solute_ion_ratio = charge_magnitude / (2 * solvent_volume * ionic_strength)
        sltcap_effective_ionic_strength = ionic_strength * (
            numpy.sqrt(1 + solute_ion_ratio * solute_ion_ratio) - solute_ion_ratio
        )

        # Compute number of ions expected for neutral solute and for SLTCAP
        n_cation_expected_neutral = int(
            numpy.round(solvent_volume * ionic_strength - total_charge / 2)
        )
        n_anion_expected_neutral = int(
            numpy.round(solvent_volume * ionic_strength + total_charge / 2)
        )
        n_cation_expected_sltcap = int(
            numpy.round(
                solvent_volume
                * ionic_strength
                * numpy.sqrt(1 + solute_ion_ratio * solute_ion_ratio)
                - total_charge / 2
            )
        )
        n_anion_expected_sltcap = int(
            numpy.round(
                solvent_volume
                * ionic_strength
                * numpy.sqrt(1 + solute_ion_ratio * solute_ion_ratio)
                + total_charge / 2
            )
        )

        print(
            f"Solute-to-ion ratio |Q| / (2 e V_w C_0): {solute_ion_ratio:f}"
            "\nEffective ionic strength: "
            f"{sltcap_effective_ionic_strength.m_as(unit.molar):f} M"
            "\nExpected number of ions for neutral solute: "
            f"{int(numpy.round(n_cation_expected_neutral)):d} Na+ "
            f"{int(numpy.round(n_anion_expected_neutral)):d} Cl-"
            "\nExpected number of ions for SLTCAP: "
            f"{int(numpy.round(n_cation_expected_sltcap)):d} Na+ "
            f"{int(numpy.round(n_anion_expected_sltcap)):d} Cl-"
        )

        # Add ions using modeller with effective ionic strength to achieve
        # expected number of ion pairs from SLTCAP
        modeller._addIons(
            force_field,
            n_water,
            water_positions,
            ionicStrength=sltcap_effective_ionic_strength.to_openmm(),
            neutralize=True,
        )

    elif total_charge != 0:
        # Neutralize if there are no bulk ions
        modeller._addIons(
            force_field,
            n_water,
            water_positions,
            ionicStrength=ionic_strength.to_openmm(),
            neutralize=True,
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

    print(f"Actual number of ions added: {n_cation} Na+ {n_anion} Cl-")

    # Write solvated system to PDB file
    write_pdb(solvated_pdb_file, modeller.topology, modeller.positions)


def assign_parameters(
    simulation_platform: str,
    nonbonded_cutoff: unit.Quantity,
    vdw_switch_width: unit.Quantity,
    protonated_pdb_file: str,
    solvated_pdb_file: str,
    parametrized_system: str,
    water_model: str,
    force_field_file: str,
    water_model_file: str = None,
    hydrogen_mass: unit.Quantity = 3.0 * unit.dalton,
    setup_prefix: str = None,
):
    """
    Assign parameters to the solvated topology.

    Parameters
    ----------
    simulation_platform
        Simulation platform for file exporting.
    nonbonded_cutoff
        The cutoff for the Lennard-Jones potential and PME direct space
        summation.
    vdw_switch_width
        The distance from the nonbonded cutoff at which to apply the
        switching function.
    protonated_pdb_file
        The path to the protonated PDB with initial coordinates.
    solvated_pdb_file
        The path to the solvated PDB with ions.
    parametrized_system
        The path to write the parametrized system.
    water_model
        The name of the water model used to parametrize the water.
    force_field_file
        The path to the force field to parametrize the system.
    water_model_file
        The path to the force field containing the water model.
    hydrogen_mass
        The mass of solute hydrogen atoms for hydrogen mass repartitioning.
    """

    # Check arguments
    if not nonbonded_cutoff.is_compatible_with(unit.nanometer):
        raise ValueError("nonbonded_cutoff does not have units of Length")
    if not vdw_switch_width.is_compatible_with(unit.nanometer):
        raise ValueError("vdw_switch_width does not have units of Length")
    if not hydrogen_mass.is_compatible_with(unit.dalton):
        raise ValueError("hydrogen_mass does not have units of Mass")

    solvated_pdb = app.PDBFile(solvated_pdb_file)
    openmm_topology = solvated_pdb.topology

    # Check bond between OXT and HO for protonated C terminus
    for chain in openmm_topology.chains():
        for residue in chain.residues():
            # Check for both OXT and HO atoms in the same residue
            oxt = None
            ho = None
            for atom in residue.atoms():
                if atom.name == "OXT":
                    oxt = atom
                elif atom.name == "HO":
                    ho = atom
            if oxt is not None and ho is not None:
                # Check that there is a bond between these atoms
                for bond in openmm_topology.bonds():
                    if (bond[0] == oxt and bond[1] == ho) or (
                        bond[0] == ho and bond[1] == oxt
                    ):
                        break
                else:
                    # Create a bond between these atoms
                    openmm_topology.addBond(oxt, ho)

    # Set up force field
    ff_extension = Path(force_field_file).suffix
    if ff_extension == ".offxml":
        ff_type = "smirnoff"
        ff_class = ForceField

    elif ff_extension == ".pt":
        ff_type = "espaloma"

        # Make a force field with the water model, then add espaloma using a
        # template generator from openmmforcefields
        ff_class = app.ForceField
        espaloma_model_file = force_field_file
        force_field_file = water_model_file
        water_model_file = None

    else:
        ff_type = "openmm"
        ff_class = app.ForceField

    if water_model_file is None:
        force_field = ff_class(force_field_file)
        print(f"Force field read from\n    {force_field_file}")

    else:
        force_field = ff_class(force_field_file, water_model_file)
        print(
            f"Force field read from\n    {force_field_file}"
            f"\n    and {water_model_file}"
        )

    # Create the parametrized system from the solvated topology
    if ff_type in {"smirnoff", "espaloma"}:
        # Get list of OpenFF Molecules in solute OpenMM topology
        openff_solute_topology = Topology.from_pdb(protonated_pdb_file)
        openff_solute_molecules = openff_solute_topology.unique_molecules

    if ff_type == "smirnoff":
        # Get unique molecules (solute, water, sodium ion, chloride ion) needed
        # for openff.toolkit.topology.Topology.from_openmm()
        unique_molecules = [
            *openff_solute_molecules,
            Molecule.from_smiles("O"),
            Molecule.from_smiles("[Na+1]"),
            Molecule.from_smiles("[Cl-1]"),
        ]

        if water_model in {"opc", "tip4p-fb"}:
            # Create a new OpenMM topology without water virtual sites
            no_vsite_topology = app.Topology()
            no_vsite_topology.setPeriodicBoxVectors(
                openmm_topology.getPeriodicBoxVectors()
            )
            no_vsite_atoms = dict()

            for chain in openmm_topology.chains():
                no_vsite_chain = no_vsite_topology.addChain(chain.id)

                for residue in chain.residues():
                    no_vsite_residue = no_vsite_topology.addResidue(
                        residue.name,
                        no_vsite_chain,
                        residue.id,
                        residue.insertionCode,
                    )

                    for atom in residue.atoms():
                        if residue.name != "HOH" or atom.name in {"O", "H1", "H2"}:
                            no_vsite_atom = no_vsite_topology.addAtom(
                                atom.name, atom.element, no_vsite_residue
                            )
                            no_vsite_atoms[atom] = no_vsite_atom

            # Include bonds only between non-virtual site atoms
            for bond in openmm_topology.bonds():
                if bond[0] in no_vsite_atoms and bond[1] in no_vsite_atoms:
                    no_vsite_topology.addBond(
                        no_vsite_atoms[bond[0]], no_vsite_atoms[bond[1]]
                    )

            openff_topology = Topology.from_openmm(
                no_vsite_topology,
                unique_molecules=unique_molecules,
            )

        else:
            openff_topology = Topology.from_openmm(
                openmm_topology,
                unique_molecules=unique_molecules,
            )

        openmm_system = force_field.create_openmm_system(openff_topology)

        # Repartition hydrogen mass to bonded heavy atom
        if simulation_platform == "openmm":
            openmm_hydrogen_mass = hydrogen_mass.to_openmm()
            # Manually change hydrogen masses in OpenMM system. Taken from
            # https://github.com/openmm/openmm/blob/f30d716ace8331003c5115bdfa9e03341a757878/wrappers/python/openmm/app/forcefield.py#L1249
            _hydrogen = app.element.hydrogen
            for atom1, atom2 in modeller.topology.bonds():
                if atom1.element == _hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if (
                    atom2.element == _hydrogen
                    and atom1.element not in {_hydrogen, None}
                    and atom2.residue.name != "HOH"
                ):
                    transfer_mass = (
                        openmm_hydrogen_mass
                        - openmm_system.getParticleMass(atom2.index)
                    )
                    heavy_mass = (
                        openmm_system.getParticleMass(atom1.index) - transfer_mass
                    )
                    openmm_system.setParticleMass(atom2.index, openmm_hydrogen_mass)
                    openmm_system.setParticleMass(atom1.index, heavy_mass)

    else:
        if ff_type == "espaloma":
            # Set up template generator with espaloma model
            espaloma_generator = EspalomaTemplateGenerator(
                molecules=openff_solute_molecules,
                forcefield=espaloma_model_file,
            )
            force_field.registerTemplateGenerator(espaloma_generator.generator)

            # Get OpenMM topology of solute with one residue per molecule
            openmm_solute_topology = openff_solute_topology.to_openmm()
            espaloma_topology = app.Topology()
            espaloma_topology.setPeriodicBoxVectors(
                openmm_topology.getPeriodicBoxVectors()
            )
            solute_atoms = dict()
            chain_index = 0

            for chain in openmm_solute_topology.chains():
                new_chain = espaloma_topology.addChain(chain.id)
                residue_name = f"XX{chain_index:01d}"
                chain_index += 1
                resid = "1"
                new_residue = espaloma_topology.addResidue(
                    residue_name, new_chain, resid
                )

                for residue in chain.residues():
                    for atom in residue.atoms():
                        new_atom = espaloma_topology.addAtom(
                            atom.name, atom.element, new_residue, atom.id
                        )
                        solute_atoms[atom] = new_atom

            for bond in openmm_solute_topology.bonds():
                if bond[0] in solute_atoms and bond[1] in solute_atoms:
                    espaloma_topology.addBond(
                        solute_atoms[bond[0]], solute_atoms[bond[1]]
                    )

            # Copy solvent topology from solvated system topology
            solute_chain_indices = set(
                chain.index for chain in espaloma_topology.chains()
            )
            solvent_atoms = dict()

            for chain in openmm_topology.chains():
                if chain.index in solute_chain_indices:
                    continue
                new_chain = espaloma_topology.addChain(chain.id)
                for residue in chain.residues():
                    new_residue = espaloma_topology.addResidue(
                        residue.name, new_chain, residue.id
                    )
                    for atom in residue.atoms():
                        new_atom = espaloma_topology.addAtom(
                            atom.name, atom.element, new_residue, atom.id
                        )
                        solvent_atoms[atom] = new_atom

            for bond in openmm_topology.bonds():
                if bond[0] in solvent_atoms and bond[1] in solvent_atoms:
                    espaloma_topology.addBond(
                        solvent_atoms[bond[0]], solvent_atoms[bond[1]]
                    )

            openmm_topology = espaloma_topology

        # Create parametrized system
        switch_distance = nonbonded_cutoff - vdw_switch_width
        if simulation_platform == "openmm":
            openmm_system = force_field.createSystem(
                openmm_topology,
                nonbondedMethod=app.PME,
                nonbondedCutoff=nonbonded_cutoff.to_openmm(),
                constraints=app.HBonds,
                rigidWater=True,
                hydrogenMass=hydrogen_mass.to_openmm(),
                switchDistance=switch_distance.to_openmm(),
            )
        elif simulation_platform == "gmx":
            openmm_system = force_field.createSystem(
                openmm_topology,
                nonbondedMethod=app.PME,
                nonbondedCutoff=nonbonded_cutoff.to_openmm(),
                rigidWater=False,
                hydrogenMass=hydrogen_mass.to_openmm(),
                switchDistance=switch_distance.to_openmm(),
            )

    if (
        ff_type == "smirnoff" and simulation_platform == "openmm"
    ) or ff_type != "smirnoff":
        # Validate total charge of solvated system
        for force in openmm_system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                nonbonded_force = force
                break

        else:
            raise ValueError("The ForceField does not specify a NonbondedForce")

        solvated_total_charge = openmm_unit.Quantity(0, openmm_unit.elementary_charge)
        for i in range(nonbonded_force.getNumParticles()):
            solvated_total_charge += nonbonded_force.getParticleParameters(i)[0]

        solvated_total_charge = int(
            numpy.round(
                solvated_total_charge.value_in_unit(openmm_unit.elementary_charge)
            )
        )

        if solvated_total_charge != 0:
            raise ValueError(
                f"Total charge of solvated system is {solvated_total_charge:d}"
            )

    if simulation_platform == "openmm":
        # Write OpenMM system to XML file
        write_xml(parametrized_system, openmm_system)

    elif simulation_platform == "gmx":
        # TODO write GMX files with Interchange
        import parmed as pmd

        # Write GROMACS files
        struct = pmd.openmm.load_topology(
            openmm_topology, openmm_system, xyz=modeller.positions
        )
        hmass = pmd.tools.HMassRepartition(struct, hydrogen_mass.m_as(unit.dalton))
        hmass.execute()

        struct.save(str(setup_prefix) + ".gro")
        struct.save(parametrized_system)

        # Add position restraints file to topology
        itp_file = f"{setup_prefix.name}_posre.itp"
        setup_split = str(setup_prefix).split("/")
        match_string = "[ moleculetype ]"
        insert_string = f'#ifdef POSRES\n#include "{itp_file}"\n#endif\n'

        mol = 0
        with open(parametrized_system, "r+") as fd:
            contents = fd.readlines()
            # Handle last line to prevent IndexError
            if match_string in contents[-1]:
                contents.append(insert_string)
            else:
                for index, line in enumerate(contents):
                    if match_string in line:
                        if mol == 1 and insert_string not in contents[index - 1]:
                            contents.insert(index - 1, insert_string)
                            break
                        else:
                            mol = 1
            fd.seek(0)
            fd.writelines(contents)
        print("GROMACS Files Printed")


def minimize_openmm(
    parametrized_system: str,
    solvated_pdb_file: str,
    minimized_coords_file: str,
    restraint_energy_constant: unit.Quantity,
    force_tolerance: unit.Quantity,
):
    """
    Minimize energy of solvated system with Cartesian restraints on non-hydrogen
    solute atoms using OpenMM.

    Parameters
    ----------
    parametrized_system
        The path to the parametrized system (OpenMM XML or GMX TOP).
    solvated_pdb_file
        The path to the solvated PDB with ions.
    minimized_coords_file
        The path to write the minimized coords (OpenMM PDB or GMX GRO).
    restraint_energy_constant
        Energy constant for Cartesian restraints (units Energy Length^-2).
    force_tolerance
        Force tolerance for minimization (units Energy Length^-1).
    """

    # Check units of arguments
    if not restraint_energy_constant.is_compatible_with(
        unit.kilojoule_per_mole / unit.nanometer**2
    ):
        raise ValueError(
            "restraint_energy_constant does not have units of Energy Length^-2"
        )
    if not force_tolerance.is_compatible_with(unit.kilojoule_per_mole / unit.nanometer):
        raise ValueError("force_tolerance does not have units of Energy Length^-1")

    k = restraint_energy_constant.m_as(unit.kilocalories_per_mole / unit.angstrom**2)

    print(
        f"Minimizing energy with Cartesian restraint energy constant {k:.4f} "
        "kcal mol^-1 angstrom^-2"
    )

    # Load OpenMM system and solvated PDB
    openmm_system = read_xml(parametrized_system)
    solvated_pdb = app.PDBFile(solvated_pdb_file)

    # Create Cartesian restraints on non-hydrogen solute atoms
    cartesian_restraint = openmm.CustomExternalForce(
        "k * periodicdistance(x, y, z, x0, y0, z0)^2"
    )
    openmm_system.addForce(cartesian_restraint)
    cartesian_restraint.addGlobalParameter("k", restraint_energy_constant.to_openmm())
    cartesian_restraint.addPerParticleParameter("x0")
    cartesian_restraint.addPerParticleParameter("y0")
    cartesian_restraint.addPerParticleParameter("z0")

    _hydrogen = openmm.app.element.hydrogen
    _sodium_residue = app.element.sodium.symbol.upper()
    _chlorine_residue = app.element.chlorine.symbol.upper()
    solvent_residue_names = ["HOH", _sodium_residue, _chlorine_residue]
    for chain in solvated_pdb.topology.chains():
        for residue in chain.residues():
            if residue.name not in solvent_residue_names:
                for atom in residue.atoms():
                    if atom.element != _hydrogen:
                        cartesian_restraint.addParticle(
                            atom.index, solvated_pdb.positions[atom.index]
                        )

    # Set up minimization and print initial energy
    integrator = openmm.VerletIntegrator(1.0 * openmm_unit.femtosecond)
    simulation = app.Simulation(
        solvated_pdb.topology,
        openmm_system,
        integrator,
        openmm.Platform.getPlatformByName("CUDA"),
    )
    simulation.context.setPositions(solvated_pdb.positions)
    initial_state = simulation.context.getState(getEnergy=True)
    initial_energy = initial_state.getPotentialEnergy().value_in_unit(
        openmm_unit.kilocalorie_per_mole
    )
    print(f"Initial energy of solvated system: {initial_energy:.4f} kcal mol^-1")

    # Run minimization with Cartesian restraints and print final energy
    simulation.minimizeEnergy()
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = final_state.getPotentialEnergy().value_in_unit(
        openmm_unit.kilocalorie_per_mole
    )
    print(f"Final energy of minimized system: {final_energy:.4f} kcal mol^-1")

    # Write minimized coordinates to PDB
    write_pdb(
        minimized_coords_file,
        simulation.topology,
        final_state.getPositions(),
    )


def minimize_gmx(
    parametrized_system: str,
    solvated_pdb_file: str,
    minimized_coords_file: str,
    setup_prefix: str,
    gmx_executable: str,
    nonbonded_cutoff: unit.Quantity,
    vdw_switch_width: unit.Quantity,
    restraint_energy_constant: unit.Quantity,
    force_tolerance: unit.Quantity,
):
    """
    Minimize energy of solvated system with Cartesian restraints on non-hydrogen
    solute atoms using GROMACS.

    Parameters
    ----------
    parametrized_system
        The path to the parametrized system (OpenMM XML or GMX TOP).
    solvated_pdb_file
        The path to the solvated PDB with ions.
    minimized_coords_file
        The path to write the minimized coords (OpenMM PDB or GMX GRO).
    gmx_executable
        Name of GROMACS executable to pass to subprocess.
    nonbonded_cutoff
        The cutoff for the Lennard-Jones potential and PME direct space
        summation.
    vdw_switch_width
        The distance from the nonbonded cutoff at which to apply the
        switching function.
    restraint_energy_constant
        Energy constant for Cartesian restraints (units Energy Length^-2).
    force_tolerance
        Force tolerance for minimization (units Energy Length^-1).
    """

    import subprocess

    # Check units of arguments
    if not nonbonded_cutoff.is_compatible_with(unit.nanometer):
        raise ValueError("nonbonded_cutoff does not have units of Length")
    if not vdw_switch_width.is_compatible_with(unit.nanometer):
        raise ValueError("vdw_switch_width does not have units of Length")
    if not restraint_energy_constant.is_compatible_with(
        unit.kilojoule_per_mole / unit.nanometer**2
    ):
        raise ValueError(
            "restraint_energy_constant does not have units of Energy Length^-2"
        )
    if not force_tolerance.is_compatible_with(unit.kilojoule_per_mole / unit.nanometer):
        raise ValueError("force_tolerance does not have units of Energy Length^-1")

    nonbonded_cutoff = nonbonded_cutoff.m_as(unit.nanometer)
    switch_distance = (nonbonded_cutoff - vdw_switch_width).m_as(unit.nanometer)
    restraint_energy_constant = restraint_energy_constant.m_as(
        unit.kilojoule_per_mole / unit.nanometer
    )
    force_tolerance = force_tolerance.m_as(unit.kilojoule_per_mole / unit.nanometer)

    # Create MDP file
    mdp_file = str(setup_prefix) + "-minimized.mdp"
    with open(mdp_file, "w") as mdp_file_w:
        mdp_file_w.write(
            "define=-DPOSRES\n"
            "integrator = l-bfgs\n"
            f"emtol = {force_tolerance}\n"
            "emstep = 0.01\n"
            "nsteps = -1\n"
            "nstlist = 1\n"
            "constraint_algorithm = lincs\n"
            "lincs-warnangle = 45\n"
            "constraints = h-bonds\n"
            "lincs_iter = 2\n"
            "lincs_order = 4\n"
            "cutoff-scheme = Verlet\n"
            "nstlist = 40\n"
            "vdwtype = cut-off\n"
            "vdw-modifier = potential-switch\n"
            f"rvdw = {self.nonbonded_cutoff}\n"
            f"rvdw-switch = {self.switch_distance}\n"
            "coulombtype = PME\n"
            f"rcoulomb = {self.nonbonded_cutoff}\n"
            "pme_order = 4\n"
            "fourierspacing = 0.16\n"
            "pbc = xyz\n"
            "DispCorr = EnerPres\n"
        )

    # Create position restraints file for backbone atoms
    restr = subprocess.Popen(
        [
            gmx_executable,
            "genrestr",
            "-f",
            solvated_pdb_file,
            "-o",
            f"{setup_prefix}_posre.itp",
            "-fc",
            restraint_energy_constant,
        ],
        stdin=subprocess.PIPE,
    )
    restr.communicate(b"4\n")
    restr.wait()

    # Generate TPR for Energy Minimization
    out_prefix = str(setup_prefix) + "-minimized"
    grompp = subprocess.run(
        [
            gmx_executable,
            "grompp",
            "-f",
            mdp_file,
            "-p",
            parametrized_system,
            "-c",
            solvated_pdb_file,
            "-r",
            solvated_pdb_file,
            "-o",
            out_prefix,
        ]
    )

    # Run Energy Minimization
    run = subprocess.run([gmx_executable, "mdrun", "-deffnm", out_prefix])
