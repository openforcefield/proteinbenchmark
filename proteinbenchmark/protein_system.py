from pathlib import Path

from proteinbenchmark.analysis import (
    align_trajectory,
    assign_dihedral_clusters,
    compute_fraction_helix,
    compute_h_bond_scalar_couplings,
    compute_scalar_couplings,
    measure_dihedrals,
    measure_h_bond_geometries,
)
from proteinbenchmark.gmx_simulation import GMXSimulation
from proteinbenchmark.openmm_simulation import OpenMMSimulation
from proteinbenchmark.simulation_parameters import *
from proteinbenchmark.system_setup import (
    assign_parameters,
    build_initial_coordinates,
    minimize_gmx,
    minimize_openmm,
    solvate,
)
from proteinbenchmark.utilities import exists_and_not_empty, merge_csvs


class ProteinBenchmarkSystem:
    """
    A class representing a benchmark system with a force field, water model, and
    thermodynamic state (pressure, temperature, pH, and ionic strength).
    """

    def __init__(
        self,
        result_directory: str,
        target_name: str,
        target_parameters: dict,
        force_field_name: str,
        water_model_name: str,
        force_field_file: str,
        water_model_file: str = None,
        simulation_platform: str = "openmm",
        gmx_executable: str = None,
    ):
        """
        Initializes the ProteinBenchmarkSystem object with target parameters.

        Parameters
        ----------
        result_directory
            The path to the top level directory where results will be stored.
        target_name
            The name of this benchmark target.
        target_parameters
            A dictionary of parameters, including thermodynamic state
            (temperature, pressure, pH, and ionic strength) and build method
            (initial PDB or amino acid sequence).
        force_field_name
            The name of the force field for this benchmark.
        water_model_name
            The water model for this benchmark.
        force_field_file
            The name of the file containing the force field parameters.
        water_model_file
            The name of the file containing the water model parameters.
        simuation_platform
            Simulation platform from which to run energy minimization,
            equilibration, and production simulations.
        gmx_executable
            Name of GROMACS executable to pass to subprocess.
        """

        self.target_name = target_name
        self.target_parameters = target_parameters
        self.force_field = force_field_name
        self.water_model = water_model_name
        self.force_field_file = force_field_file
        self.water_model_file = water_model_file
        self.simulation_platform = simulation_platform.lower()
        self.gmx_executable = gmx_executable

        # Check thermodynamic state
        for quantity in ["pressure", "temperature", "ph", "ionic_strength"]:
            if quantity not in self.target_parameters:
                raise ValueError(
                    f"benchmark_targets for target {target_name} must contain "
                    f'"{quantity}"'
                )

        # Set up nonbonded parameters. Do this here so these values are
        # persistent.
        if "nonbonded_cutoff" in self.target_parameters:
            self.nonbonded_cutoff = self.target_parameters["nonbonded_cutoff"]
        else:
            self.nonbonded_cutoff = NONBONDED_CUTOFF

        if "vdw_switch_width" in self.target_parameters:
            self.vdw_switch_width = self.target_parameters["vdw_switch_width"]
        else:
            self.vdw_switch_width = VDW_SWITCH_WIDTH

        # Create a directory to store results for this benchmark system
        self.system_name = f"{target_name}-{force_field_name}"
        self.base_path = Path(result_directory, self.system_name)

        # File paths for setup
        self.setup_dir = Path(self.base_path, "setup")
        self.setup_prefix = Path(self.setup_dir, self.system_name)
        self.initial_pdb = f"{self.setup_prefix}-initial.pdb"
        self.protonated_pdb = f"{self.setup_prefix}-protonated.pdb"

        if self.simulation_platform == "openmm":
            self.minimized_coords = f"{self.setup_prefix}-minimized.pdb"
            self.parametrized_system = f"{self.setup_prefix}-openmm-system.xml"
        elif self.simulation_platform == "gmx":
            self.minimized_coords = f"{self.setup_prefix}-minimized.gro"
            self.parametrized_system = f"{self.setup_prefix}-gmx-system.top"
        else:
            raise ValueError(
                f"simulation_platform must be one of:\n    openmm\n    gmx"
            )

    def setup(self):
        """
        Build initial coordinates, solvate, and minimize energy. This should be
        deterministic and needs to be run once for all replicas.
        """

        # Create the setup directory if it doesn't already exist
        self.setup_dir.mkdir(parents=True, exist_ok=True)

        solvated_pdb = f"{self.setup_prefix}-solvated.pdb"

        # Build initial coordinates
        if not exists_and_not_empty(self.protonated_pdb):
            print(f"Building initial coordinates for system {self.system_name}")

            if "initial_pdb" in self.target_parameters:
                # Copy initial PDB to results directory
                Path(self.initial_pdb).write_text(
                    self.target_parameters["initial_pdb"].read_text()
                )

                build_initial_coordinates(
                    build_method="pdb",
                    ph=self.target_parameters["ph"],
                    initial_pdb=self.initial_pdb,
                    protonated_pdb=self.protonated_pdb,
                )

            elif "aa_sequence" in self.target_parameters:
                if "build_method" in self.target_parameters:
                    build_method = self.target_parameters["build_method"]
                else:
                    build_method = "extended"

                if "nterm_cap" in self.target_parameters:
                    nterm_cap = self.target_parameters["nterm_cap"]
                else:
                    nterm_cap = None

                if "cterm_cap" in self.target_parameters:
                    cterm_cap = self.target_parameters["cterm_cap"]
                else:
                    cterm_cap = None

                build_initial_coordinates(
                    build_method=build_method,
                    ph=self.target_parameters["ph"],
                    initial_pdb=self.initial_pdb,
                    protonated_pdb=self.protonated_pdb,
                    aa_sequence=self.target_parameters["aa_sequence"],
                    nterm_cap=nterm_cap,
                    cterm_cap=cterm_cap,
                )

            else:
                raise ValueError(
                    f"benchmark_targets for target {self.target_name} must "
                    'contain one of "aa_sequence" or "initial_pdb"'
                )

        # Solvate and add ions using Amber ff14SB as a reference force field
        if not exists_and_not_empty(solvated_pdb):
            print(f"Solvating system {self.system_name}")

            # Get parameters for solvation and constructing OpenMM system
            if "solvent_padding" in self.target_parameters:
                solvent_padding = self.target_parameters["solvent_padding"]
            elif self.target_parameters["target_type"] != "folded":
                solvent_padding = DISORDERED_SOLVENT_PADDING
            else:
                solvent_padding = SOLVENT_PADDING

            solvate(
                ionic_strength=self.target_parameters["ionic_strength"],
                protonated_pdb_file=self.protonated_pdb,
                solvated_pdb_file=solvated_pdb,
                water_model=self.water_model,
                solvent_padding=solvent_padding,
            )

        # Construct parametrized system
        if not exists_and_not_empty(self.parametrized_system):
            print(f"Parametrizing system {self.system_name}")

            assign_parameters(
                simulation_platform=self.simulation_platform,
                nonbonded_cutoff=self.nonbonded_cutoff,
                vdw_switch_width=self.vdw_switch_width,
                protonated_pdb_file=self.protonated_pdb,
                solvated_pdb_file=solvated_pdb,
                parametrized_system=self.parametrized_system,
                water_model=self.water_model,
                force_field_file=self.force_field_file,
                water_model_file=self.water_model_file,
                setup_prefix=self.setup_prefix,
            )

        # Minimize energy of solvated system with Cartesian restraints on
        # non-hydrogen solute atoms
        if not exists_and_not_empty(self.minimized_coords):
            print(f"Minimizing energy for system {self.system_name}")

            if "restraint_energy_constant" in self.target_parameters:
                restraint_energy_constant = self.target_parameters[
                    "restraint_energy_constant"
                ]
            else:
                restraint_energy_constant = RESTRAINT_ENERGY_CONSTANT

            if "force_tolerance" in self.target_parameters:
                force_tolerance = self.target_parameters["force_tolerance"]
            else:
                force_tolerance = FORCE_TOLERANCE

            if self.simulation_platform == "openmm":
                minimize_openmm(
                    parametrized_system=self.parametrized_system,
                    solvated_pdb_file=solvated_pdb,
                    minimized_coords_file=self.minimized_coords,
                    restraint_energy_constant=restraint_energy_constant,
                    force_tolerance=force_tolerance,
                )

            elif self.simulation_platform == "gmx":
                minimize_gmx(
                    parametrized_system=self.parametrized_system,
                    solvated_pdb_file=solvated_pdb,
                    minimized_coords_file=self.minimized_coords,
                    setup_prefix=self.setup_prefix,
                    gmx_executable=self.gmx_executable,
                    nonbonded_cutoff=self.nonbonded_cutoff,
                    vdw_switch_width=self.vdw_switch_width,
                    restraint_energy_constant=restraint_energy_constant,
                    force_tolerance=force_tolerance,
                )

        print(f"Setup complete for system {self.system_name}")

    def run_simulations(self, replica: int = 1):
        """Equilibrate and run production trajectories for one replica."""

        # Create a directory for this replica if it doesn't already exist
        replica_dir = Path(self.base_path, f"replica-{replica:d}")
        replica_dir.mkdir(parents=True, exist_ok=True)

        replica_prefix = Path(replica_dir, self.system_name)
        setup_dir = Path(self.base_path, "setup")
        setup_prefix = Path(setup_dir, self.system_name)

        equil_prefix = f"{replica_prefix}-equilibration"
        prod_prefix = f"{replica_prefix}-production"

        # Saved state from the end of the equilibration simulation
        if self.simulation_platform == "openmm":
            equilibrated_state = f"{equil_prefix}-1.xml"
        elif self.simulation_platform == "gmx":
            equilibrated_state = f"{equil_prefix}.cpt"

        # Equilibrate at constant pressure and temperature
        if not exists_and_not_empty(equilibrated_state):
            print(f"Running NPT equilibration for system {self.system_name}")

            # Get parameters for equilibration simulation
            if "equil_timestep" in self.target_parameters:
                equil_timestep = self.target_parameters["equil_timestep"]
            else:
                equil_timestep = EQUIL_TIMESTEP

            if "equil_traj_length" in self.target_parameters:
                equil_traj_length = self.target_parameters["equil_traj_length"]
            else:
                equil_traj_length = EQUIL_TRAJ_LENGTH

            if "equil_frame_length" in self.target_parameters:
                equil_frame_length = self.target_parameters["equil_frame_length"]
            else:
                equil_frame_length = EQUIL_FRAME_LENGTH

            if "equil_langevin_friction" in self.target_parameters:
                equil_langevin_friction = self.target_parameters[
                    "equil_langevin_friction"
                ]
            else:
                equil_langevin_friction = EQUIL_LANGEVIN_FRICTION

            # Initialize the equilibration simulation
            if self.simulation_platform == "openmm":
                # Get OpenMM Specific Parameters
                if "equil_barostat_frequency" in self.target_parameters:
                    equil_barostat_frequency = self.target_parameters[
                        "equil_barostat_frequency"
                    ]
                else:
                    equil_barostat_frequency = EQUIL_OPENMM_BAROSTAT_FREQUENCY

                equilibration_dcd = f"{equil_prefix}.dcd"
                equilibration_state_data = f"{equil_prefix}.out"
                equilibration_checkpoint = f"{equil_prefix}.chk"

                equilibration_simulation = OpenMMSimulation(
                    openmm_system_file=self.parametrized_system,
                    initial_pdb_file=self.minimized_coords,
                    dcd_reporter_file=equilibration_dcd,
                    state_reporter_file=equilibration_state_data,
                    checkpoint_file=equilibration_checkpoint,
                    save_state_prefix=equil_prefix,
                    temperature=self.target_parameters["temperature"].to_openmm(),
                    pressure=self.target_parameters["pressure"].to_openmm(),
                    langevin_friction=equil_langevin_friction.to_openmm(),
                    barostat_frequency=equil_barostat_frequency,
                    timestep=equil_timestep.to_openmm(),
                    traj_length=equil_traj_length.to_openmm(),
                    frame_length=equil_frame_length.to_openmm(),
                    checkpoint_length=equil_traj_length.to_openmm(),
                    save_state_length=equil_traj_length.to_openmm(),
                )

            else:
                # Get GROMACS Specific Parameters
                if "equil_barostat_constant" in self.target_parameters:
                    equil_barostat_constant = self.target_parameters[
                        "equil_barostat_constant"
                    ]
                else:
                    equil_barostat_constant = EQUIL_GMX_BAROSTAT_CONSTANT

                equilibration_simulation = GMXSimulation(
                    gmx_executable=self.gmx_executable,
                    parametrized_system=self.parametrized_system,
                    initial_coords_file=self.minimized_coords,
                    output_prefix=equil_prefix,
                    nonbonded_cutoff=self.nonbonded_cutoff,
                    vdw_switch_width=self.vdw_switch_width,
                    temperature=self.target_parameters["temperature"],
                    pressure=self.target_parameters["pressure"],
                    langevin_friction=equil_langevin_friction,
                    barostat_constant=equil_barostat_constant,
                    timestep=equil_timestep,
                    traj_length=equil_traj_length,
                    frame_length=equil_frame_length,
                )

            # Run equilibration
            equilibration_simulation.start_from_pdb()

        print(f"Running NPT production for system {self.system_name}")

        # Get parameters for production simulation
        if "timestep" in self.target_parameters:
            timestep = self.target_parameters["timestep"]
        else:
            timestep = TIMESTEP

        if "traj_length" in self.target_parameters:
            traj_length = self.target_parameters["traj_length"]
        elif self.target_parameters["target_type"] == "peptide":
            traj_length = PEPTIDE_TRAJ_LENGTH
        elif self.target_parameters["target_type"] == "folded":
            traj_length = FOLDED_TRAJ_LENGTH
        elif self.target_parameters["target_type"] == "disordered":
            traj_length = DISORDERED_TRAJ_LENGTH
        else:
            raise ValueError(
                f"benchmark_targets for target {self.target_name} must "
                f'contain "traj_length" or "target_type" must be one of '
                '"peptide", "folded", or "disordered".'
            )

        if "frame_length" in self.target_parameters:
            frame_length = self.target_parameters["frame_length"]
        else:
            frame_length = FRAME_LENGTH

        if "langevin_friction" in self.target_parameters:
            langevin_friction = self.target_parameters["langevin_friction"]
        else:
            langevin_friction = LANGEVIN_FRICTION

        # Initialize the production simulation
        if self.simulation_platform == "openmm":
            # Get OpenMM Specific Parameters
            if "checkpoint_length" in self.target_parameters:
                checkpoint_length = self.target_parameters["checkpoint_length"]
            else:
                checkpoint_length = CHECKPOINT_LENGTH

            if "save_state_length" in self.target_parameters:
                save_state_length = self.target_parameters["save_state_length"]
            else:
                save_state_length = SAVE_STATE_LENGTH

            if "barostat_frequency" in self.target_parameters:
                barostat_frequency = self.target_parameters["barostat_frequency"]
            else:
                barostat_frequency = OPENMM_BAROSTAT_FREQUENCY

            production_dcd = f"{prod_prefix}.dcd"
            production_state_data = f"{prod_prefix}.out"
            production_checkpoint = f"{prod_prefix}.chk"

            production_simulation = OpenMMSimulation(
                openmm_system_file=self.parametrized_system,
                initial_pdb_file=self.minimized_coords,
                dcd_reporter_file=production_dcd,
                state_reporter_file=production_state_data,
                checkpoint_file=production_checkpoint,
                save_state_prefix=prod_prefix,
                temperature=self.target_parameters["temperature"].to_openmm(),
                pressure=self.target_parameters["pressure"].to_openmm(),
                langevin_friction=langevin_friction.to_openmm(),
                barostat_frequency=barostat_frequency,
                timestep=timestep.to_openmm(),
                traj_length=traj_length.to_openmm(),
                frame_length=frame_length.to_openmm(),
                checkpoint_length=checkpoint_length.to_openmm(),
                save_state_length=save_state_length.to_openmm(),
            )

        else:
            # Get GROMACS Specific Parameters
            if "barostat_constant" in self.target_parameters:
                barostat_constant = self.target_parameters["barostat_constant"]
            else:
                barostat_constant = GMX_BAROSTAT_CONSTANT

            production_checkpoint = f"{prod_prefix}.cpt"

            production_simulation = GMXSimulation(
                gmx_executable=self.gmx_executable,
                parametrized_system=self.parametrized_system,
                initial_coords_file=equil_prefix,
                output_prefix=prod_prefix,
                nonbonded_cutoff=self.nonbonded_cutoff,
                vdw_switch_width=self.vdw_switch_width,
                temperature=self.target_parameters["temperature"],
                pressure=self.target_parameters["pressure"],
                langevin_friction=langevin_friction,
                barostat_constant=barostat_constant,
                timestep=timestep,
                traj_length=traj_length,
                frame_length=frame_length,
            )

        # Run production
        if not exists_and_not_empty(production_checkpoint):
            # Start production simulation, initializing positions and
            # velocities to the final state from the equilibration simulation
            production_simulation.start_from_save_state(equilibrated_state)

        else:
            # Resume from a previous production checkpoint
            production_simulation.resume_from_checkpoint()

    def analyze_observables(self, replica: int = 1):
        """Process trajectories and estimate observables."""

        analysis_dir = Path(self.base_path, "analysis")
        analysis_dir.mkdir(parents=True, exist_ok=True)

        analysis_prefix = Path(analysis_dir, f"{self.system_name}-{replica}")

        reimaged_topology = f"{analysis_prefix}-reimaged.pdb"
        reimaged_trajectory = f"{analysis_prefix}-reimaged.dcd"

        if "frame_length" in self.target_parameters:
            frame_length = self.target_parameters["frame_length"]
        else:
            frame_length = FRAME_LENGTH

        # Align production trajectory
        if not exists_and_not_empty(reimaged_topology):
            print(f"Aligning production trajectory for system {self.system_name}")

            replica_dir = Path(self.base_path, f"replica-{replica:d}")
            replica_prefix = Path(replica_dir, self.system_name)

            if self.simulation_platform == "openmm":
                traj_path = f"{replica_prefix}-production.dcd"
                output_selection = 'chainid == "A"'
            elif self.simulation_platform == "gmx":
                traj_path = f"{replica_dir}/traj.xtc"
                output_selection = 'resname != "HOH" && resname != "NA"'

            align_trajectory(
                topology_path=self.minimized_coords,
                trajectory_path=traj_path,
                output_prefix=f"{analysis_prefix}-reimaged",
                output_selection=output_selection,
                align_selection='name == "CA"',
                reference_path=self.initial_pdb,
            )

        # Measure dihedrals
        dihedrals = f"{analysis_prefix}-dihedrals.dat"

        if not exists_and_not_empty(dihedrals):
            print(f"Measuring dihedrals for system {self.system_name}")

            fragment_index = measure_dihedrals(
                topology_path=reimaged_topology,
                trajectory_path=reimaged_trajectory,
                frame_length=frame_length,
                output_path=dihedrals,
            )

            if fragment_index > 0:
                merge_csvs(dihedrals)

        # Measure hydrogen bond geometries
        h_bond_geometries = f"{analysis_prefix}-hydrogen-bond-geometries.dat"

        if not exists_and_not_empty(h_bond_geometries):
            print(
                "Measuring hydrogen bond geometries for system " f"{self.system_name}"
            )

            fragment_index = measure_h_bond_geometries(
                topology_path=self.protonated_pdb,
                trajectory_path=reimaged_trajectory,
                frame_length=frame_length,
                output_path=h_bond_geometries,
            )

            if fragment_index > 0:
                merge_csvs(h_bond_geometries)

        # Dihedral cluster assignments
        dihedral_clusters = f"{analysis_prefix}-dihedral-clusters.dat"

        if not exists_and_not_empty(dihedral_clusters):
            print(f"Assigning dihedral clusters for system {self.system_name}")

            assign_dihedral_clusters(
                dihedrals_path=dihedrals,
                output_path=dihedral_clusters,
            )

        # Compute observables
        target_observables = self.target_parameters["observables"]

        # Scalar couplings
        scalar_couplings = f"{analysis_prefix}-scalar-couplings.dat"

        if "scalar_couplings" in target_observables and not exists_and_not_empty(
            scalar_couplings
        ):
            print(f"Computing scalar couplings for system {self.system_name}")

            data = target_observables["scalar_couplings"]["observable_path"]

            compute_scalar_couplings(
                observable_path=data,
                dihedrals_path=dihedrals,
                output_path=scalar_couplings,
            )

        # Hydrogen bond interresidue scalar couplings
        h_bond_scalar_couplings = f"{analysis_prefix}-h-bond-scalar-couplings.dat"

        if (
            "h_bond_scalar_couplings" in target_observables
            and not exists_and_not_empty(h_bond_scalar_couplings)
        ):
            print(
                "Computing hydrogen bond scalar couplings for system "
                f"{self.system_name}"
            )

            data = target_observables["h_bond_scalar_couplings"]["observable_path"]

            compute_h_bond_scalar_couplings(
                observable_path=data,
                h_bond_geometries_path=h_bond_geometries,
                output_path=h_bond_scalar_couplings,
            )

        # Fraction helix
        fraction_helix = f"{analysis_prefix}-fraction-helix.dat"

        if "fraction_helix" in target_observables and not exists_and_not_empty(
            fraction_helix
        ):
            print(f"Computing fraction helix for system {self.system_name}")

            data = target_observables["fraction_helix"]["observable_path"]

            compute_fraction_helix(
                observable_path=data,
                dihedral_clusters_path=dihedral_clusters,
                h_bond_geometries_path=h_bond_geometries,
                output_path=fraction_helix,
            )
