import numpy
from openmm import app, unit
import openmm
from pathlib import Path
from proteinbenchmark.openmm_simulation import OpenMMSimulation
from proteinbenchmark.system_setup import (
    build_initial_coordinates,
    solvate,
    minimize,
)
from proteinbenchmark.simulation_parameters import *
from proteinbenchmark.utilities import exists_and_not_empty


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
        """

        self.target_name = target_name
        self.target_parameters = target_parameters
        self.force_field = force_field_name
        self.water_model = water_model_name
        self.force_field_file = force_field_file
        self.water_model_file = water_model_file

        # Check thermodynamic state
        for quantity in ['pressure', 'temperature', 'ph', 'ionic_strength']:
            if quantity not in self.target_parameters:

                raise ValueError(
                    f'benchmark_targets for target {target_name} must contain '
                    f'"{quantity}"'
                )

        pressure = target_parameters['pressure'].value_in_unit(unit.atmospheres)
        temperature = target_parameters['temperature'].value_in_unit(
            unit.kelvin
        )
        ph = target_parameters['ph']
        ionic_strength = target_parameters['ionic_strength'].value_in_unit(
            unit.molar
        )

        self.system_name = (
            f'{target_name}-{force_field_name}-{water_model_name}-'
            f'{pressure:.2f}-{temperature:.2f}-{ph:.2f}-{ionic_strength:.2f}'
        )

        # Create a directory to store results for this benchmark system
        self.base_path = Path(result_directory, self.system_name)

        # File paths for setup
        self.setup_dir = Path(self.base_path, 'setup')
        self.setup_prefix = Path(self.setup_dir, self.system_name)
        self.minimized_pdb = f'{self.setup_prefix}-minimized.pdb'
        self.openmm_system = f'{self.setup_prefix}-openmm-system.xml'


    def setup(self):
        """
        Build initial coordinates, solvate, and minimize energy. This should be
        deterministic and needs to be run once for all replicas.
        """

        # Create the setup directory if it doesn't already exist
        self.setup_dir.mkdir(parents = True, exist_ok = True)

        initial_pdb = f'{self.setup_prefix}-initial.pdb'
        protonated_pdb = f'{self.setup_prefix}-protonated.pdb'
        solvated_pdb = f'{self.setup_prefix}-solvated.pdb'

        # Build initial coordinates
        if not exists_and_not_empty(protonated_pdb):

            print(f'Building initial coordinates for target {self.target_name}')

            if 'initial_pdb' in self.target_parameters:

                # Copy initial PDB to results directory
                Path(initial_pdb).write_text(
                    self.target_parameters['initial_pdb'].read_text()
                )

                build_initial_coordinates(
                    build_method = 'pdb',
                    ph = self.target_parameters['ph'],
                    initial_pdb = initial_pdb,
                    protonated_pdb = protonated_pdb,
                )

            elif 'aa_sequence' in self.target_parameters:

                build_initial_coordinates(
                    build_method = 'extended',
                    ph = self.target_parameters['ph'],
                    initial_pdb = initial_pdb,
                    protonated_pdb = protonated_pdb,
                    aa_sequence = self.target_parameters['aa_sequence'],
                )

            else:

                raise ValueError(
                    f'benchmark_targets for target {self.target_name} must '
                    'contain one of "aa_sequence" or "initial_pdb"'
                )

        # Solvate, add ions, and construct OpenMM system
        if not exists_and_not_empty(self.openmm_system):

            print(f'Solvating target {self.target_name}')

            # Get parameters for solvation and constructing OpenMM system
            if 'solvent_padding' in self.target_parameters:
                solvent_padding = self.target_parameters['solvent_padding']
            elif self.target_parameters['target_type'] == 'disordered':
                solvent_padding = DISORDERED_SOLVENT_PADDING
            else:
                solvent_padding = SOLVENT_PADDING

            if 'nonbonded_cutoff' in self.target_parameters:
                nonbonded_cutoff = self.target_parameters['nonbonded_cutoff']
            else:
                nonbonded_cutoff = NONBONDED_CUTOFF

            if 'vdw_switch_width' in self.target_parameters:
                vdw_switch_width = self.target_parameters['vdw_switch_width']
            else:
                vdw_switch_width = VDW_SWITCH_WIDTH

            solvate(
                solvent_padding = solvent_padding,
                ionic_strength = self.target_parameters['ionic_strength'],
                nonbonded_cutoff = nonbonded_cutoff,
                vdw_switch_width = vdw_switch_width,
                protonated_pdb_file = protonated_pdb,
                solvated_pdb_file = solvated_pdb,
                openmm_system_xml = self.openmm_system,
                water_model = self.water_model,
                force_field_file = self.force_field_file,
                water_model_file = self.water_model_file,
            )

        # Minimize energy of solvated system with Cartesian restraints on
        # non-hydrogen solute atoms
        if not exists_and_not_empty(self.minimized_pdb):

            print(f'Minimizing energy for target {self.target_name}')

            if 'restraint_energy_constant' in self.target_parameters:

                restraint_energy_constant = (
                    self.target_parameters['restraint_energy_constant']
                )

            else:
                restraint_energy_constant = RESTRAINT_ENERGY_CONSTANT

            minimize(
                restraint_energy_constant = restraint_energy_constant,
                openmm_system_xml = self.openmm_system,
                solvated_pdb_file = solvated_pdb,
                minimized_pdb_file = self.minimized_pdb,
            )

        print(f'Setup complete for target {self.target_name}')


    def run_simulations(self, replica: int = 1):
        """Equilibrate and run production trajectories for one replica."""

        # Create a directory for this replica if it doesn't already exist
        replica_dir = Path(self.base_path, f'replica-{replica:d}')
        replica_dir.mkdir(parents = True, exist_ok = True)

        replica_prefix = Path(replica_dir, self.system_name)
        equil_prefix = f'{replica_prefix}-equilibration'
        prod_prefix = f'{replica_prefix}-production'

        # Serialized OpenMM state from the end of the equilibration simulation
        equilibrated_state = f'{equil_prefix}-1.xml'

        # Equilibrate at constant pressure and temperature
        if not exists_and_not_empty(equilibrated_state):

            print(f'Running NPT equilibration for target {self.target_name}')

            # Get parameters for equilibration simulation
            if 'equil_langevin_friction' in self.target_parameters:

                equil_langevin_friction = (
                    self.target_parameters['equil_langevin_friction']
                )

            else:
                equil_langevin_friction = EQUIL_LANGEVIN_FRICTION

            if 'equil_barostat_frequency' in self.target_parameters:

                equil_barostat_frequency = (
                    self.target_parameters['equil_barostat_frequency']
                )

            else:
                equil_barostat_frequency = EQUIL_BAROSTAT_FREQUENCY

            if 'equil_timestep' in self.target_parameters:
                equil_timestep = self.target_parameters['equil_timestep']
            else:
                equil_timestep = EQUIL_TIMESTEP

            if 'equil_traj_length' in self.target_parameters:
                equil_traj_length = self.target_parameters['equil_traj_length']
            else:
                equil_traj_length = EQUIL_TRAJ_LENGTH

            if 'equil_frame_length' in self.target_parameters:

                equil_frame_length = (
                    self.target_parameters['equil_frame_length']
                )

            else:
                equil_frame_length = EQUIL_FRAME_LENGTH

            # Initialize the equilibration simulation
            equilibration_dcd = f'{equil_prefix}.dcd'
            equilibration_state_data = f'{equil_prefix}.out'
            equilibration_checkpoint = f'{equil_prefix}.chk'

            equilibration_simulation = OpenMMSimulation(
                openmm_system_file = self.openmm_system,
                initial_pdb_file = self.minimized_pdb,
                dcd_reporter_file = equilibration_dcd,
                state_reporter_file = equilibration_state_data,
                checkpoint_file = equilibration_checkpoint,
                save_state_prefix = equil_prefix,
                temperature = self.target_parameters['temperature'],
                pressure = self.target_parameters['pressure'],
                langevin_friction = equil_langevin_friction,
                barostat_frequency = equil_barostat_frequency,
                timestep = equil_timestep,
                traj_length = equil_traj_length,
                frame_length = equil_frame_length,
                checkpoint_length = equil_traj_length,
                save_state_length = equil_traj_length,
            )

            # Run equilibration
            equilibration_simulation.start_from_pdb()

        print(f'Running NPT production for target {self.target_name}')

        # Get parameters for production simulation
        if 'langevin_friction' in self.target_parameters:
            langevin_friction = self.target_parameters['langevin_friction']
        else:
            langevin_friction = LANGEVIN_FRICTION

        if 'barostat_frequency' in self.target_parameters:
            barostat_frequency = (self.target_parameters['barostat_frequency'])
        else:
            barostat_frequency = BAROSTAT_FREQUENCY

        if 'timestep' in self.target_parameters:
            timestep = self.target_parameters['timestep']
        else:
            timestep = TIMESTEP

        if 'traj_length' in self.target_parameters:
            traj_length = self.target_parameters['traj_length']
        elif self.target_parameters['target_type'] == 'peptide':
            traj_length = PEPTIDE_TRAJ_LENGTH
        elif self.target_parameters['target_type'] == 'folded':
            traj_length = FOLDED_TRAJ_LENGTH
        elif self.target_parameters['target_type'] == 'disordered':
            traj_length = DISORDERED_TRAJ_LENGTH

        else:

            raise ValueError(
                f'benchmark_targets for target {self.target_name} must '
                f'contain "traj_length" or "target_type" must be one of '
                '"peptide", "folded", or "disordered".'
            )

        if 'frame_length' in self.target_parameters:
            frame_length = self.target_parameters['frame_length']
        else:
            frame_length = FRAME_LENGTH

        if 'checkpoint_length' in self.target_parameters:
            checkpoint_length = self.target_parameters['checkpoint_length']
        else:
            checkpoint_length = CHECKPOINT_LENGTH

        if 'save_state_length' in self.target_parameters:
            save_state_length = self.target_parameters['save_state_length']
        else:
            save_state_length = SAVE_STATE_LENGTH

        # Initialize the production simulation
        production_dcd = f'{prod_prefix}.dcd'
        production_state_data = f'{prod_prefix}.out'
        production_checkpoint = f'{prod_prefix}.chk'

        production_simulation = OpenMMSimulation(
            openmm_system_file = self.openmm_system,
            initial_pdb_file = self.minimized_pdb,
            dcd_reporter_file = production_dcd,
            state_reporter_file = production_state_data,
            checkpoint_file = production_checkpoint,
            save_state_prefix = prod_prefix,
            temperature = self.target_parameters['temperature'],
            pressure = self.target_parameters['pressure'],
            langevin_friction = langevin_friction,
            barostat_frequency = barostat_frequency,
            timestep = timestep,
            traj_length = traj_length,
            frame_length = frame_length,
            checkpoint_length = checkpoint_length,
            save_state_length = save_state_length,
        )

        # Run production
        if not exists_and_not_empty(production_checkpoint):

            # Start production simulation, initializing positions and velocities
            # to the final state from the equilibration simulation
            production_simulation.start_from_save_state(equilibrated_state)

        else:

            # Resume from a previous production checkpoint
            production_simulation.resume_from_checkpoint()


    def analyze(self):
        """Process trajectories and estimate observables."""

        pass

