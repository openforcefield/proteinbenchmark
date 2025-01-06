from pathlib import Path

import numpy
import openmm
from openmm import app, unit

from proteinbenchmark.utilities import exists_and_not_empty, read_xml


class OpenMMSimulation:
    """A class representing a simulation in OpenMM."""

    def __init__(
        self,
        openmm_system_file: str,
        initial_pdb_file: str,
        dcd_reporter_file: str,
        state_reporter_file: str,
        checkpoint_file: str,
        save_state_prefix: str,
        temperature: unit.Quantity,
        pressure: unit.Quantity,
        langevin_friction: unit.Quantity,
        barostat_frequency: int,
        timestep: unit.Quantity,
        traj_length: unit.Quantity,
        frame_length: unit.Quantity,
        checkpoint_length: unit.Quantity,
        save_state_length: unit.Quantity,
    ):
        """
        Initializes the simulation parameters and checks units.

        Parameters
        ----------
        openmm_system_file
            The path to the parametrized OpenMM system as a serialized XML.
        initial_pdb_file
            Path to PDB file used to set initial coordinates.
        dcd_reporter_file
            Path to DCD file to write trajectory coordinates.
        state_reporter_file
            Path to file to write state data, e.g. energies and temperature.
        checkpoint_file
            Path to file to write binary checkpoints.
        save_state_prefix
            Path prefix for files to write serialized simulation states.
        temperature
            The target temperature of the Langevin thermostat.
        pressure
            The target pressure of the Monte Carlo barostat.
        langevin_friction
            The collision frequency of the Langevin integrator.
        barostat_frequency
            The number of steps between attempted pressure changes for the Monte
            Carlo barostat.
        timestep
            The timestep for the Langevin integrator.
        traj_length
            The length of the trajectory (N_steps = traj_length / timestep).
        frame_length
            The amount of time between writing coordinates and state data to
            disk.
        checkpoint_length
            The amount of time between writing binary checkpoints to disk.
        save_state_length
            The amount of time between writing serialized simulation states to
            disk.
        """

        self.openmm_system_file = openmm_system_file
        self.initial_pdb_file = initial_pdb_file
        self.dcd_reporter_file = dcd_reporter_file
        self.state_reporter_file = state_reporter_file
        self.checkpoint_file = checkpoint_file
        self.save_state_prefix = save_state_prefix

        # Check units of arguments
        if not temperature.unit.is_compatible(unit.kelvin):
            raise ValueError("temperature does not have units of Temperature")
        if not pressure.unit.is_compatible(unit.atmosphere):
            raise ValueError("pressure does not have units of Mass Length^-1 Time^-2")
        if not langevin_friction.unit.is_compatible(unit.picosecond**-1):
            raise ValueError("langevin_friction does not have units of Time^-1")
        if not timestep.unit.is_compatible(unit.picosecond):
            raise ValueError("timestep does not have units of Time")
        if not traj_length.unit.is_compatible(unit.picosecond):
            raise ValueError("traj_length does not have units of Time")
        if not frame_length.unit.is_compatible(unit.picosecond):
            raise ValueError("frame_length does not have units of Time")
        if not checkpoint_length.unit.is_compatible(unit.picosecond):
            raise ValueError("checkpoint_length does not have units of Time")
        if not save_state_length.unit.is_compatible(unit.picosecond):
            raise ValueError("save_state_length does not have units of Time")

        self.temperature = temperature
        self.pressure = pressure
        self.langevin_friction = langevin_friction
        self.barostat_frequency = barostat_frequency
        self.timestep = timestep
        self.n_steps = int(numpy.round(traj_length / timestep))
        self.output_frequency = int(numpy.round(frame_length / timestep))
        self.checkpoint_frequency = int(numpy.round(checkpoint_length / timestep))
        self.save_state_frequency = int(numpy.round(save_state_length / timestep))

        print(
            "Running simulation with\n    temperature "
            f"{temperature.value_in_unit(unit.kelvin)} K"
            f"\n    pressure {pressure.value_in_unit(unit.atmosphere):.3f} atm"
            f"\n    langevin_friction "
            f"{langevin_friction.value_in_unit(unit.picosecond**-1):.1f} ps^-1"
            f"\n    barostat_frequency {barostat_frequency:d} steps"
            f"\n    timestep {timestep.value_in_unit(unit.femtosecond):.1f} fs"
            f"\n    n_steps {self.n_steps:d} steps"
            f"\n    output_frequency {self.output_frequency:d} steps"
            f"\n    checkpoint_frequency {self.checkpoint_frequency:d} steps"
            f"\n    save_state_frequency {self.save_state_frequency:d} steps"
        )

    def setup_simulation(
        self,
        return_pdb: bool = False,
    ):
        """
        Set up an OpenMM simulation with a Langevin integrator and a Monte Carlo
        barostat.

        Parameters
        ----------
        return_pdb
            Return OpenMM PDBFile as well as OpenMM Simulation.
        """

        # Load OpenMM system and initial PDB
        openmm_system = read_xml(self.openmm_system_file)
        initial_pdb = app.PDBFile(self.initial_pdb_file)

        # Set up BAOAB Langevin integrator from openmmtools with VRORV splitting
        integrator = openmm.LangevinMiddleIntegrator(
            self.temperature,
            self.langevin_friction,
            self.timestep,
        )

        # Set up Monte Carlo barostat
        if self.pressure.value_in_unit(unit.atmosphere) > 0:
            openmm_system.addForce(
                openmm.MonteCarloBarostat(
                    self.pressure,
                    self.temperature,
                    self.barostat_frequency,
                )
            )

        # Create simulation
        # only use cuda if available
        try:
            platform = openmm.Platform.getPlatformByName("CUDA")
            platform_properties = {"Precision": "mixed"}
        except Exception:
            platform = openmm.Platform.getPlatformByName("Reference")
            platform_properties = None
        simulation = app.Simulation(
            initial_pdb.topology,
            openmm_system,
            integrator,
            platform,
            platform_properties
        )

        if return_pdb:
            return simulation, initial_pdb
        else:
            return simulation

    def start_from_pdb(self):
        """
        Start a new simulation initializing positions from a PDB and velocities
        to random samples from a Boltzmann distribution at the simulation
        temperature.
        """

        # Create an OpenMM simulation
        simulation, initial_pdb = self.setup_simulation(return_pdb=True)

        # Initialize positions from the topology PDB
        simulation.context.setPositions(initial_pdb.positions)

        # Initialize velocities to random samples from a Boltzmann distribution
        simulation.context.setVelocitiesToTemperature(self.temperature)

        # Run dynamics
        self.run_dynamics(simulation, append=False)

    def start_from_save_state(
        self,
        save_state_file: str,
    ):
        """
        Start a new simulation initializing positions and velocities from a
        serialized OpenMM state.

        Parameters
        ----------
        save_state_file
            Path to the serialized simulation state.
        """

        # Create an OpenMM simulation
        simulation = self.setup_simulation()

        # Load the serialized state
        if not exists_and_not_empty(save_state_file):
            raise ValueError(
                f"Serialized OpenMM state file {save_state_file} does not "
                "exist or is empty."
            )

        simulation.loadState(save_state_file)

        # Set step number and simulation time to zero
        simulation.context.setStepCount(0)
        simulation.context.setTime(0.0 * unit.picosecond)

        # Run dynamics
        self.run_dynamics(simulation, append=False)

    def resume_from_checkpoint(self):
        """
        Resume an existing OpenMM simulation from a binary checkpoint. The state
        data and DCD reporter files will be truncated to the expected number of
        frames from the checkpoint.
        """

        import mdtraj
        from mdtraj.formats.dcd import DCDTrajectoryFile
        from mdtraj.utils import in_units_of

        # Create an OpenMM simulation
        simulation = self.setup_simulation()

        # Load the checkpoint
        if not exists_and_not_empty(self.checkpoint_file):
            raise ValueError(
                f"Checkpoint file {self.checkpoint_file} does not exist or is " "empty."
            )

        simulation.loadCheckpoint(self.checkpoint_file)

        # Check whether the simulation has already finished
        if simulation.currentStep == self.n_steps:
            return

        # Get expected number of frames based on current step from checkpoint
        # If the state data or DCD reporters have additional frames written,
        # truncate them to the expected number of frames
        expected_frame_count = int(simulation.currentStep / self.output_frequency)

        # Check number of frames in state data reporter file
        with open(self.state_reporter_file, "r") as state_reporter:
            # Subtract one for header line
            state_reporter_frames = sum(1 for _ in state_reporter) - 1

        if state_reporter_frames < expected_frame_count:
            raise ValueError(
                f"The state data reporter file has {state_reporter_frames:d} "
                f"frames but {expected_number_of_frames:d} were expected."
            )

        elif state_reporter_frames > expected_frame_count:
            # Write to a temporary file so that we don't have to read the entire
            # state reporter file into memory
            tmp_file = f"{self.state_reporter_file}.tmp"

            with open(self.state_reporter_file, "r") as input_state_data:
                with open(tmp_file, "w") as output_state_data:
                    # Write header line
                    output_state_data.write(input_state_data.readline())

                    # Write frames up to the expected number from the checkpoint
                    frame_index = 0
                    while frame_index < expected_frame_count:
                        frame_index += 1
                        output_state_data.write(input_state_data.readline())

            # Overwrite the state reporter file with the truncated temporary
            # file
            Path(tmp_file).rename(self.state_reporter_file)

        # Check number of frames in DCD reporter file
        mdtraj_top = mdtraj.load_topology(self.initial_pdb_file)
        dcd_frames = 0
        for traj in mdtraj.iterload(self.dcd_reporter_file, top=mdtraj_top):
            dcd_frames += len(traj)

        if dcd_frames < expected_frame_count:
            raise ValueError(
                f"The DCD reporter file has {dcd_frames:d} frames but "
                f"{expected_number_of_frames:d} were expected."
            )

        elif dcd_frames > expected_frame_count:
            # Write to a temporary file so that we don't have to read the entire
            # DCD file into memory
            tmp_file = f"{self.dcd_reporter_file}.tmp"

            with DCDTrajectoryFile(self.dcd_reporter_file, "r") as input_dcd:
                with DCDTrajectoryFile(tmp_file, "w") as output_dcd:
                    # Write frames up to the expected number from the checkpoint
                    frame_index = 0
                    while frame_index < expected_frame_count:
                        frame_index += 1
                        frame = input_dcd.read_as_traj(mdtraj_top, n_frames=1)

                        output_dcd.write(
                            xyz=in_units_of(
                                frame.xyz,
                                frame._distance_unit,
                                output_dcd.distance_unit,
                            ),
                            cell_lengths=in_units_of(
                                frame.unitcell_lengths,
                                frame._distance_unit,
                                output_dcd.distance_unit,
                            ),
                            cell_angles=frame.unitcell_angles[0],
                        )

            # Overwrite the state reporter file with the truncated temporary
            # file
            Path(tmp_file).rename(self.dcd_reporter_file)

        # Resume dynamics with the checkpointed simulation
        self.run_dynamics(simulation, append=True)

    def run_dynamics(
        self,
        simulation: app.Simulation,
        append: bool = False,
    ):
        """
        Run dynamics for a simulation and write output.

        Parameters
        ----------
        simulation
            An OpenMM Simulation object.
        append
            Append to DCD and state data reporters instead of overwriting them.
        """

        # Set up reporters for DCD trajectory coordinates, state data, and
        # binary checkpoints
        dcd_reporter = app.DCDReporter(
            self.dcd_reporter_file, self.output_frequency, append=append
        )
        state_reporter = app.StateDataReporter(
            self.state_reporter_file,
            self.output_frequency,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            speed=True,
            separator=" ",
            append=append,
        )
        checkpoint_reporter = app.CheckpointReporter(
            self.checkpoint_file, self.checkpoint_frequency
        )

        simulation.reporters.extend([dcd_reporter, state_reporter, checkpoint_reporter])

        # Get current index of serialized simulation state files
        save_state_index = 0
        current_saved_states = []
        if append:
            save_state_dir = Path(self.save_state_prefix).parent
            glob_prefix = Path(self.save_state_prefix).name

            for save_state_file in save_state_dir.glob(f"{glob_prefix}-*.xml"):
                current_saved_states.append(save_state_file)
                file_index = int(save_state_file.stem.split("-")[-1])
                if file_index > save_state_index:
                    save_state_index = file_index

        # Run dynamics until the desired number of steps is reached
        n_saved_states = 200
        # sort saved states files
        current_saved_states.sort(key=lambda x: int(x.stem.split("-")[-1]))
        if len(current_saved_states) > n_saved_states:
            for save_state_file in current_saved_states[:-n_saved_states]:
                save_state_file.unlink()

        while simulation.currentStep < self.n_steps:
            steps_remaining = self.n_steps - simulation.currentStep
            steps_to_take = min(self.save_state_frequency, steps_remaining)

            # Run dynamics
            simulation.step(steps_to_take)

            # Write serialized simulation state
            save_state_index += 1
            save_state_file = f"{self.save_state_prefix}-{save_state_index:05d}.xml"
            simulation.saveState(save_state_file)
            
            if save_state_index > n_saved_states:
                current_saved_states.pop(0).unlink()
