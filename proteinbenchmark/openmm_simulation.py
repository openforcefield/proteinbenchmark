import dataclasses
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Literal, Self, overload

import numpy
import openmm
import openmm.app
import openmm.unit
import openmmtools.mcmc
import openmmtools.multistate
from openff.toolkit import Quantity
from openff.units import ensure_quantity
from openmm import (
    app,
    unit,
)

from proteinbenchmark.utilities import exists_and_not_empty, read_xml, write_xml


class FrameCountMismatchError(Exception):
    pass


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
            f"\n    save_state_frequency {self.save_state_frequency:d} steps",
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
                ),
            )

        # Create simulation
        simulation = app.Simulation(
            initial_pdb.topology,
            openmm_system,
            integrator,
            openmm.Platform.getPlatformByName("CUDA"),
            {"Precision": "mixed"},
        )

        if return_pdb:
            return simulation, initial_pdb
        else:
            return simulation

    def check_frame_count(self, simulation: app.Simulation):
        """
        Check the number of frames in the state data reporter and the DCD
        reporter match the expected number of frames from the OpenMM simulation.

        Parameters
        ---------
        simulation
            An OpenMM Simulation object.
        """

        import mdtraj
        from mdtraj.formats.dcd import DCDTrajectoryFile
        from mdtraj.utils import in_units_of

        # Get expected number of frames based on current step from checkpoint
        # If the state data or DCD reporters have additional frames written,
        # truncate them to the expected number of frames
        expected_frame_count = int(simulation.currentStep / self.output_frequency)

        # Check number of frames in state data reporter file
        with open(self.state_reporter_file, "r") as state_reporter:
            # Subtract one for header line
            state_reporter_frames = sum(1 for _ in state_reporter) - 1

        if state_reporter_frames < expected_frame_count:
            raise FrameCountMismatchError(
                f"The state data reporter file has {state_reporter_frames:d} "
                f"frames but {expected_frame_count:d} were expected.",
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
            raise FrameCountMismatchError(
                f"The DCD reporter file has {dcd_frames:d} frames but "
                f"{expected_frame_count:d} were expected.",
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
        resume: bool = False,
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
                "exist or is empty.",
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

        # Create an OpenMM simulation
        simulation = self.setup_simulation()

        # Load the checkpoint
        if not exists_and_not_empty(self.checkpoint_file):
            raise ValueError(
                f"Checkpoint file {self.checkpoint_file} does not exist or is empty.",
            )

        simulation.loadCheckpoint(self.checkpoint_file)

        # Check whether the simulation has already finished
        if simulation.currentStep == self.n_steps:
            return

        # Check frame counts between simulation, state data reporter, and DCD
        # reporter
        try:
            self.check_frame_count(simulation)

        except FrameCountMismatchError:
            # Load the last serialized saved state XML
            save_state_index = 0
            save_state_dir = Path(self.save_state_prefix).parent
            glob_prefix = Path(self.save_state_prefix).name

            for save_state_file in save_state_dir.glob(f"{glob_prefix}-*.xml"):
                file_index = int(save_state_file.stem.split("-")[-1])
                if file_index > save_state_index:
                    save_state_index = file_index

            save_state_file = f"{self.save_state_prefix}-{save_state_index}.xml"
            simulation.loadState(save_state_file)

            # Truncate both the state data reporter and DCD reporter to the
            # number of frames in the saved state
            self.check_frame_count(simulation)

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
            self.dcd_reporter_file,
            self.output_frequency,
            append=append,
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
            self.checkpoint_file,
            self.checkpoint_frequency,
        )

        simulation.reporters.extend([dcd_reporter, state_reporter, checkpoint_reporter])

        # Get current index of serialized simulation state files
        save_state_index = 0
        if append:
            save_state_dir = Path(self.save_state_prefix).parent
            glob_prefix = Path(self.save_state_prefix).name

            for save_state_file in save_state_dir.glob(f"{glob_prefix}-*.xml"):
                file_index = int(save_state_file.stem.split("-")[-1])
                if file_index > save_state_index:
                    save_state_index = file_index

        # Run dynamics until the desired number of steps is reached
        while simulation.currentStep < self.n_steps:
            steps_remaining = self.n_steps - simulation.currentStep
            steps_at_next_save_state = int(
                numpy.ceil((simulation.currentStep + 1) / self.save_state_frequency)
                * self.save_state_frequency,
            )
            steps_to_next_save_state = steps_at_next_save_state - simulation.currentStep
            steps_to_take = min(steps_to_next_save_state, steps_remaining)

            # Run dynamics
            simulation.step(steps_to_take)

            # Write serialized simulation state
            save_state_index += 1
            save_state_file = f"{self.save_state_prefix}-{save_state_index}.xml"
            simulation.saveState(save_state_file)


@dataclasses.dataclass
class OpenMMHrexEnsemble:
    """An H-REx ensemble of simulations with Gibbs sampling using openmmtools"""

    base_simulation: OpenMMSimulation
    system_fn_ladder: list[Path]
    steps_between_exchange_attempts: int

    @classmethod
    def construct_rest2(
        cls,
        n_replicas: int,
        base_simulation: OpenMMSimulation,
        tempered_atom_idcs: Iterable[int],
        steps_between_exchange_attempts: int,
        max_effective_temperature: unit.Quantity | Quantity = unit.Quantity(
            600.0,
            "kelvin",
        ),
    ) -> Self:
        t_0 = base_simulation.temperature
        t_max = ensure_quantity(max_effective_temperature, "openmm")
        effective_temperatures: list[openmm.unit.Quantity] = [
            t_0 * (t_max / t_0) ** (m / (n_replicas - 1)) for m in range(n_replicas)
        ]
        assert effective_temperatures[0] == t_0
        assert effective_temperatures[-1] == t_max
        assert len(effective_temperatures) == n_replicas

        tempered_atom_idcs = set(tempered_atom_idcs)

        ladder: list[Path] = []

        for t_m, m in zip(
            effective_temperatures,
            range(n_replicas),
            strict=True,
        ):
            openmm_system = read_xml(base_simulation.openmm_system_file)
            if m != 0:
                for force in openmm_system.getForces():
                    _rest2_scale_force(
                        force,
                        tempered_atom_idcs=tempered_atom_idcs,
                        base_temperature=t_0,
                        effective_temperature=t_m,
                    )

            base_system_path = Path(base_simulation.openmm_system_file)
            rung_system_fn = base_system_path.with_name(
                f"{base_system_path.stem}-rung{m}-{t_m.value_in_unit(unit.kelvin):3.2f}k.openmm.xml",
            )
            write_xml(rung_system_fn, openmm_system)
            ladder.append(rung_system_fn)

        return cls(
            base_simulation=base_simulation,
            system_fn_ladder=ladder,
            steps_between_exchange_attempts=steps_between_exchange_attempts,
        )

    @overload
    def setup_simulation(
        self,
        return_pdb: Literal[True],
    ) -> tuple[openmmtools.multistate.ReplicaExchangeSampler, openmm.app.PDBFile]: ...

    @overload
    def setup_simulation(
        self,
        return_pdb: Literal[False] = False,
    ) -> openmmtools.multistate.ReplicaExchangeSampler: ...

    def setup_simulation(
        self,
        return_pdb: bool = False,
    ) -> (
        openmmtools.multistate.ReplicaExchangeSampler
        | tuple[openmmtools.multistate.ReplicaExchangeSampler, openmm.app.PDBFile]
    ):
        """
        Set up an OpenMM simulation with a Langevin integrator and a Monte Carlo
        barostat.

        Parameters
        ----------
        return_pdb
            Return OpenMM PDBFile as well as OpenMM Simulation.
        """
        if self.base_simulation.n_steps % self.steps_between_exchange_attempts != 0:
            raise ValueError(
                "steps_between_exchange_attempts must evenly divide n_steps",
            )
        if (
            self.base_simulation.checkpoint_frequency
            % self.steps_between_exchange_attempts
            != 0
        ):
            raise ValueError(
                "steps_between_exchange_attempts must evenly divide checkpoint_frequency",
            )
        if (
            self.base_simulation.output_frequency % self.steps_between_exchange_attempts
            != 0
        ):
            raise ValueError(
                "steps_between_exchange_attempts must evenly divide output_frequency",
            )

        # Load OpenMM systems and initial PDB
        openmm_systems = [read_xml(fn) for fn in self.system_fn_ladder]
        initial_pdb = app.PDBFile(self.base_simulation.initial_pdb_file)

        # Set up BAOAB Langevin integrator from openmmtools with VRORV splitting
        integrator = openmm.LangevinMiddleIntegrator(
            self.base_simulation.temperature,
            self.base_simulation.langevin_friction,
            self.base_simulation.timestep,
        )

        # Set up Monte Carlo barostat
        if self.base_simulation.pressure.value_in_unit(unit.atmosphere) > 0:
            for system in openmm_systems:
                system.addForce(
                    openmm.MonteCarloBarostat(
                        self.base_simulation.pressure,
                        self.base_simulation.temperature,
                        self.base_simulation.barostat_frequency,
                    ),
                )

        # Create sampler
        sampler = openmmtools.multistate.ReplicaExchangeSampler(
            mcmc_moves=openmmtools.mcmc.IntegratorMove(
                integrator=integrator,
                n_steps=self.steps_between_exchange_attempts,
            ),
            replica_mixing_scheme="swap-all",
            number_of_iterations=(
                self.base_simulation.n_steps // self.steps_between_exchange_attempts
            ),
        )

        # Configure reporters
        if (
            self.base_simulation.dcd_reporter_file
            == self.base_simulation.state_reporter_file
        ):
            analysis_frequency = self.steps_between_exchange_attempts
            storage = openmmtools.multistate.MultiStateReporter(
                storage=self.base_simulation.dcd_reporter_file,
                checkpoint_storage=self.base_simulation.checkpoint_file,
                checkpoint_interval=self.base_simulation.checkpoint_frequency
                // analysis_frequency,
                position_interval=self.base_simulation.output_frequency
                // analysis_frequency,
                velocity_interval=self.base_simulation.output_frequency
                // analysis_frequency,
            )
        else:
            raise ValueError("DCD and state reporter files must have same name")

        # Initialize with states
        sampler.create(
            thermodynamic_states=[
                openmmtools.states.ThermodynamicState(
                    system=system,
                    temperature=integrator.getTemperature(),
                )
                for system in openmm_systems
            ],
            sampler_states=[
                openmmtools.states.SamplerState(
                    positions=initial_pdb.getPositions(),
                    box_vectors=initial_pdb.getTopology().getPeriodicBoxVectors(),
                )
                for system in openmm_systems
            ],
            storage=storage,
        )
        assert sampler.n_replicas == len(openmm_systems)

        # Specify platform and mixed precision
        try:
            platform = openmm.Platform.getPlatformByName("CUDA")
        except openmm.OpenMMException:
            platform = openmm.Platform.getPlatformByName("HIP")
        context_cache_dict = dict(
            platform=platform,
            platform_properties={"Precision": "mixed"},
            capacity=None,
            time_to_live=None,
        )
        sampler.energy_context_cache = openmmtools.cache.ContextCache(
            **context_cache_dict,
        )
        sampler.sampler_context_cache = openmmtools.cache.ContextCache(
            **context_cache_dict,
        )

        if return_pdb:
            return sampler, initial_pdb
        else:
            return sampler

    def start_from_storage(
        self,
        storage: str,
    ) -> openmmtools.multistate.ReplicaExchangeSampler:
        sampler = openmmtools.multistate.ReplicaExchangeSampler.from_storage(storage)
        assert isinstance(sampler, openmmtools.multistate.ReplicaExchangeSampler)
        assert sampler.n_replicas == len(self.system_fn_ladder)
        n_iterations = sampler.number_of_iterations
        assert len(sampler.mcmc_moves) == 1
        assert sampler.mcmc_moves[0].n_steps == self.steps_between_exchange_attempts
        assert (
            n_iterations * self.steps_between_exchange_attempts
            == self.base_simulation.n_steps
        )
        return sampler


def _rest2_scale_force(
    force: openmm.Force,
    *,
    tempered_atom_idcs: set[int],
    base_temperature: openmm.unit.Quantity,
    effective_temperature: openmm.unit.Quantity,
    scale_non_torsion_bonded_interactions: bool = True,
) -> openmm.Force:
    beta_0 = 1 / base_temperature
    beta_m = 1 / effective_temperature
    match force:
        case (
            openmm.CMMotionRemover()
            | openmm.AndersenThermostat()
            | openmm.MonteCarloAnisotropicBarostat()
            | openmm.MonteCarloBarostat()
            | openmm.MonteCarloFlexibleBarostat()
            | openmm.MonteCarloMembraneBarostat()
        ):
            # No scaling for non-force Forces
            pass
        case openmm.NonbondedForce():
            # solute-solute should be scaled by:                             beta_m/beta_0
            # solute-solvent should be scaled by:                            sqrt(beta_m/beta_0)
            # solvent-solvent should not be scaled
            # electrostatics are of form q1*q2*<other stuff>
            # scaling solute charges by sqrt(beta_m/beta_0) gives
            #   solute-solute:   sqrt(beta_m/beta_0) * sqrt(beta_m/beta_0) = beta_m/beta_0
            #   solute-solvent:                    sqrt(beta_m/beta_0) * 1 = sqrt(beta_m/beta_0)
            # LJs are of form sqrt(eps1*eps2)*<other stuff>
            # scaling solute epsilons by beta_m/beta_0 gives
            #   solute-solute:         sqrt(beta_m/beta_0 * beta_m/beta_0) = beta_m/beta_0
            #   solute-solvent:                    sqrt(beta_m/beta_0 * 1) = sqrt(beta_m/beta_0)
            # solvent parameters are unchanged so solvent-solvent interactions are unscaled
            for i in tempered_atom_idcs:
                charge, sigma, epsilon = force.getParticleParameters(i)
                charge_tempered = math.sqrt(beta_m / beta_0) * charge
                epsilon_tempered = (beta_m / beta_0) * epsilon
                force.setParticleParameters(i, charge_tempered, sigma, epsilon_tempered)
        case openmm.PeriodicTorsionForce():
            # Energies are proportional to K
            # each parameter set in the force fully defines the relevent particles
            # therefore we just multiply K by the desired scale for the parameter
            for i in range(force.getNumTorsions()):
                p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(i)
                torsion_particles = {p1, p2, p3, p4}
                if torsion_particles.isdisjoint(tempered_atom_idcs):
                    # no particles are tempered, this is solvent-solvent
                    k_tempered = k
                elif torsion_particles.issubset(tempered_atom_idcs):
                    # all particles are tempered, this is solute-solute
                    k_tempered = k * beta_m / beta_0
                else:
                    # some but not all particles are tempered, this is solute-solvent
                    k_tempered = k * math.sqrt(beta_m / beta_0)
                force.setTorsionParameters(
                    i,
                    p1,
                    p2,
                    p3,
                    p3,
                    periodicity,
                    phase,
                    k_tempered,
                )
        case openmm.HarmonicAngleForce() if scale_non_torsion_bonded_interactions:
            # Energies are proportional to K
            # each parameter set in the force fully defines the relevent particles
            # therefore we just multiply K by the desired scale for the parameter
            for i in range(force.getNumAngles()):
                p1, p2, p3, angle, k = force.getAngleParameters(i)
                angle_particles = {p1, p2, p3}
                if angle_particles.isdisjoint(tempered_atom_idcs):
                    # no particles are tempered, this is solvent-solvent
                    k_tempered = k
                elif angle_particles.issubset(tempered_atom_idcs):
                    k_tempered = k * beta_m / beta_0
                else:
                    k_tempered = k * math.sqrt(beta_m / beta_0)
                force.setAngleParameters(i, p1, p2, p3, angle, k_tempered)
        case openmm.HarmonicBondForce() if scale_non_torsion_bonded_interactions:
            # Energies are proportional to K
            # each parameter set in the force fully defines the relevent particles
            # therefore we just multiply K by the desired scale for the parameter
            for i in range(force.getNumBonds()):
                p1, p2, length, k = force.getBondParameters(i)
                bond_particles = {p1, p2}
                if bond_particles.isdisjoint(tempered_atom_idcs):
                    # no particles are tempered, this is solvent-solvent
                    k_tempered = k
                elif bond_particles.issubset(tempered_atom_idcs):
                    # all particles are tempered, this is solute-solute
                    k_tempered = k * beta_m / beta_0
                else:
                    # some but not all particles are tempered, this is solute-solvent
                    k_tempered = k * math.sqrt(beta_m / beta_0)
                force.setBondParameters(i, p1, p2, length, k_tempered)
        case openmm.HarmonicAngleForce() | openmm.HarmonicBondForce():
            # User has asked not to scale these forces, so leave them alone
            assert not scale_non_torsion_bonded_interactions
            pass
        case _:
            raise ValueError(f"Cannot scale force {force}")
    return force
