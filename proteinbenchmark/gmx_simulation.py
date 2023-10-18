import subprocess
from pathlib import Path

# import gmxapi as gmx
import numpy
from openff.units import unit


class GMXSimulation:
    """A class representing a simulation in GROMACS."""

    def __init__(
        self,
        gmx_executable: str,
        parametrized_system: str,
        initial_coords_file: str,
        output_prefix: str,
        nonbonded_cutoff: unit.Quantity,
        vdw_switch_width: unit.Quantity,
        temperature: unit.Quantity,
        pressure: unit.Quantity,
        langevin_friction: unit.Quantity,
        barostat_constant: unit.Quantity,
        timestep: unit.Quantity,
        traj_length: unit.Quantity,
        frame_length: unit.Quantity,
    ):
        """
        Initializes the simulation parameters and checks units.

        Parameters
        ----------

        gmx_executable
            Name of GROMACS executable to pass to subprocess.
        parametrized_system
            The path to the parametrized system as a top file.
        initial_coords_file
            Path to GRO file used to set initial coordinates.
        output_prefix
            Path prefix for GMX output files.
        nonbonded_cutoff
            The cutoff for the Lennard-Jones potential and PME direct space
            summation.
        vdw_switch_width
            The distance from the nonbonded cutoff at which to apply the
            switching function.
        temperature
            The target temperature of the Langevin thermostat.
        pressure
            The target pressure of the Monte Carlo barostat.
        langevin_friction
            The collision frequency of the Langevin integrator.
        barostat_constant
            The time constant for pressure coupling
        timestep
            The timestep for the Langevin integrator.
        traj_length
            The length of the trajectory (N_steps = traj_length / timestep).
        frame_length
            The amount of time between writing coordinates and state data to
            disk.
        restaints_present
            Should position restraints be applied?
        """

        self.gmx_executable = gmx_executable
        self.parametrized_system = parametrized_system
        self.initial_coords_file = initial_coords_file
        self.output_prefix = output_prefix

        # Check units of arguments
        if not nonbonded_cutoff.is_compatible_with(unit.nanometer):
            raise ValueError("nonbonded_cutoff does not have units of Length")
        if not vdw_switch_width.is_compatible_with(unit.nanometer):
            raise ValueError("vdw_switch_width does not have units of Length")
        if not temperature.is_compatible_with(unit.kelvin):
            raise ValueError("temperature does not have units of Temperature")
        if not pressure.is_compatible_with(unit.atmosphere):
            raise ValueError("pressure does not have units of Mass Length^-1 Time^-2")
        if not langevin_friction.unit.is_compatible(unit.picosecond**-1):
            raise ValueError("langevin_friction does not have units of Time^-1")
        if not timestep.is_compatible_with(unit.picosecond):
            raise ValueError("timestep does not have units of Time")
        if not traj_length.is_compatible_with(unit.picosecond):
            raise ValueError("traj_length does not have units of Time")
        if not frame_length.is_compatible_with(unit.picosecond):
            raise ValueError("frame_length does not have units of Time")

        # Ensure proper units for output
        self.nonbonded_cutoff = nonbonded_cutoff.m_as(unit.nanometer)
        self.switch_distance = (nonbonded_cutoff - vdw_switch_width).m_as(
            unit.nanometer
        )
        self.temperature = temperature.m_as(unit.kelvin)  # Kelvin
        self.pressure = pressure.m_as(unit.bar)  # atm
        self.thermostat_constant = (1.0 / langevin_friction).m_as(unit.picosecond)  # ps
        self.barostat_constant = barostat_constant.m_as(unit.picosecond)  # ps
        self.timestep = timestep.m_as(unit.picosecond)  # fs
        self.n_steps = int(numpy.round(traj_length / timestep))  # steps
        self.output_frequency = int(numpy.round(frame_length / timestep))  # steps

    def setup_simulation(self, continuation: bool = False):
        """
        Set up a GROMACS simulation by writing the MDP file.

        Parameters
        ----------
        continuation
            Continue from an existing simulation.
        """

        # create mdp file
        mdp_file = self.output_prefix + ".mdp"

        # Write settings to MDP File
        with open(mdp_file, "w") as mdp_file_w:
            if continuation:
                mdp_file_w.write("continuation=yes\ngen_vel = no\n")

            else:
                mdp_file_w.write(
                    "continuation=no\n"
                    "gen_vel = yes\n"
                    f"gen_temp = {self.temperature}\n"
                    "gen_seed = -1\n"
                )

            mdp_file_w.write(
                f"integrator = sd\n"
                f"nsteps={self.n_steps}\n"
                f"dt = {self.timestep}\n"
                f"nstenergy = {self.output_frequency}\n"
                f"nstlog = {self.output_frequency}\n"
                f"nstxout-compressed = {self.output_frequency}\n"
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
                f"tau_t = {self.thermostat_constant}\n"
                f"ref_t = {self.temperature}\n"
                "pbc = xyz\n"
                "DispCorr = EnerPres\n"
                "pcoupl = C-rescale\n"
                "pcoupltype = isotropic\n"
                f"tau_p = {self.barostat_constant}\n"
                f"ref_p = {self.pressure}\n"
                "compressibility = 4.5e-5\n"
                "refcoord_scaling = com\n"
            )

        return mdp_file

    def start_from_pdb(self):
        """
        Start a new simulation initializing positions from a PDB and velocities
        to random samples from a Boltzmann distribution at the simulation
        temperature.
        """

        # Generate TPR
        grompp = subprocess.run(
            [
                self.gmx_executable,
                "grompp",
                "-f",
                mdp_file,
                "-p",
                self.parametrized_system,
                "-c",
                self.initial_coords_file,
                "-maxwarn",
                "2",
                "-o",
                str(self.output_prefix),
            ]
        )

        # Run simulation
        simulation = subprocess.run(
            [
                self.gmx_executable,
                "mdrun",
                "-deffnm",
                str(self.output_prefix),
                "-ntmpi",
                "1",
            ]
        )
        # simulation = gmx.mdrun(
        #    input=grompp.output.file['-o'], runtime_args={'-ntmpi': '1'}
        # )

    def start_from_save_state(
        self,
        checkpoint_file: str,
    ):
        """
        Start a new simulation initializing positions and velocities from a
        GMX checkpoint.

        Parameters
        ----------
        checkpoint_file
            Path to the checkpoint file.
        """

        # Generate MDP
        mdp_file = self.setup_simulation(continuation=True)

        # Generate TPR
        grompp = subprocess.run(
            [
                self.gmx_executable,
                "grompp",
                "-f",
                mdp_file,
                "-p",
                self.parametrized_system,
                "-c",
                self.initial_coords_file,
                "-t",
                checkpoint_file,
                "-maxwarn",
                "2",
                "-o",
                str(self.output_prefix),
            ]
        )

        # Run Simulation
        simulation = subprocess.run(
            [
                self.gmx_executable,
                "mdrun",
                "-deffnm",
                str(self.output_prefix),
                "-ntomp",
                "1",
            ]
        )
        # simulation = gmx.mdrun(
        #    input=grompp.output.file['-o'], runtime_args={'-ntomp': '1'}
        # )

    def resume_from_checkpoint(self):
        """
        Resume an existing simulation from a GMX checkpoint.
        """

        # Run Simulation
        simulation = subprocess.run(
            [
                self.gmx_executable,
                "mdrun",
                "-deffnm",
                str(self.output_prefix),
                "-cpi",
                f"{self.output_prefix}.cpt",
                "-ntomp",
                "1",
            ]
        )
