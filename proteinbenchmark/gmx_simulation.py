import subprocess
from pathlib import Path

# import gmxapi as gmx
import numpy
from openff.units import unit

from proteinbenchmark.utilities import exists_and_not_empty, read_xml


class GMXSimulation:
    """A class representing a simulation in GROMACS."""

    def __init__(
        self,
        gmx_executable: str,
        parametrized_system: str,
        initial_coords_file: str,
        save_state_prefix: str,
        setup_prefix: str,
        temperature: unit.Quantity,
        pressure: unit.Quantity,
        barostat_constant: unit.Quantity,
        thermostat_constant: unit.Quantity,
        timestep: unit.Quantity,
        traj_length: unit.Quantity,
        frame_length: unit.Quantity,
        restraints_present: bool,
        load_state_prefix: str = None,
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
        save_state_prefix
            Path prefix for files to write serialized simulation states.
        temperature
            The target temperature of the Langevin thermostat.
        pressure
            The target pressure of the Monte Carlo barostat.
        barostat_constant
            The time constant for pressure coupling
        thermostat_constant
            The time constant for temperature coupling
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
        self.save_state_prefix = save_state_prefix
        self.setup_prefix = setup_prefix
        self.load_prefix = load_state_prefix

        # Check units of arguments
        if not temperature.is_compatible_with(unit.kelvin):
            raise ValueError("temperature does not have units of Temperature")

        if not pressure.is_compatible_with(unit.atmosphere):
            raise ValueError("pressure does not have units of Mass Length^-1 Time^-2")

        if not timestep.is_compatible_with(unit.picosecond):
            raise ValueError("timestep does not have units of Time")
        if not traj_length.is_compatible_with(unit.picosecond):
            raise ValueError("traj_length does not have units of Time")
        if not frame_length.is_compatible_with(unit.picosecond):
            raise ValueError("frame_length does not have units of Time")

        # Ensure proper units for output
        self.temperature = temperature.m_as(unit.kelvin)  # Kelvin
        self.pressure = pressure.m_as(unit.atmosphere)  # atm
        self.thermostat_constant = thermostat_constant.m_as(unit.picosecond)  # ps
        self.barostat_constant = barostat_constant.m_as(unit.picosecond)  # ps
        self.timestep = timestep.m_as(unit.picosecond)  # fs
        self.n_steps = int(numpy.round(traj_length / timestep))  # steps
        self.output_frequency = int(numpy.round(frame_length / timestep))  # steps
        self.restraints_present = restraints_present  # str

    def setup_simulation(self):
        """
        Set up a GROMACS simulation.
        """

        # create mdp file
        mdp_file = self.save_state_prefix + ".mdp"

        # Write settings to MDP File
        with open(mdp_file, "w") as mdp_file_w:
            # Apply position restraints for equilibration
            if self.restraints_present == "NPT":
                mdp_file_w.write("define=-DPOSRES\n" + "continuation=no\n")
            else:
                mdp_file_w.write("continuation=yes\n")

            mdp_file_w.write(
                f"integrator = md\nnsteps={self.n_steps}\n"
                f"dt = {self.timestep}\nnstenergy = 5000\nnstlog = 5000\n"
                f"nstxout-compressed = {self.output_frequency}\n"
                "constraint_algorithm = lincs\nlincs-warnangle = 45\n"
                "constraints = h-bonds\nlincs_iter = 2\nlincs_order = 4\n"
                "cutoff-scheme = Verlet\nns_type = grid\nnstlist = 40\n"
                "rlist = 0.9\nvdwtype = cutoff\nvdw-modifier = force-switch\n"
                "rvdw-switch = 0.88\nrvdw = 0.9\ncoulombtype = PME\n"
                "rcoulomb = 0.9\npme_order = 4\nfourierspacing = 0.16\n"
                "tcoupl = V-rescale\ntc-grps = System\n"
                f"tau_t = {self.thermostat_constant}\n"
                f"ref_t = {self.temperature}\npbc = xyz\nDispCorr = EnerPres\n"
                f"pcoupltype = isotropic\ntau_p = {self.barostat_constant}\n"
                f"ref_p = {self.pressure}\ncompressibility = 4.5e-5\n"
                "refcoord_scaling = com\n"
            )
            if self.restraints_present == "NPT":
                mdp_file_w.write(
                    f"gen_vel = yes\ngen_temp = {self.temperature}\n"
                    "gen_seed = -1\npcoupl = C-rescale\n"
                )
            else:
                mdp_file_w.write("gen_vel = no\npcoupl = Parrinello-Rahman\n")

        return mdp_file

    def run(self):
        """
        Start a new simulation initializing positions from a PDB and velocities
        to random samples from a Boltzmann distribution at the simulation
        temperature.
        """

        # Generate MDP
        mdp_file = self.setup_simulation()

        # Generate TPR
        if self.restraints_present == "NPT":
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
                    "-r",
                    self.initial_coords_file,
                    "-maxwarn",
                    "2",
                    "-o",
                    str(self.save_state_prefix),
                ]
            )

        else:
            print(self.gmx_executable)
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
                    f"{self.load_prefix}.cpt",
                    "-maxwarn",
                    "2",
                    "-o",
                    str(self.save_state_prefix),
                ]
            )

        # Run Simulation
        if self.restraints_present == "NPT":
            simulation = subprocess.run(
                [
                    self.gmx_executable,
                    "mdrun",
                    "-deffnm",
                    str(self.save_state_prefix),
                    "-ntmpi",
                    "1",
                ]
            )
            # simulation = gmx.mdrun(
            #    input=grompp.output.file['-o'], runtime_args={'-ntmpi': '1'}
            # )
        else:
            simulation = subprocess.run(
                [
                    self.gmx_executable,
                    "mdrun",
                    "-deffnm",
                    str(self.save_state_prefix),
                    "-ntomp",
                    "1",
                ]
            )
            # simulation = gmx.mdrun(
            #    input=grompp.output.file['-o'], runtime_args={'-ntomp': '1'}
            # )

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
        # Run Simulation
        simulation = subprocess.run(
            [
                self.gmx_executable,
                "mdrun",
                "-deffnm",
                str(self.save_state_prefix),
                "-cpi",
                save_state_file,
                "-ntomp",
                "1",
            ]
        )
