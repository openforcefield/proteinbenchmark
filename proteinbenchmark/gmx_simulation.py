from pathlib import Path

import numpy
import openmm
from openmm import app, unit
from openmmtools.integrators import LangevinIntegrator

from proteinbenchmark.utilities import exists_and_not_empty, read_xml


class GMXSimulation:
    """A class representing a simulation in GROMACS."""

    def __init__(
        self,
        initial_pdb_file: str,
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
        restraints_present: bool,
    ):
        """
        Initializes the simulation parameters and checks units.

        Parameters
        ----------
        initial_pdb_file
            Path to PDB file used to set initial coordinates.
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
        restaints_present
            Should position restraints be applied?
        """

        self.initial_pdb_file = initial_pdb_file
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
        self.restraints_present = restraints_present

    def setup_simulation(
        self,
    ):
        """
        Set up an GROMACS simulation with a Langevin integrator and a Monte Carlo
        barostat.

        Parameters
        ----------
        return_pdb
            Return PDBFile as well as print TPR
        """
        
        #Create MDP File
        mdpfile = self.save_state_prefix + '-min.mdp'
        mdpfile_w = open(mdpfile, 'w')
        if self.restraints_present == 'NVT' or self.restraints_present == 'NPT':
            mdpfile_w.write('define = -DPOSRES')
        mdpfile_w.write('integrator = bd\n' + 'nsteps=' + str(self.n_steps) + '\n' + 'dt = ' + str(self.timestep) + '\n' + 'nstenergy = ' + str(self.save_state_frequency) + '\n' + 'nstlog = ' + str(self.save_state_frequency) + '\n'
                        'nstxout-compressed = ' + str(self.save_state_frequency) + '\n' + 'continuation = no\n' + 'constraint_algorithm = lincs\n' + 'constraints = h-bond\n' + 'lincs_iter = 1\n' + 'lincs_order = 4' + 'cuttoff-scheme = Verlet\n'
                         'ns_type = grid\n' + 'rlist = 1.0\n' + 'vdwtype = cutoff\n' + 'vdw-modifier = force-switch\n' + 'rvdw-switch = 1.0\n' + 'rvdw = 1.0\n' + 'coulombtype = PME\n' + 'rcoulomb = 1.0\n'  + 'pme_order = 4\n' + 'fourierspacing = 0.16\n' 
                         'bd_fric = ' + str(self.langevin_friction) + '\n' + 'ref_t = ' + str(self.temperature) + '\n' + 'pbc = xyz\n')
        if self.restraints_present == 'NVT':
            mdpfile_w.write('gen_vel = yes\n' + 'gen_temp = ' + str(self.temperature) + '\n' + 'gen_seed = -1')
        else:
            mdpfile_w.write('pcoupl = Parrinello-Rahman\n' + 'pcoupltype = isotropic\n' + 'tau_p = ' + str(self.barostat_frequency) + '\n' + 'ref_p = ' + str(self.pressure) + '\n' + 'compressibility = 4.5e-5' + 'refcoord_scaling = com' + 'gen_vel = no')
        mdpfile_w.close()

        #Generate TPR
        if self.restraints_present == 'NVT':
            out_tprfile = str(self.save_state_prefix) + '-nvt.tpr'
            grompp = gmx.commandline_operation('gmx', 'grompp',
                                   input_files={
                                       '-f': mdpfile,
                                       '-p': str(setup_prefix) + '.top',
                                       '-c': str(setup_prefix) + '-min.gro',
                                       '-r': str(setup_prefix) + '-min.gro',
                                       '-maxwarn': '2',
                                   },
                                   output_files={'-o': out_tprfile})
        elif self.restraints_present == 'NPT':
            out_tprfile = str(self.save_state_prefix) + '-npt.tpr'
            grompp = gmx.commandline_operation('gmx', 'grompp',
                                   input_files={
                                       '-f': mdpfile,
                                       '-p': str(setup_prefix) + '.top',
                                       '-c': str(self.save_state_prefix) + '-nvt.gro',
                                       '-r': str(self.save_state_prefix) + '-nvt.gro',
                                       '-t': str(self.save_state_prefix) + '-nvt.cpt',
                                       '-maxwarn': '2',
                                   },
                                   output_files={'-o': out_tprfile})
        else:
            out_tprfile = str(self.save_state_prefix) + '-npt.tpr'
            grompp = gmx.commandline_operation('gmx', 'grompp',
                                   input_files={
                                       '-f': mdpfile,
                                       '-p': str(setup_prefix) + '.top',
                                       '-c': str(self.save_state_prefix) + '-npt.gro',
                                       '-r': str(self.save_state_prefix) + '-npt.gro',
                                       '-t': str(self.save_state_prefix) + '-npt.cpt',
                                       '-maxwarn': '2',
                                   },
                                   output_files={'-o': out_tprfile})

        
        #Check for errors on TPR File Generation
        if grompp.output.returncode.result() != 0:
            raise Exception(grompp.output.stderr.result())
        
        return grompp

    def start_from_pdb(self):
        """
        Start a new simulation initializing positions from a PDB and velocities
        to random samples from a Boltzmann distribution at the simulation
        temperature.
        """

        # Create an OpenMM simulation
        grompp = self.setup_simulation()

        # Run Simulation
        simulation = gmx.mdrun(input=grompp)
        simulation.run

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
        grompp = self.setup_simulation()

        # Run Simulation
        simulation = gmx.mdrun(input=grompp, runtime_args={'-cpi': save_state_file})
        simulation.run
       
