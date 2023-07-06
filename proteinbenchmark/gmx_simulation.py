from pathlib import Path

import gmxapi as gmx
import numpy
import openmm
from openmm import app, unit
from openmmtools.integrators import LangevinIntegrator
import os

from proteinbenchmark.utilities import exists_and_not_empty, read_xml


class GMXSimulation:
    """A class representing a simulation in GROMACS."""

    def __init__(
        self,
        initial_pdb_file: str,
        save_state_prefix: str,
        setup_prefix:str,
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
        self.save_state_prefix = save_state_prefix
        self.setup_prefix=setup_prefix

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
        
        #create mdp file
        state_dir = str(self.save_state_prefix).rsplit('/', 1)
        if self.restraints_present == 'NVT':
            mdpfile = state_dir[1] + '-nvt.mdp'
        elif self.restraints_present == 'NPT':
            mdpfile = state_dir[1] + '-npt.mdp'
        else:
            mdpfile = state_dir[1] + '-md.mdp'
        #convert timestep from fs to ps
        time_step = str(float(str(self.timestep).split()[0])/1000)
        
        #Write settings to MDP File
        mdpfile_w = open(mdpfile, 'w')
        if self.restraints_present == 'NVT' or self.restraints_present == 'NPT':
            mdpfile_w.write('define=-DPOSRES\n' + 'nsteps=' + str(int(float(self.n_steps)/2)) + '\n')
        else:
            mdpfile_w.write('nsteps=' + str(self.n_steps) + '\n')

        if self.restraints_present == 'NVT':
            mdpfile_w.write('continuation=no\n')
        else:
            mdpfile_w.write('continuation=yes\n')
        mdpfile_w.write('integrator = md\n' + 'dt = ' + str(time_step) + '\n' + 'nstenergy = ' + '5000' + '\n' + 'nstlog = ' + '5000' + '\n'
                        'nstxout-compressed = ' + '5000' + '\n' + 'constraint_algorithm = lincs\n' + 'constraints = h-bonds\n' + 'lincs_iter = 2\n' + 'lincs_order = 4\n' + 'cutoff-scheme = Verlet\n'
                         'ns_type = grid\n' + 'nstlist = 20\n' + 'rlist = 0.9\n' + 'vdwtype = cutoff\n' + 'vdw-modifier = force-switch\n' + 'rvdw-switch = 0.88\n' + 'rvdw = 0.9\n' + 'coulombtype = PME\n' + 'rcoulomb = 1.0\n'  + 'pme_order = 4\n' + 'fourierspacing = 0.16\n' 
                         + 'tcoupl = V-rescale\n' + 'tc-grps = System\n' + 'tau_t = 0.2\n' + 'ref_t = ' + str(self.temperature).split()[0] + '\n' + 'pbc = xyz\n' + 'DispCorr = EnerPres\n')
        if self.restraints_present == 'NVT':
            mdpfile_w.write('gen_vel = yes\n' + 'gen_temp = ' + str(self.temperature).split()[0] + '\n' + 'gen_seed = -1')
        else:
            mdpfile_w.write('pcoupl = Parrinello-Rahman\n' + 'pcoupltype = isotropic\n' + 'tau_p = ' + str(self.barostat_frequency) + '\n' + 'ref_p = ' + str(self.pressure) + '\n' + 'compressibility = 4.5e-5\n' + 'refcoord_scaling = com\n' + 'gen_vel = no\n')
        mdpfile_w.close()
        
        return mdpfile

    def run(self):
        """
        Start a new simulation initializing positions from a PDB and velocities
        to random samples from a Boltzmann distribution at the simulation
        temperature.
        """
        #Change Working directory
        state_dir = str(self.save_state_prefix).rsplit('/', 1)
        if self.restraints_present == 'NVT' or self.restraints_present == 'NPT':
            os.chdir(state_dir[0])
            if not os.path.exists(self.restraints_present):
                os.mkdir(self.restraints_present)
            os.chdir(self.restraints_present)
        else:
            os.chdir(state_dir[0])
        print('Working Directory for Energy Minimization: ' + str(os.getcwd()))
        
        #Generate MDP
        mdpfile = self.setup_simulation()

        #Generate TPR
        setup_dir = str(self.setup_prefix).rsplit('/', 1)
        if self.restraints_present == 'NVT':
            out_tprfile = '../' + state_dir[1] + '-nvt.tpr'
            grompp = gmx.commandline_operation('gmx', 'grompp',
                                   input_files={
                                       '-f': '../' + mdpfile,
                                       '-p': '../../../setup/' + setup_dir[1] + '.top',
                                       '-c': '../../../setup/mdrun_0_i0_0/confout.gro',
                                       '-r': '../../../setup/mdrun_0_i0_0/confout.gro',
                                       '-maxwarn': '2',
                                   },
                                   output_files={'-o': out_tprfile})
        
        elif self.restraints_present == 'NPT':
            out_tprfile = '../' + state_dir[1] + '-npt.tpr'
            grompp = gmx.commandline_operation('gmx', 'grompp',
                                   input_files={
                                       '-f': '../' + mdpfile,
                                       '-p': '../../../setup/' + setup_dir[1] + '.top',
                                       '-c': '../../NVT/mdrun_0_i0_0/confout.gro',
                                       '-r': '../../NVT/mdrun_0_i0_0/confout.gro',
                                       '-t': '../../NVT/mdrun_0_i0_0/state.cpt',
                                       '-maxwarn': '2',
                                   },
                                   output_files={'-o': out_tprfile})
        
        else:
            out_tprfile = state_dir[1] + '-md.tpr'
            grompp = gmx.commandline_operation('gmx', 'grompp',
                                   input_files={
                                       '-f': '../' + mdpfile,
                                       '-p': '../../setup/' + setup_dir[1] + '.top',
                                       '-c': '../NPT/mdrun_0_i0_0/confout.gro',
                                       '-t': '../NPT/mdrun_0_i0_0/state.cpt',
                                       '-maxwarn': '2',
                                   },
                                   output_files={'-o': out_tprfile})

        #Check for errors on TPR File Generation
        if grompp.output.returncode.result() != 0:
            print(grompp.output.stderr.result())
        print('TPR check')

        # Run Simulation
        #simulation = gmx.mdrun(input=grompp, runtime_args={'-ntomp':1, '-ntmpi':2})
        simulation = gmx.mdrun(input=grompp.output.file['-o'])
        simulation.run()

        os.chdir('../../../')
        if self.restraints_present == 'NVT' or self.restraints_present == 'NPT':
            os.chdir('../')

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
       
