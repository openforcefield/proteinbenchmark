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
        barostat_constant: unit.Quantity,
        thermostat_constant: unit.Quantity,
        timestep: unit.Quantity,
        traj_length: unit.Quantity,
        frame_length: unit.Quantity,
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

        self.initial_pdb_file = initial_pdb_file
        self.save_state_prefix = save_state_prefix
        self.setup_prefix=setup_prefix

        # Check units of arguments
        if not temperature.unit.is_compatible(unit.kelvin):
            raise ValueError("temperature does not have units of Temperature")

        if not pressure.unit.is_compatible(unit.atmosphere):
            raise ValueError("pressure does not have units of Mass Length^-1 Time^-2")

        if not timestep.unit.is_compatible(unit.picosecond):
            raise ValueError("timestep does not have units of Time")
        if not traj_length.unit.is_compatible(unit.picosecond):
            raise ValueError("traj_length does not have units of Time")
        if not frame_length.unit.is_compatible(unit.picosecond):
            raise ValueError("frame_length does not have units of Time")

        #Ensure proper units for output
        self.temperature = temperature.value_in_unit(unit.kelvin) #Kelvin
        self.pressure = pressure.value_in_unit(unit.atmosphere) #atm
        self.thermostat_constant = thermostat_constant.value_in_unit(unit.picosecond) #ps
        self.barostat_constant = barostat_constant.value_in_unit(unit.picosecond) #ps
        self.timestep = timestep.value_in_unit(unit.picosecond) #fs
        self.n_steps = int(numpy.round(traj_length.value_in_unit(unit.picosecond) / timestep.value_in_unit(unit.picosecond))) #steps
        self.output_frequency = int(numpy.round(frame_length.value_in_unit(unit.picosecond) / timestep.value_in_unit(unit.picosecond))) #steps
        self.restraints_present = restraints_present #str

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
        if self.restraints_present == 'NPT':
            mdpfile = state_dir[1] + '.mdp'
        else:
            mdpfile = state_dir[1] + '.mdp'
        
        #Write settings to MDP File
        mdpfile_w = open(mdpfile, 'w')
        if self.restraints_present == 'NPT': #Apply position restraints for equilibration
            mdpfile_w.write('define=-DPOSRES\n' + 'continuation=no\n')
        else:
            mdpfile_w.write('continuation=yes\n')

        mdpfile_w.write('integrator = md\n' + f'nsteps={self.n_steps}\n' + f'dt = {self.timestep}\n' + 'nstenergy = 5000\n' + 'nstlog = 5000\n' + 
                        f'nstxout-compressed = {self.output_frequency}\n' + 'constraint_algorithm = lincs\n' + 'constraints = h-bonds\n' + 
                        'lincs_iter = 2\n' + 'lincs_order = 4\n' + 'cutoff-scheme = Verlet\n' + 'ns_type = grid\n' + 'nstlist = 40\n' + 
                        'rlist = 0.9\n' + 'vdwtype = cutoff\n' + 'vdw-modifier = force-switch\n' + 'rvdw-switch = 0.88\n' + 'rvdw = 0.9\n' + 
                        'coulombtype = PME\n' + 'rcoulomb = 0.9\n'  + 'pme_order = 4\n' + 'fourierspacing = 0.16\n' + 'tcoupl = V-rescale\n' + 
                        'tc-grps = System\n' + f'tau_t = {self.thermostat_constant}\n' + f'ref_t = {self.temperature}\n' + 'pbc = xyz\n' + 
                        'DispCorr = EnerPres\n' + 'pcoupltype = isotropic\n' + f'tau_p = {self.barostat_constant}\n' + 
                        f'ref_p = {self.pressure}\n' + 'compressibility = 4.5e-5\n' + 'refcoord_scaling = com\n')
        if self.restraints_present == 'NPT':
            mdpfile_w.write('gen_vel = yes\n' + f'gen_temp = {self.temperature}\n' + 'gen_seed = -1\n' + 'pcoupl = C-rescale\n')
        else:
            mdpfile_w.write('gen_vel = no\n' + 'pcoupl = Parrinello-Rahman\n')
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
        if self.restraints_present == 'NPT':
            os.chdir(state_dir[0])
            if not os.path.exists(self.restraints_present):
                os.mkdir(self.restraints_present)
            os.chdir(self.restraints_present)
        else:
            os.chdir(state_dir[0])
        print(f'Working Directory for Energy Minimization: {os.getcwd()}')
        
        #Generate MDP
        mdpfile = self.setup_simulation()

        #Generate TPR
        setup_dir = str(self.setup_prefix).rsplit('/', 1)
        if self.restraints_present == 'NPT':
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

        # Run Simulation
        if self.restraints_present == 'NPT':
            simulation = gmx.mdrun(input=grompp.output.file['-o'], runtime_args={'-ntmpi': '1'})
        else:
            simulation = gmx.mdrun(input=grompp.output.file['-o'], runtime_args={'-ntomp': '1'})
        simulation.run()

        os.chdir('../../../')
        if self.restraints_present == 'NPT':
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
       
