"""Default parameter values for MD simulations."""
from openff.units import unit

# Default system setup parameters
SOLVENT_PADDING = 1.2 * unit.nanometer
NONBONDED_CUTOFF = 0.9 * unit.nanometer
VDW_SWITCH_WIDTH = 0.1 * unit.nanometer

# Use larger solvent padding for peptides and disordered proteins
DISORDERED_SOLVENT_PADDING = 1.4 * unit.nanometer

# Default minimization parameters
RESTRAINT_ENERGY_CONSTANT = 1.0 * unit.kilocalorie_per_mole / unit.angstrom**2
ENERGY_TOLERANCE = 100 * unit.kilojoules_per_mole / unit.nanometer

# Default equilibration simulation parameters
EQUIL_TIMESTEP = 1.0 * unit.femtosecond
EQUIL_TRAJ_LENGTH = 1.0 * unit.nanosecond
EQUIL_FRAME_LENGTH = 10.0 * unit.picosecond
EQUIL_OPENMM_LANGEVIN_FRICTION = 5.0 / unit.picosecond
EQUIL_OPENMM_BAROSTAT_FREQUENCY = 5
EQUIL_GMX_BAROSTAT_CONSTANT = 10.0 * unit.picosecond
EQUIL_GMX_THERMOSTAT_CONSTANT = 2.0 * unit.picosecond

# Default production simulation parameters
TIMESTEP = 4.0 * unit.femtosecond
FRAME_LENGTH = 100.0 * unit.picosecond
CHECKPOINT_LENGTH = FRAME_LENGTH * 100
SAVE_STATE_LENGTH = CHECKPOINT_LENGTH * 10
OPENMM_LANGEVIN_FRICTION = 1.0 / unit.picosecond
OPENMM_BAROSTAT_FREQUENCY = 25
GMX_BAROSTAT_CONSTANT = 10.0 * unit.picosecond
GMX_THERMOSTAT_CONSTANT = 2.0 * unit.picosecond

# Production trajectory length for peptides, folded proteins, and disordered
# proteins
PEPTIDE_TRAJ_LENGTH = 500.0 * unit.nanosecond
FOLDED_TRAJ_LENGTH = 10.0 * unit.microsecond
DISORDERED_TRAJ_LENGTH = 30.0 * unit.microsecond
