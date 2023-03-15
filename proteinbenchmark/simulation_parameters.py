"""Default parameter values for MD simulations."""
from openmm import unit


# Default system setup parameters
SOLVENT_PADDING = 1.2 * unit.nanometers
NONBONDED_CUTOFF = 0.9 * unit.nanometers
VDW_SWITCH_WIDTH = 0.1 * unit.nanometers

# Use larger solvent padding for peptides and disordered proteins
DISORDERED_SOLVENT_PADDING = 1.4 * unit.nanometers

# Default minimization parameters
RESTRAINT_ENERGY_CONSTANT = 1.0 * unit.kilocalories_per_mole / unit.angstrom**2

# Default equilibration simulation parameters
EQUIL_LANGEVIN_FRICTION = 5.0 / unit.picoseconds
EQUIL_BAROSTAT_FREQUENCY = 5
EQUIL_TIMESTEP = 1.0 * unit.femtoseconds
EQUIL_TRAJ_LENGTH = 1.0 * unit.nanoseconds
EQUIL_FRAME_LENGTH = 10.0 * unit.picoseconds

# Default production simulation parameters
LANGEVIN_FRICTION = 1.0 / unit.picoseconds
BAROSTAT_FREQUENCY = 25
TIMESTEP = 2.0 * unit.femtoseconds
FRAME_LENGTH = 100.0 * unit.picoseconds
CHECKPOINT_LENGTH = FRAME_LENGTH * 10
SAVE_STATE_LENGTH = CHECKPOINT_LENGTH * 100

# Production trajectory length for peptides, folded proteins, and disordered
# proteins
PEPTIDE_TRAJ_LENGTH = 500.0 * unit.nanoseconds
FOLDED_TRAJ_LENGTH = 10.0 * unit.microseconds
DISORDERED_TRAJ_LENGTH = 30.0 * unit.microseconds
