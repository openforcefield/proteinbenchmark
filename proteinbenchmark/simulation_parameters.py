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

# Default equilibration simulation parameters
EQUIL_LANGEVIN_FRICTION = 5.0 / unit.picosecond
EQUIL_BAROSTAT_FREQUENCY = 5
EQUIL_TIMESTEP = 1.0 * unit.femtosecond
EQUIL_TRAJ_LENGTH = 1.0 * unit.nanosecond
EQUIL_FRAME_LENGTH = 10.0 * unit.picosecond

# Default production simulation parameters
LANGEVIN_FRICTION = 1.0 / unit.picosecond
BAROSTAT_FREQUENCY = 25
TIMESTEP = 4.0 * unit.femtosecond
FRAME_LENGTH = 100.0 * unit.picosecond
CHECKPOINT_LENGTH = FRAME_LENGTH * 100
SAVE_STATE_LENGTH = CHECKPOINT_LENGTH * 10

# Production trajectory length for peptides, folded proteins, and disordered
# proteins
PEPTIDE_TRAJ_LENGTH = 500.0 * unit.nanosecond
FOLDED_TRAJ_LENGTH = 10.0 * unit.microsecond
DISORDERED_TRAJ_LENGTH = 30.0 * unit.microsecond
