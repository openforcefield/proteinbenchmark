"""Default parameter values for analysis of trajectories."""
from openmm import unit


DCD_TIME_TO_PICOSECONDS = 0.04888821 * unit.picoseconds

DIHEDRAL_ATOMS = {
    'phi': {
        'atom_names': ['C', 'N', 'CA', 'C'],
        'resid_offsets': [-1, 0, 0, 0],
    },
    'psi': {
        'atom_names': ['N', 'CA', 'C', 'N'],
        'resid_offsets': [0, 0, 0, 1],
    },
}
# Karplus parameters in Hz for NMR scalar couplings from
# Best RB, Zheng W, Mittal J. (2014). J. Chem. Theory Comput. 10, 5113-5124.
BEST_KARPLUS_PARAMETERS = {
    '1j_n_ca': {
        'dihedral': 'psi',
        'delta':    0.0,
        'A':  1.70 / unit.second,
        'B': -0.98 / unit.second,
        'C':  9.51 / unit.second,
        'sigma': 0.59 / unit.second,
    },
    '2j_n_ca': {
        'dihedral': 'prev_psi',
        'delta':    0.0,
        'A': -0.66 / unit.second,
        'B': -1.52 / unit.second,
        'C':  7.85 / unit.second,
        'sigma': 0.50 / unit.second,
    },
    '3j_co_co': {
        'A':  1.36 / unit.second,
        'B': -0.93 / unit.second,
        'C':  0.60 / unit.second,
        'sigma': 0.00 / unit.second,
    },
    '3j_ha_co': {
        'dihedral': 'phi',
        'delta': 120.0,
        'A':  3.72 / unit.second,
        'B': -2.18 / unit.second,
        'C':  1.28 / unit.second,
        'sigma': 0.38 / unit.second,
    },
    '3j_hn_ca': {
        'dihedral': 'phi,prev_psi',
        'C_cos_phi': -0.23 / unit.second,
        'C_cos_psi': -0.20 / unit.second,
        'C_sin_phi':  0.07 / unit.second,
        'C_sin_psi':  0.08 / unit.second,
        'C_cos_phi_cos_psi':  0.07 / unit.second,
        'C_cos_phi_sin_psi':  0.12 / unit.second,
        'C_sin_phi_cos_psi': -0.08 / unit.second,
        'C_sin_phi_sin_psi': -0.14 / unit.second,
        'C_0':   0.54 / unit.second,
        'sigma': 0.10 / unit.second,
    },
    '3j_hn_cb': {
        'dihedral': 'phi',
        'delta':   60.0,
        'A':  3.51 / unit.second,
        'B': -0.53 / unit.second,
        'C':  0.14 / unit.second,
        'sigma': 0.25 / unit.second,
    },
    '3j_hn_co': {
        'dihedral': 'phi',
        'delta':  180.0,
        'A':  4.12 / unit.second,
        'B': -1.10 / unit.second,
        'C':  0.11 / unit.second,
        'sigma': 0.31 / unit.second,
    },
    '3j_hn_ha': {
        'dihedral': 'phi',
        'delta':  -60.0,
        'A':  7.97 / unit.second,
        'B': -1.26 / unit.second,
        'C':  0.63 / unit.second,
        'sigma': 0.42 / unit.second,
    },
}
