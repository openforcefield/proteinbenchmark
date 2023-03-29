"""Default parameter values for analysis of trajectories."""
from openmm import unit


DCD_TIME_TO_PICOSECONDS = 0.04888821 * unit.picoseconds

DIHEDRAL_ATOMS = {
    residue: {
        'phi': {
            'atom_names': ['C', 'N', 'CA', 'C'],
            'resid_offsets': [-1, 0, 0, 0],
        },
        'psi': {
            'atom_names': ['N', 'CA', 'C', 'N'],
            'resid_offsets': [0, 0, 0, 1],
        },
        'omega': {
            'atom_names': ['CA', 'C', 'N', 'CA'],
            'resid_offsets': [0, 0, 1, 1],
        },
    }
    for residue in [
        'ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYS', 'CYX', 'GLH', 'GLN', 'GLU',
        'GLY', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE',
        'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    ]
}
DIHEDRAL_ATOMS['ARG']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ASH']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ASN']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ASP']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['CYS']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'SG']}
DIHEDRAL_ATOMS['CYX']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'SG']}
DIHEDRAL_ATOMS['GLH']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['GLN']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['GLU']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['HID']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['HIE']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['HIP']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ILE']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG1']}
DIHEDRAL_ATOMS['LEU']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['LYN']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['LYS']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['MET']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['PHE']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['PRO']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['SER']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'OG']}
DIHEDRAL_ATOMS['THR']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'OG1']}
DIHEDRAL_ATOMS['TRP']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['TYR']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['VAL']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG1']}
DIHEDRAL_ATOMS['ARG']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['ASH']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'OD1']}
DIHEDRAL_ATOMS['ASN']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'OD1']}
DIHEDRAL_ATOMS['ASP']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'OD1']}
DIHEDRAL_ATOMS['GLH']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['GLN']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['GLU']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['HID']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'ND1']}
DIHEDRAL_ATOMS['HIE']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'ND1']}
DIHEDRAL_ATOMS['HIP']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'ND1']}
DIHEDRAL_ATOMS['ILE']['chi2'] = {'atom_names': ['CA', 'CB', 'CG1', 'CD1']}
DIHEDRAL_ATOMS['LEU']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['LYN']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['LYS']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['MET']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'SD']}
DIHEDRAL_ATOMS['PHE']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['PRO']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['TRP']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['TYR']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['ARG']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'NE']}
DIHEDRAL_ATOMS['GLH']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'OE1']}
DIHEDRAL_ATOMS['GLN']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'OE1']}
DIHEDRAL_ATOMS['GLU']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'OE1']}
DIHEDRAL_ATOMS['LYN']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'CE']}
DIHEDRAL_ATOMS['LYS']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'CE']}
DIHEDRAL_ATOMS['MET']['chi3'] = {'atom_names': ['CB', 'CG', 'SD', 'CE']}
DIHEDRAL_ATOMS['ARG']['chi4'] = {'atom_names': ['CG', 'CD', 'NE', 'CZ']}
DIHEDRAL_ATOMS['LYN']['chi4'] = {'atom_names': ['CG', 'CD', 'CE', 'NZ']}
DIHEDRAL_ATOMS['LYS']['chi4'] = {'atom_names': ['CG', 'CD', 'CE', 'NZ']}
DIHEDRAL_ATOMS['ARG']['chi5'] = {'atom_names': ['CD', 'NE', 'CZ', 'NH1']}

# Clusters on Ramachandran map from
# Hollingsworth SA, Karplus PA. (2010). BioMol Concepts 1, 271-283.
RAMACHANDRAN_CLUSTERS = {
    '\alpha': {'phi': [285, 315], 'psi': [-60, -30]},
    '\beta': {'phi': [180, 270], 'psi': [105, 195]},
    '\gamma': {'phi': [65, 105], 'psi': [-90, -30]},
    '\delta': {'phi': [225, 315], 'psi': [-60, 45]},
    '\varepsilon': {'phi': [40, 185], 'psi': [120, 240]},
    '\zeta': {'phi': [195, 255], 'psi': [45, 105]},
    'P_{II}': {'phi': [270, 315], 'psi': [120, 195]},
    "\gamma'": {'phi': [255, 300], 'psi': [45, 105]},
    "\delta'": {'phi': [30, 120], 'psi': [-15, 75]},
    "P_{II}'": {'phi': [45, 120], 'psi': [165, 240]},
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
        'dihedral': 'phi',
        'delta':    0.0,
        'A':  1.36 / unit.second,
        'B': -0.93 / unit.second,
        'C':  0.60 / unit.second,
        'sigma': 0.22 / unit.second,
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

# Karplus parameters in Hz for NMR scalar couplings from
# Graf J, Nguyen PH, Stock G, Schwalbe H. (2007). J. Am. Chem. Soc. 129,
#     1179-1189.
GRAF_KARPLUS_PARAMETERS = {
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
        'dihedral': 'phi',
        'delta':    0.0,
        'A':  1.36 / unit.second,
        'B': -0.93 / unit.second,
        'C':  0.60 / unit.second,
        'sigma': 0.22 / unit.second,
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
        'A':  3.06 / unit.second,
        'B': -0.74 / unit.second,
        'C':  0.13 / unit.second,
        'sigma': 0.39 / unit.second,
    },
    '3j_hn_co': {
        'dihedral': 'phi',
        'delta':  180.0,
        'A':  4.29 / unit.second,
        'B': -1.01 / unit.second,
        'C':  0.00 / unit.second,
        'sigma': 0.59 / unit.second,
    },
    '3j_hn_ha': {
        'dihedral': 'phi',
        'delta':  -60.0,
        'A':  7.09 / unit.second,
        'B': -1.42 / unit.second,
        'C':  1.55 / unit.second,
        'sigma': 0.91 / unit.second,
    },
}

# Karplus parameters in Hz for NMR scalar couplings of Ace-Ala-NMe from
# Case DA, Scheurer C, Bruschweiler R. (2000). J. Am. Chem. Soc. 122,
#     10390-10397.
DFT1_KARPLUS_PARAMETERS = {
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
        'dihedral': 'phi',
        'delta':    0.0,
        'A':  2.39 / unit.second,
        'B': -1.25 / unit.second,
        'C':  0.26 / unit.second,
        'sigma': 0.22 / unit.second,
    },
    '3j_ha_co': {
        'dihedral': 'phi',
        'delta': 120.0,
        'A':  4.38 / unit.second,
        'B': -1.87 / unit.second,
        'C':  0.56 / unit.second,
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
        'A':  5.15 / unit.second,
        'B':  0.01 / unit.second,
        'C': -0.32 / unit.second,
        'sigma': 0.39 / unit.second,
    },
    '3j_hn_co': {
        'dihedral': 'phi',
        'delta':  180.0,
        'A':  5.58 / unit.second,
        'B': -1.06 / unit.second,
        'C': -0.30 / unit.second,
        'sigma': 0.59 / unit.second,
    },
    '3j_hn_ha': {
        'dihedral': 'phi',
        'delta':  -60.0,
        'A':  9.44 / unit.second,
        'B': -1.53 / unit.second,
        'C': -0.07 / unit.second,
        'sigma': 0.91 / unit.second,
    },
}

# Karplus parameters in Hz for NMR scalar couplings of Ala-Ala-NH2 from
# Case DA, Scheurer C, Bruschweiler R. (2000). J. Am. Chem. Soc. 122,
#     10390-10397.
DFT2_KARPLUS_PARAMETERS = {
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
        'dihedral': 'phi',
        'delta':   -2.56,
        'A':  2.71 / unit.second,
        'B': -0.91 / unit.second,
        'C':  0.21 / unit.second,
        'sigma': 0.22 / unit.second,
    },
    '3j_ha_co': {
        'dihedral': 'phi',
        'delta': 118.61,
        'A':  4.77 / unit.second,
        'B': -1.85 / unit.second,
        'C':  0.49 / unit.second,
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
        'delta':   58.18,
        'A':  4.58 / unit.second,
        'B': -0.36 / unit.second,
        'C': -0.31 / unit.second,
        'sigma': 0.39 / unit.second,
    },
    '3j_hn_co': {
        'dihedral': 'phi',
        'delta':  172.49,
        'A':  5.34 / unit.second,
        'B': -1.46 / unit.second,
        'C': -0.29 / unit.second,
        'sigma': 0.59 / unit.second,
    },
    '3j_hn_ha': {
        'dihedral': 'phi',
        'delta':  -64.51,
        'A':  9.14 / unit.second,
        'B': -2.28 / unit.second,
        'C': -0.29 / unit.second,
        'sigma': 0.91 / unit.second,
    },
}

