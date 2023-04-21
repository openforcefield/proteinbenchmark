"""Default parameter values for analysis of trajectories."""
from openmm import unit
import pandas
from pathlib import Path
from proteinbenchmark.utilities import package_data_directory


DCD_TIME_TO_PICOSECONDS = 0.04888821 * unit.picoseconds

# Lists of atoms that make up named dihedrals in protein residues
DIHEDRAL_ATOMS = {
    resname: {
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
        'tau': {
            'atom_names': ['N', 'CA', 'C'],
        },
    }
    for resname in [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIE', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    ]
}

DIHEDRAL_ATOMS['ARG']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ASN']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ASP']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['CYS']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'SG']}
DIHEDRAL_ATOMS['GLN']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['GLU']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['HIE']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
DIHEDRAL_ATOMS['ILE']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG1']}
DIHEDRAL_ATOMS['LEU']['chi1'] = {'atom_names': ['N', 'CA', 'CB', 'CG']}
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
DIHEDRAL_ATOMS['ASN']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'OD1']}
DIHEDRAL_ATOMS['ASP']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'OD1']}
DIHEDRAL_ATOMS['GLN']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['GLU']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['HIE']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'ND1']}
DIHEDRAL_ATOMS['ILE']['chi2'] = {'atom_names': ['CA', 'CB', 'CG1', 'CD1']}
DIHEDRAL_ATOMS['LEU']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['LYS']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['MET']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'SD']}
DIHEDRAL_ATOMS['PHE']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['PRO']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD']}
DIHEDRAL_ATOMS['TRP']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['TYR']['chi2'] = {'atom_names': ['CA', 'CB', 'CG', 'CD1']}
DIHEDRAL_ATOMS['ARG']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'NE']}
DIHEDRAL_ATOMS['GLN']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'OE1']}
DIHEDRAL_ATOMS['GLU']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'OE1']}
DIHEDRAL_ATOMS['LYS']['chi3'] = {'atom_names': ['CB', 'CG', 'CD', 'CE']}
DIHEDRAL_ATOMS['MET']['chi3'] = {'atom_names': ['CB', 'CG', 'SD', 'CE']}
DIHEDRAL_ATOMS['ARG']['chi4'] = {'atom_names': ['CG', 'CD', 'NE', 'CZ']}
DIHEDRAL_ATOMS['LYS']['chi4'] = {'atom_names': ['CG', 'CD', 'CE', 'NZ']}
DIHEDRAL_ATOMS['ARG']['chi5'] = {'atom_names': ['CD', 'NE', 'CZ', 'NH1']}

DIHEDRAL_ATOMS['ASH'] = DIHEDRAL_ATOMS['ASP']
DIHEDRAL_ATOMS['CYX'] = DIHEDRAL_ATOMS['CYS']
DIHEDRAL_ATOMS['GLH'] = DIHEDRAL_ATOMS['GLU']
DIHEDRAL_ATOMS['HID'] = DIHEDRAL_ATOMS['HIE']
DIHEDRAL_ATOMS['HIP'] = DIHEDRAL_ATOMS['HIE']
DIHEDRAL_ATOMS['LYN'] = DIHEDRAL_ATOMS['LYS']

# Clusters on Ramachandran map from
# Hollingsworth SA, Karplus PA. (2010). BioMol Concepts 1, 271-283.
# The alpha cluster is a subset of the delta cluster. Alpha must go after delta.
HOLLINGSWORTH_RAMACHANDRAN_CLUSTERS = [
    {'cluster': r'$\beta$', 'phi': [180, 270], 'psi': [105, 195]},
    {'cluster': r'$\gamma$', 'phi': [60, 105], 'psi': [-90, -30]},
    {'cluster': r'$\delta$', 'phi': [225, 315], 'psi': [-60, 45]},
    {'cluster': r'$\alpha$', 'phi': [285, 315], 'psi': [-60, -30]},
    {'cluster': r'$\varepsilon$', 'phi': [45, 180], 'psi': [120, 240]},
    {'cluster': r'$\zeta$', 'phi': [195, 255], 'psi': [45, 105]},
    {'cluster': r'$P_{II}$', 'phi': [270, 315], 'psi': [120, 195]},
    {'cluster': r"$\gamma'$", 'phi': [255, 300], 'psi': [45, 105]},
    {'cluster': r"$\delta'$", 'phi': [30, 120], 'psi': [-15, 75]},
    {'cluster': r"$P_{II}'$", 'phi': [45, 120], 'psi': [165, 240]},
]

# Rotamer library from
# Hintze BJ, Lewis SM, Richardson JS, Richardson DC. (2016). Proteins 84,
#     1177-1189.
rotamer_library_dir = Path(package_data_directory, 'rotamer-library')
HINTZE_ROTAMER_LIBRARY = dict()

for resname in [
    'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
    'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
]:

    residue_df = pandas.read_csv(Path(rotamer_library_dir, f'{resname}.csv'))

    if resname == 'HIS':
        resname = 'HIE'

    dihedral_columns = [col for col in residue_df.columns if 'chi' in col]
    residue_df = residue_df[['rotamer', 'frequency%', *dihedral_columns]]
    residue_df = residue_df.sort_values(by = 'frequency%', ignore_index = True)

    HINTZE_ROTAMER_LIBRARY[resname] = residue_df

HINTZE_ROTAMER_LIBRARY['ASH'] = HINTZE_ROTAMER_LIBRARY['ASP']
HINTZE_ROTAMER_LIBRARY['CYX'] = HINTZE_ROTAMER_LIBRARY['CYS']
HINTZE_ROTAMER_LIBRARY['GLH'] = HINTZE_ROTAMER_LIBRARY['GLU']
HINTZE_ROTAMER_LIBRARY['HID'] = HINTZE_ROTAMER_LIBRARY['HIE']
HINTZE_ROTAMER_LIBRARY['HIP'] = HINTZE_ROTAMER_LIBRARY['HIE']
HINTZE_ROTAMER_LIBRARY['LYN'] = HINTZE_ROTAMER_LIBRARY['LYS']

# Karplus parameters for backbone scalar couplings from
# Wirmer J, Schwalbe H. (2002). J. Biomol. NMR 23, 47-55.
WIRMER_KARPLUS_PARAMETERS = {
    '1j_n_ca': {
        'dihedral': 'psi',
        'delta':    0.0,
        'A':  1.70 / unit.second,
        'B': -0.98 / unit.second,
        'C':  9.51 / unit.second,
        'sigma': 0.59 / unit.second,
    },
}

# Karplus parameters for backbone scalar couplings from
# Ding K, Gronenborn AM. (2004). J. Am. Chem. Soc. 126, 6232-6233.
DING_KARPLUS_PARAMETERS = {
    '2j_n_ca': {
        'dihedral': 'prev_psi',
        'delta':    0.0,
        'A': -0.66 / unit.second,
        'B': -1.52 / unit.second,
        'C':  7.85 / unit.second,
        'sigma': 0.50 / unit.second,
    },
}

# Karplus parameters for backbone scalar couplings from
# Hennig M, Bermel W, Schwalbe H, Griesinger C. (2000). J. Am. Chem. Soc. 122,
#     6268-6277.
HENNIG_KARPLUS_PARAMETERS = {
    '3j_hn_ca': {
        'dihedral': 'phi,prev_psi',
        'cos_phi': -0.23 / unit.second,
        'cos_psi': -0.20 / unit.second,
        'sin_phi':  0.07 / unit.second,
        'sin_psi':  0.08 / unit.second,
        'cos_phi_cos_psi':  0.07 / unit.second,
        'cos_phi_sin_psi':  0.12 / unit.second,
        'sin_phi_cos_psi': -0.08 / unit.second,
        'sin_phi_sin_psi': -0.14 / unit.second,
        'C':   0.54 / unit.second,
        'sigma': 0.10 / unit.second,
    },
}

# Karplus parameters for backbone scalar couplings from
VOGELI_KARPLUS_PARAMETERS = {
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

# Karplus parameters for backbone scalar couplings from
# Hu JS, Bax AJ (1997). J. Am. Chem. Soc. 119, 6360-6368.
HU_KARPLUS_PARAMETERS = {
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
CASE_DFT1_KARPLUS_PARAMETERS = {
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

# Karplus parameters for backbone scalar couplings of Ala-Ala-NH2 from
# Case DA, Scheurer C, Bruschweiler R. (2000). J. Am. Chem. Soc. 122,
#     10390-10397.
CASE_DFT2_KARPLUS_PARAMETERS = {
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

# Karplus parameters for sidechain scalar couplings from
# Perez C, Lohr F, Ruterjans H, Schmidt JM. (2001). J. Am. Chem. Soc. 123,
#     7081-7093.
PEREZ_KARPLUS_RESIDUE_MAP = {
    residue: 'ARG,ASN,ASP,GLN,GLU,HIS,LEU,LYS,MET,PHE,PRO,TRP,TYR'
    for residue in [
        'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'HID', 'HIE', 'HIP',
        'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'TRP', 'TYR',
    ]
}

PEREZ_KARPLUS_PARAMETERS = {
    '3j_ha_hb2': {
        'dihedral': 'chi1',
        'ARG,ASN,ASP,GLN,GLU,HIS,LEU,LYS,MET,PHE,PRO,TRP,TYR': {
            'delta': -120.0,
            'A':  7.23 / unit.second,
            'B': -1.37 / unit.second,
            'C':  2.40 / unit.second,
            'sigma': 0.31 / unit.second,
        },
    },
    '3j_ha_hb3': {
        'dihedral': 'chi1',
        'ARG,ASN,ASP,GLN,GLU,HIS,LEU,LYS,MET,PHE,PRO,TRP,TYR': {
            'delta':    0.0,
            'A':  7.23 / unit.second,
            'B': -1.37 / unit.second,
            'C':  2.40 / unit.second,
            'sigma': 0.31 / unit.second,
        },
    },
}

# Karplus parameters for sidechain scalar couplings from
# Chou JJ, Case DA, Bax A. (2003). J. Am. Chem. Soc. 125, 8959-8966.
CHOU_KARPLUS_PARAMETERS = {
    '3j_co_cg1': {
        'dihedral': 'chi1',
        'VAL': {
            'delta': -115.0,
            'A':  3.42 / unit.second,
            'B': -0.59 / unit.second,
            'C':  0.17 / unit.second,
            'sigma': 0.25 / unit.second,
        },
    },
    '3j_co_cg2': {
        'dihedral': 'chi1',
        'ILE': {
            'delta':  125.0,
            'A':  3.42 / unit.second,
            'B': -0.59 / unit.second,
            'C':  0.17 / unit.second,
            'sigma': 0.25 / unit.second,
        },
        'THR': {
            'delta':  137.0,
            'A':  2.76 / unit.second,
            'B': -0.67 / unit.second,
            'C':  0.19 / unit.second,
            'sigma': 0.21 / unit.second,
        },
        'VAL': {
            'delta':    5.0,
            'A':  3.42 / unit.second,
            'B': -0.59 / unit.second,
            'C':  0.17 / unit.second,
            'sigma': 0.25 / unit.second,
        },
    },
    '3j_n_cg1': {
        'dihedral': 'chi1',
        'VAL': {
            'delta':    6.0,
            'A':  2.64 / unit.second,
            'B':  0.26 / unit.second,
            'C': -0.22 / unit.second,
            'sigma': 0.25 / unit.second,
        },
    },
    '3j_n_cg2': {
        'dihedral': 'chi1',
        'ILE': {
            'delta': -114.0,
            'A':  2.64 / unit.second,
            'B':  0.26 / unit.second,
            'C': -0.22 / unit.second,
            'sigma': 0.25 / unit.second,
        },
        'THR': {
            'delta': -113.0,
            'A':  2.01 / unit.second,
            'B':  0.21 / unit.second,
            'C': -0.12 / unit.second,
            'sigma': 0.21 / unit.second,
        },
        'VAL': {
            'delta':  126.0,
            'A':  2.64 / unit.second,
            'B':  0.26 / unit.second,
            'C': -0.22 / unit.second,
            'sigma': 0.25 / unit.second,
        },
    },
}

# 3j_n_co trans hydrogen bond scalar couplings
# Eq 12 from Barfield M. (2002). J. Am. Chem. Soc. 124, 4158-4168.
# Eq 2 from Sass HJ, Schmid FFF, Grzesiek S. (2007). J. Am. Chem. Soc. 129,
#     5898-5903.
BARFIELD_KARPLUS_PARAMETERS = {
    '3j_n_co': {
        'min_dist': 1.760 * unit.angstrom,
        'exponent': 3.2 / unit.angstrom,
        'A': 0.62 / unit.second,
        'B': 0.92 / unit.second,
        'C': 0.14 / unit.second,
        'D': -1.31 / unit.second,
        'sigma': 0.12 / unit.second,
    },
}

