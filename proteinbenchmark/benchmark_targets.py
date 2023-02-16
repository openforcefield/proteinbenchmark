"""List of benchmark targets specifying observables and thermodynamic state."""
from openmm import unit
from pathlib import Path
from proteinbenchmark.utilities import package_data_directory


# List of benchmark targets with observables and thermodynamic state (pressure,
# temperature, pH, ionic strength)
benchmark_targets = {
    'test': {
        'target_type': 'test',
        'aa_sequence': 'AAAAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 7.4,
        'ionic_strength': 0.0 * unit.molar,
        'equil_traj_length': 1.0 * unit.picosecond,
        'equil_frame_length': 100.0 * unit.femtosecond,
        'traj_length': 5.0 * unit.picosecond,
        'frame_length': 100.0 * unit.femtosecond,
        'checkpoint_length': 1.0 * unit.picosecond,
        'save_state_length': 1.0 * unit.picosecond,
    },
    'ala3': {
        'target_type': 'peptide',
        'aa_sequence': 'AAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'graf_jacs_2007',
            '2j_n_ca': 'graf_jacs_2007',
            '3j_co_co': 'graf_jacs_2007',
            '3j_ha_co': 'graf_jacs_2007',
            '3j_hn_ca': 'graf_jacs_2007',
            '3j_hn_cb': 'graf_jacs_2007',
            '3j_hn_co': 'graf_jacs_2007',
            '3j_hn_ha': 'graf_jacs_2007',
        }
    },
    'ala4': {
        'target_type': 'peptide',
        'aa_sequence': 'AAAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'graf_jacs_2007',
            '2j_n_ca': 'graf_jacs_2007',
            '3j_co_co': 'graf_jacs_2007',
            '3j_ha_co': 'graf_jacs_2007',
            '3j_hn_ca': 'graf_jacs_2007',
            '3j_hn_cb': 'graf_jacs_2007',
            '3j_hn_co': 'graf_jacs_2007',
            '3j_hn_ha': 'graf_jacs_2007',
        }
    },
    'ala5': {
        'target_type': 'peptide',
        'aa_sequence': 'AAAAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'graf_jacs_2007',
            '2j_n_ca': 'graf_jacs_2007',
            '3j_co_co': 'graf_jacs_2007',
            '3j_ha_co': 'graf_jacs_2007',
            '3j_hn_ca': 'graf_jacs_2007',
            '3j_hn_cb': 'graf_jacs_2007',
            '3j_hn_co': 'graf_jacs_2007',
            '3j_hn_ha': 'graf_jacs_2007',
        }
    },
    'asyn': {
        'target_type': 'disordered',
        'aa_sequence': (
            'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVH'
            'GVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQL'
            'GKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
        ),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 288.0 * unit.kelvin,
        'ph': 7.4,
        'ionic_strength': 0.075 * unit.molar,
        'observables': {
            'chemical_shifts': 'rao_jmb_2009',
        }
    },
    'asyn-3j': {
        'target_type': 'disordered',
        'aa_sequence': (
            'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVH'
            'GVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQL'
            'GKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
        ),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 288.0 * unit.kelvin,
        'ph': 6.0,
        'ionic_strength': 0.050 * unit.molar,
        'observables': {
            '3j_co_co': 'lee_jacs_2015',
            '3j_ha_hb': 'lee_jacs_2015',
        }
    },
    'gag': {
        'target_type': 'peptide',
        'aa_sequence': 'GAG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'hagarman_jacs_2010',
            '2j_n_ca': 'hagarman_jacs_2010',
            '3j_ha_co': 'hagarman_jacs_2010',
            '3j_hn_cb': 'hagarman_jacs_2010',
            '3j_hn_co': 'hagarman_jacs_2010',
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'gb3': {
        'target_type': 'folded',
        'initial_pdb': Path(package_data_directory, 'pdbs', 'gb3-1P7E.pdb'),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 6.5,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_hn_cb': 'vogeli_jacs_2007',
            '3j_hn_co': 'vogeli_jacs_2007',
            '3j_hn_ha': 'vogeli_jacs_2007',
            'backbone_rdc': 'ulmer_jacs_2003',
        }
    },
    'gb3-3j-ha-hb': {
        'target_type': 'folded',
        'initial_pdb': Path(package_data_directory, 'pdbs', 'gb3-1P7E.pdb'),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 5.6,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_ha_hb': 'miclet_jbnmr_2005',
            '3j_n_co': 'cornilescu_jacs_1999',
        }
    },
    'gb3-3j-nc-cg': {
        'target_type': 'folded',
        'initial_pdb': Path(package_data_directory, 'pdbs', 'gb3-1P7E.pdb'),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 6.5,
        'ionic_strength': 0.075 * unit.molar,
        'observables': {
            '3j_co_cg': 'chou_jacs_2003',
            '3j_n_cg': 'chou_jacs_2003',
        }
    },
    'gb3-backbone-S2': {
        'target_type': 'folded',
        'initial_pdb': Path(package_data_directory, 'pdbs', 'gb3-1P7E.pdb'),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 6.5,
        'ionic_strength': 0.050 * unit.molar,
        'observables': {
            'backbone_S2': 'yao_jacs_2010',
        }
    },
    'geg': {
        'target_type': 'peptide',
        'aa_sequence': 'GEG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'hagarman_jacs_2010',
            '2j_n_ca': 'hagarman_jacs_2010',
            '3j_ha_co': 'hagarman_jacs_2010',
            '3j_hn_cb': 'hagarman_jacs_2010',
            '3j_hn_co': 'hagarman_jacs_2010',
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'gfg': {
        'target_type': 'peptide',
        'aa_sequence': 'GFG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'hagarman_jacs_2010',
            '2j_n_ca': 'hagarman_jacs_2010',
            '3j_ha_co': 'hagarman_jacs_2010',
            '3j_hn_cb': 'hagarman_jacs_2010',
            '3j_hn_co': 'hagarman_jacs_2010',
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'gkg': {
        'target_type': 'peptide',
        'aa_sequence': 'GKG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 1.5,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'glg': {
        'target_type': 'peptide',
        'aa_sequence': 'GLG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'hagarman_jacs_2010',
            '2j_n_ca': 'hagarman_jacs_2010',
            '3j_ha_co': 'hagarman_jacs_2010',
            '3j_hn_cb': 'hagarman_jacs_2010',
            '3j_hn_co': 'hagarman_jacs_2010',
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'gly3': {
        'target_type': 'peptide',
        'aa_sequence': 'GGG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'graf_jacs_2007',
            '2j_n_ca': 'graf_jacs_2007',
            '3j_co_co': 'graf_jacs_2007',
            '3j_ha_co': 'graf_jacs_2007',
            '3j_hn_ca': 'graf_jacs_2007',
            '3j_hn_cb': 'graf_jacs_2007',
            '3j_hn_co': 'graf_jacs_2007',
            '3j_hn_ha': 'graf_jacs_2007',
        }
    },
    'gmg': {
        'target_type': 'peptide',
        'aa_sequence': 'GMG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 1.5,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'gsg': {
        'target_type': 'peptide',
        'aa_sequence': 'GSG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'hagarman_jacs_2010',
            '2j_n_ca': 'hagarman_jacs_2010',
            '3j_ha_co': 'hagarman_jacs_2010',
            '3j_hn_cb': 'hagarman_jacs_2010',
            '3j_hn_co': 'hagarman_jacs_2010',
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'gvg': {
        'target_type': 'peptide',
        'aa_sequence': 'GVG',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'hagarman_jacs_2010',
            '2j_n_ca': 'hagarman_jacs_2010',
            '3j_ha_co': 'hagarman_jacs_2010',
            '3j_hn_cb': 'hagarman_jacs_2010',
            '3j_hn_co': 'hagarman_jacs_2010',
            '3j_hn_ha': 'hagarman_jacs_2010',
        }
    },
    'hewl': {
        'target_type': 'folded',
        'initial_pdb': Path(
            package_data_directory, 'pdbs', 'hewl-1E8L-model-1.pdb'
        ),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 308.0 * unit.kelvin,
        'ph': 3.8,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_co_hb': 'schwalbe_prosci_2001',
            '3j_ha_hb': 'schwalbe_prosci_2001',
            'backbone_S2': 'buck_biochem_1995',
        }
    },
    'hewl-sidechain-S2': {
        'target_type': 'folded',
        'initial_pdb': Path(
            package_data_directory, 'pdbs', 'hewl-1E8L-model-1.pdb'
        ),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 308.0 * unit.kelvin,
        'ph': 4.7,
        'ionic_strength': 0.100 * unit.molar,
        'observables': {
            'sidechain_S2': 'moorman_prosci_2012',
        }
    },
    'hewl-rdc': {
        'target_type': 'folded',
        'initial_pdb': Path(
            package_data_directory, 'pdbs', 'hewl-1E8L-model-1.pdb'
        ),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 308.0 * unit.kelvin,
        'ph': 6.5,
        'ionic_strength': 0.010 * unit.molar,
        'observables': {
            'backbone_rdc': 'schwalbe_prosci_2001',
        }
    },
    'val3': {
        'target_type': 'peptide',
        'aa_sequence': 'VVV',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '1j_n_ca': 'graf_jacs_2007',
            '2j_n_ca': 'graf_jacs_2007',
            '3j_co_co': 'graf_jacs_2007',
            '3j_ha_co': 'graf_jacs_2007',
            '3j_hn_ca': 'graf_jacs_2007',
            '3j_hn_cb': 'graf_jacs_2007',
            '3j_hn_co': 'graf_jacs_2007',
            '3j_hn_ha': 'graf_jacs_2007',
        }
    },
}

# References for experimental datasets
experimental_datasets = {
    'buck_biochem_1995': {
        'references': [
            'Buck M, Boyd J, Redfield C, MacKenzie DA, Jeenes DJ, Archer DB, '
                'Dobson CM. (1995). Biochem. 34, 4041-4055.',
        ],
    },
    'chou_jacs_2003': {
        'references': [
            'Chou JJ, Case DA, Bax A. (2003). J. Am .Chem. Soc. 125, '
                '8959-8966.',
        ],
    },
    'cornilescu_jacs_1999': {
        'references': [
            'Cornilescu G, Ramirez BE, Frank MK, Clore GM, Gronenborn AM, '
                'Bax A. (1999). J. Am. Chem. Soc. 121, 6275-6279.',
        ],
    },
    'graf_jacs_2007': {
        'references': [
            'Graf J, Nguyen PH, Stock G, Schwalbe H. (2007). J. Am. Chem. Soc. '
                '129, 1179-1189.',
        ],
    },
    'hagarman_jacs_2010': {
        'references': [
            'Hagarman A, Measey TJ, Mathieu D, Schwalbe H, '
                'Schweitzer-Stenner R. (2010). J. Am. Chem. Soc. 132, 540-551.',
        ],
    },
    'lee_jacs_2015': {
        'references': [
            'Lee JH, Li F, Grishaev A, Bax A. (2015). J. Am. Chem. Soc. 137, '
                '1432-1435.',
        ],
    },
    'miclet_jbnmr_2005': {
        'references': [
            'Miclet E, Boisbouvier J, Bax A. (2005). J. Biomol. NMR 31, '
                '201-216.',
        ],
    },
    'moorman_prosci_2012': {
        'references': [
            'Moorman V, Valentine K, Wand JA. (2012). Protein Sci. 21, '
                '1066-1073.',
        ],
        'dataset_ids': ['bmrb 18304'],
    },
    'rao_jmb_2009': {
        'references': [
            'Rao JN, Kim Y, Park L, Ulmer T. (2009). J. Mol. Biol. 390, '
                '516-529.'
        ],
        'dataset_ids': ['bmrb 16300'],
    },
    'schwalbe_prosci_2001': {
        'references': [
            'Schwalbe H, Grimshaw SB, Spencer A, Buck M, Boyd J, Dobson CM, '
                'Redfield C, Smith LJ. (2001). Protein Sci. 10, 677-688.',
            'Lindorff-Larsen K, Piana S, Palmo K, Maragakis P, Klepeis JL, '
                'Dror RO, Shaw DE. (2010). Proteins 78, 1950-1958.',
        ],
        'dataset_ids': ['bmrb 4831', 'pdb 1E8L'],
    },
    'ulmer_jacs_2003': {
        'references': [
            'Ulmer TS, Ramirez BE, Delaglio F, Bax A. (2003). J. Am. Chem. '
                'Soc. 125, 9179-9191.',
        ],
        'dataset_ids': ['pdb 1P7E'],
    },
    'vogeli_jacs_2007': {
        'references': [
            'Vogeli B, Ying J, Grishaev A, Bax A. (2007). J. Am. Chem. Soc. '
                '129, 9377-9385.',
        ],
    },
    'yao_jacs_2010': {
        'references': [
            'Yao L, Grishaev A, Cornilescu G, Bax A. (2010). J. Am. Chem. Soc. '
                '132, 4295-4309.',
        ],
    },
}

