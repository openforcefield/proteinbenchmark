"""List of benchmark targets specifying observables and thermodynamic state."""
from openmm import unit
from pathlib import Path
from proteinbenchmark.utilities import package_data_directory


observable_directory = Path(package_data_directory, 'observables')
pdb_directory = Path(package_data_directory, 'pdbs')

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
    'aaqaa3': {
        'target_type': 'disordered',
        'aa_sequence': 'AAQAAAAQAAAAQAA',
        'nterm_cap': 'ace',
        'cterm_cap': 'nh2',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 274.0 * unit.kelvin,
        'ph': 7.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            'fraction_helix': {
                'experimental_datasets': 'shalongo_jacs_1994',
                'observable_path': Path(
                    observable_directory, 'aaqaa', 'aaqaa_fraction_helix.dat'
                ),
            },
        }
    },
    'aaqaa3-helix': {
        'target_type': 'disordered',
        'build_method': 'helical',
        'aa_sequence': 'AAQAAAAQAAAAQAA',
        'nterm_cap': 'ace',
        'cterm_cap': 'nh2',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 274.0 * unit.kelvin,
        'ph': 7.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            'fraction_helix': {
                'experimental_datasets': 'shalongo_jacs_1994',
                'observable_path': Path(
                    observable_directory, 'aaqaa', 'aaqaa_fraction_helix.dat'
                ),
            },
        }
    },
    'ala3': {
        'target_type': 'peptide',
        'aa_sequence': 'AAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 2.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'ala3', 'ala3_scalar_couplings.dat'
                ),
            },
        }
    },
    'ala3-neutral': {
        'target_type': 'peptide',
        'aa_sequence': 'AAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 7.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'ala3', 'ala3_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'ala4', 'ala4_scalar_couplings.dat'
                ),
            },
        }
    },
    'ala4-neutral': {
        'target_type': 'peptide',
        'aa_sequence': 'AAAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 7.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'ala4', 'ala4_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'ala5', 'ala5_scalar_couplings.dat'
                ),
            },
        }
    },
    'ala5-neutral': {
        'target_type': 'peptide',
        'aa_sequence': 'AAAAA',
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 300.0 * unit.kelvin,
        'ph': 7.0,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'ala5', 'ala5_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'gag', 'gag_scalar_couplings.dat'
                ),
            },
        }
    },
    # Differences GB3 to GB1: Q2T V6I I7L K19E E24A A29V V42E
    'gb1': {
        'target_type': 'folded',
        'initial_pdb': Path(pdb_directory, 'gb1-1PGB.pdb'),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 5.6,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_n_co': 'cornilescu_jacs_1999',
        }
    },
    'gb3': {
        'target_type': 'folded',
        'initial_pdb': Path(pdb_directory, 'gb3-1P7E.pdb'),
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
        'initial_pdb': Path(pdb_directory, 'gb3-1P7E.pdb'),
        'pressure': 1.0 * unit.atmosphere,
        'temperature': 298.0 * unit.kelvin,
        'ph': 5.6,
        'ionic_strength': 0.0 * unit.molar,
        'observables': {
            '3j_ha_hb': 'miclet_jbnmr_2005',
        }
    },
    'gb3-3j-nc-cg': {
        'target_type': 'folded',
        'initial_pdb': Path(pdb_directory, 'gb3-1P7E.pdb'),
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
        'initial_pdb': Path(pdb_directory, 'gb3-1P7E.pdb'),
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'geg', 'geg_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'gfg', 'gfg_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'gkg', 'gkg_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'glg', 'glg_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'gly3', 'gly3_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'gmg', 'gmg_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'gsg', 'gsg_scalar_couplings.dat'
                ),
            },
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
            'scalar_couplings': {
                'experimental_datasets': 'hagarman_jacs_2010',
                'observable_path': Path(
                    observable_directory, 'gvg', 'gvg_scalar_couplings.dat'
                ),
            },
        }
    },
    'hewl': {
        'target_type': 'folded',
        'initial_pdb': Path(pdb_directory, 'hewl-1E8L-model-1.pdb'),
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
        'initial_pdb': Path(pdb_directory, 'hewl-1E8L-model-1.pdb'),
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
        'initial_pdb': Path(pdb_directory, 'hewl-1E8L-model-1.pdb'),
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
            'scalar_couplings': {
                'experimental_datasets': 'graf_jacs_2007',
                'observable_path': Path(
                    observable_directory, 'val3', 'val3_scalar_couplings.dat'
                ),
            },
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
    'shalongo_jacs_1994': {
        'references': [
            'Shalongo W, Dugad L, Stellwagen E. (1994). J. Am. Chem. Soc. 116, '
                '8288-8293.',
        ],
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

