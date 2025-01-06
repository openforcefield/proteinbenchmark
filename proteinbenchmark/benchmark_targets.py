"""List of benchmark targets specifying observables and thermodynamic state."""
from pathlib import Path

from openff.units import unit

from proteinbenchmark.utilities import package_data_directory

observable_directory = Path(package_data_directory, "observables")
pdb_directory = Path(package_data_directory, "pdbs")

# List of benchmark targets with observables and thermodynamic state (pressure,
# temperature, pH, ionic strength)
benchmark_targets = {
    "test": {
        "target_type": "test",
        "aa_sequence": "AAAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 7.4,
        "ionic_strength": 0.0 * unit.molar,
        "equil_traj_length": 1.0 * unit.picosecond,
        "equil_frame_length": 100.0 * unit.femtosecond,
        "traj_length": 5.0 * unit.picosecond,
        "frame_length": 100.0 * unit.femtosecond,
        "checkpoint_length": 1.0 * unit.picosecond,
        "save_state_length": 1.0 * unit.picosecond,
    },
    "aaqaa3": {
        "target_type": "disordered",
        "aa_sequence": "AAQAAAAQAAAAQAA",
        "nterm_cap": "ace",
        "cterm_cap": "nh2",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 274.0 * unit.kelvin,
        "ph": 7.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "fraction_helix": {
                "experimental_datasets": ["shalongo_jacs_1994"],
                "observable_path": Path(
                    observable_directory,
                    "aaqaa3",
                    "aaqaa3-fraction-helix-by-residue.dat",
                ),
            },
        },
    },
    "aaqaa3-helix": {
        "target_type": "disordered",
        "build_method": "alpha",
        "aa_sequence": "AAQAAAAQAAAAQAA",
        "nterm_cap": "ace",
        "cterm_cap": "nh2",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 274.0 * unit.kelvin,
        "ph": 7.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "fraction_helix": {
                "experimental_datasets": ["shalongo_jacs_1994"],
                "observable_path": Path(
                    observable_directory,
                    "aaqaa3",
                    "aaqaa3-fraction-helix-by-residue.dat",
                ),
            },
        },
    },
    "ala3": {
        "target_type": "peptide",
        "aa_sequence": "AAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala3", "ala3-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala3-neutral": {
        "target_type": "peptide",
        "aa_sequence": "AAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 7.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala3", "ala3-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala4": {
        "target_type": "peptide",
        "aa_sequence": "AAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala4", "ala4-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala4-neutral": {
        "target_type": "peptide",
        "aa_sequence": "AAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 7.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala4", "ala4-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala5": {
        "target_type": "peptide",
        "aa_sequence": "AAAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala5", "ala5-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala5-hmr": {
        "target_type": "peptide",
        "aa_sequence": "AAAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala5", "ala5-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala5-neutral": {
        "target_type": "peptide",
        "aa_sequence": "AAAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 7.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala5", "ala5-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala6": {
        "target_type": "peptide",
        "aa_sequence": "AAAAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala6", "ala6-scalar-couplings.dat"
                ),
            },
        },
    },
    "ala7": {
        "target_type": "peptide",
        "aa_sequence": "AAAAAAA",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "ala7", "ala7-scalar-couplings.dat"
                ),
            },
        },
    },
    "asyn": {
        "target_type": "disordered",
        "aa_sequence": (
            "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVH"
            "GVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQL"
            "GKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
        ),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 288.0 * unit.kelvin,
        "ph": 7.4,
        "ionic_strength": 0.075 * unit.molar,
        "observables": {
            "chemical_shifts": "rao_jmb_2009",
        },
    },
    "asyn-3j": {
        "target_type": "disordered",
        "aa_sequence": (
            "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVH"
            "GVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQL"
            "GKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
        ),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 288.0 * unit.kelvin,
        "ph": 6.0,
        "ionic_strength": 0.050 * unit.molar,
        "observables": {
            "3j_co_co": "lee_jacs_2015",
            "3j_ha_hb": "lee_jacs_2015",
        },
    },
    "bpti": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "bpti-1PIT-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 305.0 * unit.kelvin,
        "ph": 6.2,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "backbone_rdc": "moglich_jbnmr_2002",
        },
    },
    "bpti-3j-hn-ha": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "bpti-1PIT-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 309.0 * unit.kelvin,
        "ph": 3.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "3j_hn_ha": "pardi_jmb_1984",
        },
    },
    "bpti-3j-ha-hb": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "bpti-1PIT-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 309.0 * unit.kelvin,
        "ph": 4.6,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "pardi_jmb_1984",
                    "berndt_jmb_1992",
                    "lindorff-larsen_prot_2010",
                ],
                "observable_path": Path(
                    observable_directory, "bpti", "bpti-scalar-couplings.dat"
                ),
            },
        },
    },
    "bpti-t1-relaxation": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "bpti-1PIT-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 313.0 * unit.kelvin,
        "ph": 5.1,
        "ionic_strength": 0.125 * unit.molar,
        "observables": {
            "t1_relaxation": "balasubramanian_jmr_1994",
        },
    },
    "cln025": {
        "target_type": "peptide",
        "aa_sequence": "YYDPETGTWY",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 313.0 * unit.kelvin,
        "ph": 7.0,
        "ionic_strength": 0.0 * unit.molar,
        "traj_length": 10.0 * unit.microsecond,
        "observables": {
            "fraction_folded": {
                "experimental_datasets": ["honda_jacs_2008"],
                "observable_path": Path(
                    observable_directory,
                    "cln025",
                    "cln025-fraction-folded-by-temperature.dat",
                ),
            },
        },
    },
    "gag": {
        "target_type": "peptide",
        "aa_sequence": "GAG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "gag", "gag-scalar-couplings.dat"
                ),
            },
        },
    },
    # Differences GB3 to GB1: Q2T V6I I7L K19E E24A A29V V42E
    "gb1": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb1-1PGB.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 5.6,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "3j_n_co": "cornilescu_jacs_1999_b",
        },
    },
    "gb3": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "chou_jacs_2003",
                    "miclet_jbnmr_2005",
                    "vogeli_jacs_2007",
                ],
                "observable_path": Path(
                    observable_directory, "gb3", "gb3-scalar-couplings.dat"
                ),
            },
            "h_bond_scalar_couplings": {
                "experimental_datasets": ["cornilescu_jacs_1999_b"],
                "observable_path": Path(observable_directory, "gb3", "gb3-3j-n-co.dat"),
            },
            "backbone_rdc": "ulmer_jacs_2003",
        },
    },
    "gb3-test-5us": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "chou_jacs_2003",
                    "miclet_jbnmr_2005",
                    "vogeli_jacs_2007",
                ],
                "observable_path": Path(
                    observable_directory, "gb3", "gb3-scalar-couplings.dat"
                ),
            },
            "h_bond_scalar_couplings": {
                "experimental_datasets": ["cornilescu_jacs_1999_b"],
                "observable_path": Path(observable_directory, "gb3", "gb3-3j-n-co.dat"),
            },
            "backbone_rdc": "ulmer_jacs_2003",
        },
        "traj_length": 5 * unit.microsecond,
    },
    "gb3-test": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "chou_jacs_2003",
                    "miclet_jbnmr_2005",
                    "vogeli_jacs_2007",
                ],
                "observable_path": Path(
                    observable_directory, "gb3", "gb3-scalar-couplings.dat"
                ),
            },
            "h_bond_scalar_couplings": {
                "experimental_datasets": "cornilescu_jacs_1999",
                "observable_path": Path(observable_directory, "gb3", "gb3-3j-n-co.dat"),
            },
            "backbone_rdc": "ulmer_jacs_2003",
        },
        "traj_length": 10 * unit.nanosecond,
        "equil_timestep": 0.002 * unit.picosecond,
    },
    "gb3-debug": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "chou_jacs_2003",
                    "miclet_jbnmr_2005",
                    "vogeli_jacs_2007",
                ],
                "observable_path": Path(
                    observable_directory, "gb3", "gb3-scalar-couplings.dat"
                ),
            },
            "h_bond_scalar_couplings": {
                "experimental_datasets": "cornilescu_jacs_1999",
                "observable_path": Path(observable_directory, "gb3", "gb3-3j-n-co.dat"),
            },
            "backbone_rdc": "ulmer_jacs_2003",
        },
        "traj_length": 5 * unit.microsecond,
        "checkpoint_length": 4.0 * unit.femtosecond,
        "save_state_length": 4.0 * unit.femtosecond,
    },
    "gb3-3j-ha-hb": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 5.6,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "3j_ha_hb": "miclet_jbnmr_2005",
        },
    },
    "gb3-3j-nc-cg": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.075 * unit.molar,
        "observables": {
            "3j_co_cg": "chou_jacs_2003",
            "3j_n_cg": "chou_jacs_2003",
        },
    },
    "gb3-backbone-S2": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.050 * unit.molar,
        "observables": {
            "backbone_S2": "yao_jacs_2010",
        },
    },
    "geg": {
        "target_type": "peptide",
        "aa_sequence": "GEG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "geg", "geg-scalar-couplings.dat"
                ),
            },
        },
    },
    "gfg": {
        "target_type": "peptide",
        "aa_sequence": "GFG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "gfg", "gfg-scalar-couplings.dat"
                ),
            },
        },
    },
    "gkg": {
        "target_type": "peptide",
        "aa_sequence": "GKG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 1.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "gkg", "gkg-scalar-couplings.dat"
                ),
            },
        },
    },
    "glg": {
        "target_type": "peptide",
        "aa_sequence": "GLG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "glg", "glg-scalar-couplings.dat"
                ),
            },
        },
    },
    "gly3": {
        "target_type": "peptide",
        "aa_sequence": "GGG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "gly3", "gly3-scalar-couplings.dat"
                ),
            },
        },
    },
    "gmg": {
        "target_type": "peptide",
        "aa_sequence": "GMG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 1.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "3j_hn_ha": "hagarman_jacs_2010",
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "gmg", "gmg-scalar-couplings.dat"
                ),
            },
        },
    },
    "gsg": {
        "target_type": "peptide",
        "aa_sequence": "GSG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "gsg", "gsg-scalar-couplings.dat"
                ),
            },
        },
    },
    "gvg": {
        "target_type": "peptide",
        "aa_sequence": "GVG",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["hagarman_jacs_2010"],
                "observable_path": Path(
                    observable_directory, "gvg", "gvg-scalar-couplings.dat"
                ),
            },
        },
    },
    "hewl": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "hewl-1E8L-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 308.0 * unit.kelvin,
        "ph": 3.8,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["schwalbe_prosci_2001"],
                "observable_path": Path(
                    observable_directory, "hewl", "hewl-scalar-couplings.dat"
                ),
            },
            "backbone_S2": "buck_biochem_1995",
        },
    },
    "hewl-sidechain-S2": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "hewl-1E8L-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 308.0 * unit.kelvin,
        "ph": 4.7,
        "ionic_strength": 0.100 * unit.molar,
        "observables": {
            "sidechain_S2": "moorman_prosci_2012",
        },
    },
    "hewl-rdc": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "hewl-1E8L-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 308.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.010 * unit.molar,
        "observables": {
            "backbone_rdc": "schwalbe_prosci_2001",
        },
    },
    "ubq": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 303.0 * unit.kelvin,
        "ph": 4.7,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "wang_jacs_1996",
                    "hu_jacs_1997",
                    "chou_jacs_2003",
                    "lindorff-larsen_prot_2010",
                ],
                "observable_path": Path(
                    observable_directory, "ubq", "ubq-scalar-couplings.dat"
                ),
            },
            "h_bond_scalar_couplings": {
                "experimental_datasets": ["cordier_jacs_1999"],
                "observable_path": Path(
                    observable_directory, "ubq", "ubq-3j-n-co-cordier.dat"
                ),
            },
        },
    },
    "ubq-3j-hn-ha": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 303.0 * unit.kelvin,
        "ph": 4.7,
        "ionic_strength": 0.010 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["wang_jacs_1996"],
            },
        },
    },
    "ubq-3j-n-co-cordier": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 318.0 * unit.kelvin,
        "ph": 4.6,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "h_bond_scalar_couplings": {
                "experimental_datasets": ["cordier_jacs_1999"],
                "observable_path": Path(
                    observable_directory, "ubq", "ubq-3j-n-co-cordier.dat"
                ),
            },
        },
    },
    "ubq-3j-n-co-cornilescu": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "h_bond_scalar_couplings": {
                "experimental_datasets": ["cornilescu_jacs_1999_a"],
                "observable_path": Path(
                    observable_directory, "ubq", "ubq-3j-n-co-cornilescu.dat"
                ),
            },
        },
    },
    "ubq-3j-nc-cg": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "gb3-1P7E.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 298.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.075 * unit.molar,
        "observables": {
            "3j_co_cg": "chou_jacs_2003",
            "3j_n_cg": "chou_jacs_2003",
            "3j_ha_hb": "lindorff-larsen_prot_2010",
        },
    },
    "ubq-rdc-backbone": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 303.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "backbone_rdc": {
                "experimental_datasets": ["hus_jacs_2003"],
            },
        },
    },
    "ubq-rdc-sidechain": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 303.0 * unit.kelvin,
        "ph": 6.5,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "sidechain_rdc": {
                "experimental_datasets": ["ottinger_jacs_1998"],
            },
        },
    },
    "ubq-s2-backbone": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 4.7,
        "ionic_strength": 0.010 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["tjandra_jacs_1995"],
            },
        },
    },
    "ubq-s2-sidechain": {
        "target_type": "folded",
        "initial_pdb": Path(pdb_directory, "ubq-1D3Z-model-1.pdb"),
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 303.0 * unit.kelvin,
        "ph": 5.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": [
                    "lee_jacs_1999",
                ],
            },
        },
    },
    "val3": {
        "target_type": "peptide",
        "aa_sequence": "VVV",
        "pressure": 1.0 * unit.atmosphere,
        "temperature": 300.0 * unit.kelvin,
        "ph": 2.0,
        "ionic_strength": 0.0 * unit.molar,
        "observables": {
            "scalar_couplings": {
                "experimental_datasets": ["graf_jacs_2007"],
                "observable_path": Path(
                    observable_directory, "val3", "val3-scalar-couplings.dat"
                ),
            },
        },
    },
}

# References for experimental datasets
experimental_datasets = {
    "berndt_jmb_1992": {
        "references": [
            "Berndt KD, Guntert P, Orbons LP, Wuthrich K. (1992). J. Mol. Biol. "
            "227, 757-775.",
        ],
    },
    "balasubramanian_jmr_1994": {
        "references": [
            "Balasubramanian S, Nirmala R, Beveridge DL, Bolton PH. (1994). "
            "J. Magn. Reson. 104, 240-249.",
        ],
    },
    "buck_biochem_1995": {
        "references": [
            "Buck M, Boyd J, Redfield C, MacKenzie DA, Jeenes DJ, Archer DB, "
            "Dobson CM. (1995). Biochem. 34, 4041-4055.",
        ],
    },
    "chou_jacs_2003": {
        "references": [
            "Chou JJ, Case DA, Bax A. (2003). J. Am .Chem. Soc. 125, 8959-8966.",
        ],
    },
    "cordier_jacs_1999": {
        "references": [
            "Cordier F, Grzesiek S. (1999). J. Am. Chem. Soc. 121, 1601-1602.",
        ],
    },
    "cornilescu_jacs_1999_a": {
        "references": [
            "Cornilescu G, Hu JS, Bax A. (1999). J. Am. Chem. Soc. 121, 2949-2950.",
        ],
    },
    "cornilescu_jacs_1999_b": {
        "references": [
            "Cornilescu G, Ramirez BE, Frank MK, Clore GM, Gronenborn AM, "
            "Bax A. (1999). J. Am. Chem. Soc. 121, 6275-6279.",
        ],
    },
    "graf_jacs_2007": {
        "references": [
            "Graf J, Nguyen PH, Stock G, Schwalbe H. (2007). J. Am. Chem. Soc. "
            "129, 1179-1189.",
        ],
    },
    "hagarman_jacs_2010": {
        "references": [
            "Hagarman A, Measey TJ, Mathieu D, Schwalbe H, "
            "Schweitzer-Stenner R. (2010). J. Am. Chem. Soc. 132, 540-551.",
        ],
    },
    "hu_jacs_1997": {
        "references": [
            "Hu JS, Bax A (1997). J. Am. Chem. Soc. 119, 6360-6368.",
        ],
    },
    "hus_jacs_2003": {
        "references": [
            "Hus JC, Peti W, Griesinger C, Bruschweiler R. (2003). J. Am. "
            "Chem. Soc. 125, 5596-5597.",
        ],
    },
    "lee_jacs_1999": {
        "references": [
            "Lee AL, Flynn PF, Wand AJ. (1999). J. Am. Chem. Soc. 121, 2891-2902.",
        ],
    },
    "lee_jacs_2015": {
        "references": [
            "Lee JH, Li F, Grishaev A, Bax A. (2015). J. Am. Chem. Soc. 137, "
            "1432-1435.",
        ],
    },
    "lindorff-larsen_prot_2010": {
        "references": [
            "Lindorff-Larsen K, Piana S, Palmo K, Maragakis P, Klepeis JL, "
            "Dror RO, Shaw DE. (2010). Proteins 78, 1950-1958.",
        ],
    },
    "miclet_jbnmr_2005": {
        "references": [
            "Miclet E, Boisbouvier J, Bax A. (2005). J. Biomol. NMR 31, 201-216.",
        ],
    },
    "moglich_jbnmr_2002": {
        "references": [
            "Moglich A, Wenzler M, Kramer F, Glaser SJ, Brunner E. (2002). "
            "J. Biomol. NMR 23, 211-219.",
        ],
    },
    "moorman_prosci_2012": {
        "references": [
            "Moorman V, Valentine K, Wand JA. (2012). Protein Sci. 21, 1066-1073.",
        ],
        "dataset_ids": ["bmrb 18304"],
    },
    "ottinger_jacs_1998": {
        "references": [
            "Ottinger M, Bax A. (1998). J. Am. Chem. Soc. 120, 12334-12341.",
        ],
    },
    "pardi_jmb_1984": {
        "references": [
            "Pardi A, Billeter M, Wuthrich K. (1984). J. Mol. Biol. 180, 741-751.",
        ],
    },
    "rao_jmb_2009": {
        "references": [
            "Rao JN, Kim Y, Park L, Ulmer T. (2009). J. Mol. Biol. 390, 516-529.",
        ],
        "dataset_ids": ["bmrb 16300"],
    },
    "schwalbe_prosci_2001": {
        "references": [
            "Schwalbe H, Grimshaw SB, Spencer A, Buck M, Boyd J, Dobson CM, "
            "Redfield C, Smith LJ. (2001). Protein Sci. 10, 677-688.",
        ],
        "dataset_ids": ["bmrb 4831", "pdb 1E8L"],
    },
    "shalongo_jacs_1994": {
        "references": [
            "Shalongo W, Dugad L, Stellwagen E. (1994). J. Am. Chem. Soc. 116, "
            "8288-8293.",
        ],
    },
    "tjandra_jacs_1995": {
        "references": [
            "Tjandra N, Feller SE, Pastor RW, Bax A. (1995). J. Am. Chem. Soc. "
            "117, 12562-12566.",
        ],
    },
    "ulmer_jacs_2003": {
        "references": [
            "Ulmer TS, Ramirez BE, Delaglio F, Bax A. (2003). J. Am. Chem. "
            "Soc. 125, 9179-9191.",
        ],
        "dataset_ids": ["pdb 1P7E"],
    },
    "vogeli_jacs_2007": {
        "references": [
            "Vogeli B, Ying J, Grishaev A, Bax A. (2007). J. Am. Chem. Soc. "
            "129, 9377-9385.",
        ],
    },
    "wang_jacs_1996": {
        "references": [
            "Wang AC, Bax A. (1996). J. Am. Chem. Soc. 118, 2483-2494.",
        ],
    },
    "yao_jacs_2010": {
        "references": [
            "Yao L, Grishaev A, Cornilescu G, Bax A. (2010). J. Am. Chem. Soc. "
            "132, 4295-4309.",
        ],
    },
}
