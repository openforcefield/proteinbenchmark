proteinbenchmark
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/chapincavender/proteinbenchmark/workflows/CI/badge.svg)](https://github.com/chapincavender/proteinbenchmark/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/chapincavender/proteinbenchmark/branch/main/graph/badge.svg)](https://codecov.io/gh/chapincavender/proteinbenchmark/branch/main)


This repository contains datasets of experimental observables and code for preparing, running, and analyzing benchmarks for OpenFF protein force fields.


## Installation

This repository is still under development and should not be used for production results.
To install the development version, clone this repository and then run the following in the top level directory.

```
conda env create -f devtools/conda-envs/proteinbenchmark.yaml
conda activate openff-proteinbenchmark
pip install -e
``` 


## Usage

To use this package, create an instance of the `ProteinBenchmarkSystem` class.
A benchmark system for the pentaalanine peptide using the Amber ff14SB force field and TIP3P water looks like this.

```
from proteinbenchmark import (ProteinBenchmarkSystem,
                              benchmark_targets,
                              force_fields)

benchmark_system = ProteinBenchmarkSystem(
    result_directory="results",
    target_name="ala5",
    target_parameters=benchmark_targets["ala5"],
    force_field_name="ff14sb-tip3p",
    water_model="tip3p",
    force_field_file=force_fields["ff14sb-tip3p"]["force_field_file"]
    water_model_file=force_fields["ff14sb-tip3p"]["water_model_file"]
)
```

`result_directory` is the directory path to store results.
Output will be written to the subdirectory `{result_directory}/{target_name}-{force_field_name}`.

Once constructed, set up the benchmark system by calling `benchmark_system.setup()`.
This function will:
- Build initial coordinates according to the entries in the `target_parameters` dictionary (see below)
- Solvate the system in a rhombic dodecahedron of water and add Na+ and Cl- ions
- Parametrize the solvated system with the force field and water model
- Write an OpenMM system to an XML file
- Perform an energy minimization

The output of the `setup()` function will be written to `{result_directory}/{target_name}-{force_field_name}/setup`.

If the system setup finished correctly, run equilibration and production simulations for by calling `benchmark_system.run_simulations()`.
Additional replicas can be run by passing an integer to the `replica` keyword argument, e.g. `benchmark_system.run_simulations(replica=2)`.
The output of the `run_simulations()` function will be written to `{result_directory}/{target_name}-{force_field_name}/replica-{replica}`.
If the job running this command is interrupted, it will resume from a binary checkpoint file written by default every 10 ns.

After the production simulations are finished, analyze the trajectories by calling `benchmark_system.analyze_observables(replica={replica})`.
This function will always:
- Produce a new trajectory with solvent atoms stripped and solute atoms aligned to the initial coordinates
- Measure backbone dihedrals, sidechain dihedrals, and tau angles for each residue
- Measure hydrogen bond geometries
If the `target_parameters` dictionary contains an entry with the key `observables`, then specific functions to estimate each observable will also be called.
Currently, the following observables are implemented:
- `scalar_couplings`: NMR three-bond scalar couplings
- `h_bond_scalar_couplings`: NMR interresidue hydrogen bond scalar couplings

The output of the `analyze_observables()` function will be writen to `{result_directory}/{target_name}-{force_field_name}/analysis`.


## Force fields

Force fields are defined by a name and a file path for the force field and a name and optionally a file path for a water model.
Force field files distributed in this repository are located in `proteinbenchmark/data/force-fields`.
Sample force field and water model combinations are included in a dictionary available as the top-level import `force_fields`.
For example, the dictionary for the Amber ff14SB force field with TIP3P water looks like

```
{"ff14sb-tip3p: {"force_field_file": "/path/to/proteinbenchmark/data/force-fields/nerenberg_ff14sb_c0ala_c0gly_c0val.xml",
                 "water_model": "tip3p",
                 "water_model_file": "amber/tip3p_standard.xml"}}
```


## Benchmark targets

Benchmark targets are defined by a dictionary that must be passed to the `target parameters` argument of the `ProteinBenchmarkSystem` constructor.
Sample benchmark targets are included in a dictionary available as the top-level import `benchmark_targets`.
For example, the dictionary for the pentaalanine peptide target system looks like

```
{"ala5": {"target_type": "peptide",
          "aa_sequence": "AAAAA",
          "pressure": Quantity(value=1.0, unit=atomsphere),
          "temperature": Quantity(300.0, unit=kelvin),
          "ph": 2.0,
          "ionic_strength": Quantity(value=0.0, unit=molar),
          "observables": {"scalar_couplings": {"experimental_datasets": "graf_jacs_2007",
                                               "observable_path" : "/path/to/proteinbenchmark/data/observables/ala5/ala5_scalar_couplings.dat"}}}}
```

The `target_parameters` dictionary must contain entries for `pressure`, `temperature`, `ph`, `ionic_strength`, one of `aa_sequence` or `initial_pdb`, and one of `target_type` or `traj_length`.
If `initial_pdb` is present, the initial coordinates will be built from a PDB file located at this file path.
Otherwise, a peptide will be built with the sequence `aa_sequence`.
If building from sequence, backbone dihedrals will be initialized to 180 deg if `build_method` entry is "extended" and to values in an ideal alpha helix if `build_method` is "helical".
Capped termini can be included with entries `nterm_cap` set to "ace" and `cterm_cap` set to "nme" or "nh2".
For both building from a PDB and building from sequence, protomers of titratable groups, including uncapped termini, will be set according to the pH and the initial conformer using `pmx` and `pdb2pqr`.

`target_type` can be one of "peptide", "folded", or "disordered" and is used to set default simulation parameters.
If `target_type` is absent, `traj_length` must be set to an `openmm.unit.Quantity` with units of time.

`observables` contains a dictionary of experimental observables available for this benchmark target.
Datasets of experimental observables distributed in this repository are included in `proteinbenchmark/data/observables`.

Other possible entries to the `target_parameters` dictionary are:
- `solvent_padding`: minimum distance between solute and edge of solvent box
- `nonbonded_cutoff`: nonbonded force cutoff distance for non-SMIRNOFF force fields
- `vdw_switch_width`: distance from `nonbonded_cutoff` at which the switching function is turned on for non-SMIRNOFF force fields
- `restraint_energy_constant`: energy constant for restraints on non-hydrogen solute atoms during energy minimization
- `equil_langevin_friction`: collision frequency for Langevin integrator during equilibration simulation
- `equil_barostat_frequency`: number of steps between attempted volume changes for Monte Carlo barostat during equilibration simulation
- `equil_time_step`: Time step for Langevin integrator during equilibration simulation
- `equil_traj_length`: Total simulation time for equilibration simulation
- `equil_frame_length`: Time between state data and trajectory reports during equilibration simulation
- `langevin_friction`: collision frequency for Langevin integrator during production simulation
- `barostat_frequency`: number of steps between attempted volume changes for Monte Carlo barostat during production simulation
- `time_step`: Time step for Langevin integrator during production simulation
- `traj_length`: Total simulation time for production simulation
- `frame_length`: Time between state data and trajectory reports during production simulation
- `checkpoint_length`: Time between binary checkpoint reports during production simulation
- `save_state_length`: Time between writes to serialized state XML during production simulation

If these values are not present in the `target_parameters` dictionary, then default values will be used from `proteinbenchmark/simulation_parameters`.


## Copyright

Copyright (c) 2022, Open Force Field Initiative, Chapin E. Cavender


## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

