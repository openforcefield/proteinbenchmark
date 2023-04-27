"""List of force fields and water models."""
from pathlib import Path
from proteinbenchmark.utilities import package_data_directory


# List of force fields with force field XML file, water model, and water model
# XML file
force_fields = {
    'ff14sb-opc': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'nerenberg_ff14sb_c0ala_c0gly_c0val.xml',
        ),
        'water_model': 'opc',
    },
    'ff14sb-tian-opc': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'tian_ff14sb_c0ala.xml',
        ),
        'water_model': 'opc',
    },
    'ff14sb-tian-tip3p': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'tian_ff14sb_c0ala.xml',
        ),
        'water_model': 'tip3p',
    },
    'ff14sb-tip3p': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'nerenberg_ff14sb_c0ala_c0gly_c0val.xml',
        ),
        'water_model': 'tip3p',
    },
    'null-0.0.1-tip3p': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'Protein-Null-0.0.1.offxml',
        ),
        'water_model': 'tip3p',
        'water_model_file': None,
    },
    'null-0.0.2-opc': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'Protein-Null-0.0.2-NH2.offxml',
        ),
        'water_model': 'opc',
        'water_model_file': Path(
            package_data_directory,
            'force-fields',
            'opc-1.0.0.offxml',
        ),
    },
    'null-0.0.2-tip3p': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'Protein-Null-0.0.2-NH2.offxml',
        ),
        'water_model': 'tip3p',
        'water_model_file': None,
    },
    'specific-0.0.1-tip3p': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'Protein-Specific-0.0.1.offxml',
        ),
        'water_model': 'tip3p',
        'water_model_file': None,
    },
    'specific-0.0.2-opc': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'Protein-Specific-0.0.2-NH2.offxml',
        ),
        'water_model': 'opc',
        'water_model_file': Path(
            package_data_directory,
            'force-fields',
            'opc-1.0.0.offxml',
        ),
    },
    'specific-0.0.2-tip3p': {
        'force_field_file': Path(
            package_data_directory,
            'force-fields',
            'Protein-Specific-0.0.2-NH2.offxml',
        ),
        'water_model': 'tip3p',
        'water_model_file': None,
    },
}

# Add implicit water model files
water_model_files = {
    'tip3p': 'amber/tip3p_standard.xml',
    'opc': 'amber/opc_standard.xml'
}

for force_field_name, ff_parameters in force_fields.items():

    if 'water_model_file' not in ff_parameters:

        water_model_file = water_model_files[ff_parameters['water_model']]
        ff_parameters['water_model_file'] = water_model_file

