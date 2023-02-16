"""Utilities for file input and output."""
import openmm
from openmm import app
from pathlib import Path


package_data_directory = Path(Path(__file__).parent.absolute(), 'data')


def exists_and_not_empty(file_name):
    """Returns True if file_name exists and is not empty."""

    path = Path(file_name)
    return path.exists() and path.stat().st_size > 0


def read_xml(xml_file_name):
    """Read an OpenMM system from an XML file."""

    with open(xml_file_name, 'r') as xml_file:
        system = openmm.XmlSerializer.deserialize(xml_file.read())

    return system


def remove_model_lines(pdb_file_name, remove_end = False):
    """Remove lines beginning with "MODEL " or "ENDMDL" from a PDB file."""

    # List of strings to exclude
    exclusion_list = ['MODEL ', 'ENDMDL']
    if remove_end:
        exclusion_list.append('END')

    # Copy lines not beginning with elements of exclusion_list to new file
    pdb_path = Path(pdb_file_name)
    new_file_name = f'{pdb_path.stem}-tmp.pdb'
    new_path = Path(new_file_name)
    with open(pdb_path, 'r') as pdb_file:
        with open(new_path, 'w') as new_file:
            for line in pdb_file:
                if line.strip('\n')[:6] not in exclusion_list:
                    new_file.write(line)

    # Replace PDB file with new file
    new_path.replace(pdb_path)


def write_pdb(pdb_file_name, topology, positions):
    """Write a topology and set of positions to a PDB file."""

    with open(pdb_file_name, 'w') as pdb_file:
        app.PDBFile.writeFile(topology, positions, pdb_file)


def write_xml(xml_file_name, openmm_system):
    """Write an OpenMM system to an XML file."""

    with open(xml_file_name, 'w') as xml_file:
        xml_file.write(openmm.XmlSerializer.serialize(openmm_system))


