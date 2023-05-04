"""Utilities for file input and output."""
from pathlib import Path

import openmm
import pandas
from openmm import app

package_data_directory = Path(Path(__file__).parent.absolute(), "data")


def exists_and_not_empty(file_name):
    """Returns True if file_name exists and is not empty."""

    path = Path(file_name)
    return path.exists() and path.stat().st_size > 0


def list_of_dicts_to_csv(list_of_dicts, csv_path):
    """Convert a list of dicts to a pandas DataFrame and then write to csv."""

    df = pandas.DataFrame(list_of_dicts)
    df.to_csv(csv_path)


def merge_csvs(csv_prefix):
    """Merge multiple CSV files into one CSV."""

    parent_dir = Path(csv_prefix).parent
    glob_prefix = Path(csv_prefix).name

    # Get list of CSVs
    file_indices = list()

    for csv_file in parent_dir.glob(f"{glob_prefix}-*"):
        file_indices.append(int(csv_file.suffix.split("-")[-1]))

    # Write merged CSV
    df = pandas.DataFrame()

    for i in sorted(file_indices):
        csv_file = Path(f"{csv_prefix}-{i}")
        df = pandas.concat([df, pandas.read_csv(csv_file, index_col=0)])

    df.to_csv(Path(csv_prefix))

    # Delete original CSVs
    for i in sorted(file_indices):
        Path(f"{csv_prefix}-{i}").unlink()


def read_xml(xml_file_name):
    """Read an OpenMM system from an XML file."""

    with open(xml_file_name, "r") as xml_file:
        system = openmm.XmlSerializer.deserialize(xml_file.read())

    return system


def remove_model_lines(pdb_file_name, remove_end=False):
    """Remove lines beginning with "MODEL " or "ENDMDL" from a PDB file."""

    # List of strings to exclude
    exclusion_list = ["MODEL ", "ENDMDL"]
    if remove_end:
        exclusion_list.append("END")

    # Copy lines not beginning with elements of exclusion_list to new file
    pdb_path = Path(pdb_file_name)
    new_file_name = f"{pdb_path.stem}-tmp.pdb"
    new_path = Path(new_file_name)
    with open(pdb_path, "r") as pdb_file:
        with open(new_path, "w") as new_file:
            for line in pdb_file:
                if line.strip("\n")[:6] not in exclusion_list:
                    new_file.write(line)

    # Replace PDB file with new file
    new_path.replace(pdb_path)


def write_pdb(pdb_file_name, topology, positions):
    """Write a topology and set of positions to a PDB file."""

    with open(pdb_file_name, "w") as pdb_file:
        app.PDBFile.writeFile(topology, positions, pdb_file)


def write_xml(xml_file_name, openmm_system):
    """Write an OpenMM system to an XML file."""

    with open(xml_file_name, "w") as xml_file:
        xml_file.write(openmm.XmlSerializer.serialize(openmm_system))
