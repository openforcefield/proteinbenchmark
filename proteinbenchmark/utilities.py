"""Utilities for file input and output."""

from operator import itemgetter
from pathlib import Path

import loos
import numpy
import openmm
import pandas
from openmm import app

package_data_directory = Path(Path(__file__).parent.absolute(), "data")


def copy_internal_coords(ra, rb, rc, rd, ta, tb, tc):
    """
    Given four reference atoms and three target atoms, return the coordinates of
    a fourth target atom so that it has the same distance, angle, and dihedral
    as the reference atoms.
    """

    distance = rd.coords().distance(rc.coords())
    angle = numpy.deg2rad(loos.angle(rd, rc, rb))
    dihedral = numpy.deg2rad(loos.torsion(rd, rc, rb, ra))
    r_cb = tb.coords() - tc.coords()
    r_cb /= r_cb.length()
    r_ba = ta.coords() - tb.coords()
    r_ba /= r_ba.length()
    ba_cross_cb = r_ba.cross(r_cb)
    ba_cross_cb /= ba_cross_cb.length()
    z = r_cb
    y = z.cross(ba_cross_cb)
    y /= y.length()
    x = y.cross(z)
    x /= x.length()
    dx = distance * numpy.sin(angle) * numpy.sin(dihedral)
    dy = distance * numpy.sin(angle) * numpy.cos(dihedral)
    dz = distance * numpy.cos(angle)
    r_cd = numpy.dot(
        numpy.array([[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]]),
        numpy.array([dx, dy, dz]),
    )
    return numpy.array([tc.coords()[0], tc.coords()[1], tc.coords()[2]]) + r_cd


def exists_and_not_empty(file_name):
    """Returns True if file_name exists and is not empty."""

    path = Path(file_name)
    return path.exists() and path.stat().st_size > 0


def extract_noe_upper_distances(input_path: str, output_path: str):
    """Extract NOE upper distance boundaries from an NMR STAR file."""

    read_noe_distances = False
    constraint_type = None

    noe_upper_distances = list()

    with open(input_path, "r") as input_file:
        for line_index, line in enumerate(input_file):
            fields = line.split()

            if len(fields) == 0:
                continue

            if fields[0] == "_Gen_dist_constraint_list.Constraint_type":
                constraint_type = fields[1]

            elif read_noe_distances and fields[0] == "stop_":
                read_noe_distances = False

            # Only read the first restraint for each entry
            elif read_noe_distances and fields[1] == "1":
                # TODO: are the field indices guaranteed to be the same for all
                # STAR files?
                resid_i = int(fields[36])
                resname_i = fields[37]
                atom_i = fields[51]
                resid_j = int(fields[43])
                resname_j = fields[44]
                atom_j = fields[59]
                distance = float(fields[28])

                if resid_j < resid_i:
                    resid_i, resid_j = resid_j, resid_i
                    resname_i, resname_j = resname_j, resname_i
                    atom_i, atom_j = atom_j, atom_i

                noe_upper_distances.append(
                    (resid_i, resname_i, atom_i, resid_j, resname_j, atom_j, distance)
                )

            elif (
                fields[0] == "_Gen_dist_constraint.Gen_dist_constraint_list_ID"
                and constraint_type == "NOE"
            ):
                read_noe_distances = True

    # Sort by (resid_i, atom_i, resid_j, atom_j) by sorting one-at-a-time in
    # reverse order. For atom names, sort by Greek character in the second
    # position of the atom string
    atom_sort_order = ["N", "A", "B", "G", "D", "E", "Z", "H"]
    atom_sort_dict = {atom: index for index, atom in enumerate(atom_sort_order)}

    noe_upper_distances.sort(key=lambda t: atom_sort_dict[t[5][1]])
    noe_upper_distances.sort(key=itemgetter(3))
    noe_upper_distances.sort(key=lambda t: atom_sort_dict[t[2][1]])
    noe_upper_distances.sort(key=itemgetter(0))

    pseudoatom_name_map = dict()

    for resname in [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]:
        pseudoatom_name_map[resname] = dict()

        if resname != "PRO":
            pseudoatom_name_map[resname]["HN"] = "H"

        if resname == "GLY":
            pseudoatom_name_map[resname]["HA#"] = "PA"
        else:
            pseudoatom_name_map[resname]["HA"] = "HA"

        if resname == "ALA":
            pseudoatom_name_map[resname]["HB#"] = "MB"
        elif resname == "ILE":
            pseudoatom_name_map[resname]["HB"] = "HB"
            pseudoatom_name_map[resname]["HB#"] = "HB"
        elif resname in {"THR", "VAL"}:
            pseudoatom_name_map[resname]["HB"] = "HB"
        else:
            pseudoatom_name_map[resname]["HB2"] = "HB2"
            pseudoatom_name_map[resname]["HB1"] = "HB3"
            pseudoatom_name_map[resname]["HB#"] = "PB"

    pseudoatom_name_map["ARG"]["HG#"] = "PG"
    pseudoatom_name_map["ARG"]["HD#"] = "PD"
    pseudoatom_name_map["ASN"]["HD#"] = "ND2"
    pseudoatom_name_map["GLN"]["HG#"] = "PG"
    pseudoatom_name_map["GLU"]["HG#"] = "PG"
    pseudoatom_name_map["HIS"]["HD2"] = "HD2"
    pseudoatom_name_map["HIS"]["HE1"] = "HE1"
    pseudoatom_name_map["ILE"]["HG1#"] = "PG1"
    pseudoatom_name_map["ILE"]["HG2#"] = "MG2"
    pseudoatom_name_map["ILE"]["HG*"] = "IG"
    pseudoatom_name_map["ILE"]["HD#"] = "MD1"
    pseudoatom_name_map["LEU"]["HG"] = "HG"
    pseudoatom_name_map["LEU"]["HD*"] = "QD"
    pseudoatom_name_map["LYS"]["HG#"] = "PG"
    pseudoatom_name_map["LYS"]["HD#"] = "PD"
    pseudoatom_name_map["LYS"]["HE#"] = "PE"
    pseudoatom_name_map["MET"]["HG#"] = "PG"
    pseudoatom_name_map["MET"]["HE#"] = "ME"
    pseudoatom_name_map["PHE"]["HD*"] = "RD"
    pseudoatom_name_map["PHE"]["HE*"] = "RE"
    pseudoatom_name_map["PHE"]["HZ"] = "HZ"
    pseudoatom_name_map["PRO"]["HG#"] = "PG"
    pseudoatom_name_map["PRO"]["HD#"] = "PD"
    pseudoatom_name_map["THR"]["HG2#"] = "MG2"
    pseudoatom_name_map["TRP"]["HD1"] = "HD1"
    pseudoatom_name_map["TRP"]["HE1"] = "HE1"
    pseudoatom_name_map["TRP"]["HE3"] = "HE3"
    pseudoatom_name_map["TRP"]["HZ2"] = "HZ2"
    pseudoatom_name_map["TRP"]["HZ3"] = "HZ3"
    pseudoatom_name_map["TRP"]["HH2"] = "HH2"
    pseudoatom_name_map["TYR"]["HD*"] = "RD"
    pseudoatom_name_map["TYR"]["HE*"] = "RE"
    pseudoatom_name_map["VAL"]["HG1#"] = "MG1"
    pseudoatom_name_map["VAL"]["HG2#"] = "MG2"
    pseudoatom_name_map["VAL"]["HG*"] = "QG"

    pseudoatom_corrections = {
        "H": 0.0,
        "P": 1.0,
        "M": 1.5,
        "N": 1.0,
        "I": 1.0,
        "R": 2.4,
        "Q": 2.4,
    }

    # Write NOE upper distances to output path
    with open(output_path, "w") as output_file:
        output_file.write(
            "Observable         Resid_i Resname_i Atom_i Resid_j Resname_j "
            "Atom_j Experiment"
        )

        for (
            resid_i,
            resname_i,
            atom_i,
            resid_j,
            resname_j,
            atom_j,
            distance,
        ) in noe_upper_distances:
            try:
                pseudoatom_i = pseudoatom_name_map[resname_i][atom_i]
            except KeyError:
                raise ValueError(
                    f"Atom {atom_i} of residue {resname_i} {resid_i} is not "
                    "in the pseudoatom map."
                )

            try:
                pseudoatom_j = pseudoatom_name_map[resname_j][atom_j]
            except KeyError:
                raise ValueError(
                    f"Atom {atom_j} of residue {resname_j} {resid_j} is not "
                    "in the pseudoatom map."
                )

            raw_distance = distance
            for pseudoatom in (pseudoatom_i, pseudoatom_j):
                if pseudoatom != "ME":
                    raw_distance -= pseudoatom_corrections[pseudoatom[0]]

            output_file.write(
                f"\nNOE_upper_distance {resid_i:3d}     {resname_i:3s}       "
                f"{pseudoatom_i:4s}   {resid_j:3d}     {resname_j:3s}       "
                f"{pseudoatom_j:4s}   {raw_distance:4.1f}"
            )


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
