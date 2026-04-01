"""Utilities for file input and output."""

from operator import itemgetter
from pathlib import Path

import loos
import numpy
import openmm
import pandas
from openmm import app
from openff.toolkit import Molecule

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


def enforce_iupac_branch_numbering(offmol: Molecule, conformer_index: int = 0):
    """
    Ensure that a molecule with a HierarchyScheme for canonical amino acids
    follows the IUPAC rules for branch numbering in tetrahedral configurations
    (https://iupac.qmul.ac.uk/misc/ppep1.html#222) for a particular conformer.
    """

    if "residues" not in offmol.hierarchy_schemes:
        return

    conformer = offmol.conformers[conformer_index]
    new_conformer_indices = numpy.arange(conformer.shape[0])

    # Construct a dictionary that maps residue_number and atom name to molecule
    # atom index
    atom_dict = {
        int(residue.residue_number): {
            atom.metadata["atom_name"]: atom.molecule_atom_index
            for atom in residue.atoms
        }
        for residue in offmol.residues
    }
        
    for residue in offmol.residues:
        resnum = int(residue.residue_number)

        # Alpha methylene HA2 and HA3
        if residue.residue_name == "GLY":
            n = atom_dict[resnum]["N"]
            ca = atom_dict[resnum]["CA"]
            c = atom_dict[resnum]["C"]
            ha2 = atom_dict[resnum]["HA2"]
            ha3 = atom_dict[resnum]["HA3"]

            # Select GLY HA2 and HA3 based on phi, but use psi for first residue
            if resnum == 1:
                n_next = atom_dict[resnum + 1]["N"]

                psi = measure_dihedral(conformer, n, ca, c, n_next)
                ha2_torsion = measure_dihedral(conformer, ha2, ca, c, n_next)
                ha3_torsion = measure_dihedral(conformer, ha3, ca, c, n_next)

                if (ha3_torsion - psi) % 360 > (ha2_torsion - psi) % 360:
                    new_conformer_indices[ha2] = ha3
                    new_conformer_indices[ha3] = ha2

            else:
                c_prev = atom_dict[resnum - 1]["C"]

                phi = measure_dihedral(conformer, c_prev, n, ca, c)
                ha2_torsion = measure_dihedral(conformer, c_prev, n, ca, ha2)
                ha3_torsion = measure_dihedral(conformer, c_prev, n, ca, ha3)

                if (ha2_torsion - phi) % 360 > (ha3_torsion - phi) % 360:
                    new_conformer_indices[ha2] = ha3
                    new_conformer_indices[ha3] = ha2

        # Beta methylene HB2 and HB3
        if residue.residue_name in {
            "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "TYR", "TRP",
        }:
            n = atom_dict[resnum]["N"]
            ca = atom_dict[resnum]["CA"]
            cb = atom_dict[resnum]["CB"]
            hb2 = atom_dict[resnum]["HB2"]
            hb3 = atom_dict[resnum]["HB3"]

            if residue.residue_name == "CYS":
                ag = atom_dict[resnum]["SG"]
            elif residue.residue_name == "SER":
                ag = atom_dict[resnum]["OG"]
            else:
                ag = atom_dict[resnum]["CG"]

            chi1 = measure_dihedral(conformer, n, ca, cb, ag)
            hb2_torsion = measure_dihedral(conformer, n, ca, cb, hb2)
            hb3_torsion = measure_dihedral(conformer, n, ca, cb, hb3)

            if (hb2_torsion - chi1) % 360 > (hb3_torsion - chi1) % 360:
                new_conformer_indices[hb2] = hb3
                new_conformer_indices[hb3] = hb2

        # Gamma methylene HG2 and HG3
        if residue.residue_name in {"ARG", "GLN", "GLU", "LYS", "MET", "PRO"}:
            ca = atom_dict[resnum]["CA"]
            cb = atom_dict[resnum]["CB"]
            cg = atom_dict[resnum]["CG"]
            hg2 = atom_dict[resnum]["HG2"]
            hg3 = atom_dict[resnum]["HG3"]

            if residue.residue_name == "MET":
                ad = atom_dict[resnum]["SD"]
            else:
                ad = atom_dict[resnum]["CD"]

            chi2 = measure_dihedral(conformer, ca, cb, cg, ad)
            hg2_torsion = measure_dihedral(conformer, ca, cb, cg, hg2)
            hg3_torsion = measure_dihedral(conformer, ca, cb, cg, hg3)

            if (hg2_torsion - chi2) % 360 > (hg3_torsion - chi2) % 360:
                new_conformer_indices[hg2] = hg3
                new_conformer_indices[hg3] = hg2

        # Delta methylene HD2 and HD3
        if residue.residue_name in {"ARG", "LYS", "PRO"}:
            cb = atom_dict[resnum]["CB"]
            cg = atom_dict[resnum]["CG"]
            cd = atom_dict[resnum]["CD"]
            hd2 = atom_dict[resnum]["HD2"]
            hd3 = atom_dict[resnum]["HD3"]

            if residue.residue_name == "ARG":
                ae = atom_dict[resnum]["NE"]
            elif residue.residue_name == "PRO":
                ae = atom_dict[resnum]["N"]
            else:
                ae = atom_dict[resnum]["CE"]

            chi3 = measure_dihedral(conformer, cb, cg, cd, ae)
            hd2_torsion = measure_dihedral(conformer, cb, cg, cd, hd2)
            hd3_torsion = measure_dihedral(conformer, cb, cg, cd, hd3)

            if (hd2_torsion - chi3) % 360 > (hd3_torsion - chi3) % 360:
                new_conformer_indices[hd2] = hd3
                new_conformer_indices[hd3] = hd2

        # Epsilon methylene HE2 and HE3
        if residue.residue_name == "LYS":
            cg = atom_dict[resnum]["CG"]
            cd = atom_dict[resnum]["CD"]
            ce = atom_dict[resnum]["CE"]
            nz = atom_dict[resnum]["NZ"]
            he2 = atom_dict[resnum]["HE2"]
            he3 = atom_dict[resnum]["HE3"]

            chi4 = measure_dihedral(conformer, cg, cd, ce, nz)
            he2_torsion = measure_dihedral(conformer, cg, cd, ce, he2)
            he3_torsion = measure_dihedral(conformer, cg, cd, ce, he3)

            if (he2_torsion - chi4) % 360 > (he3_torsion - chi4) % 360:
                new_conformer_indices[he2] = he3
                new_conformer_indices[he3] = he2

        # Beta isopropyl CG1 and CG2
        if residue.residue_name == "VAL":
            n = atom_dict[resnum]["N"]
            ca = atom_dict[resnum]["CA"]
            cb = atom_dict[resnum]["CB"]
            hb = atom_dict[resnum]["HB"]
            cg1 = atom_dict[resnum]["CG1"]
            cg2 = atom_dict[resnum]["CG2"]

            chi1 = measure_dihedral(conformer, n, ca, cb, cg1)
            cg2_torsion = measure_dihedral(conformer, n, ca, cb, cg2)
            hb_torsion = measure_dihedral(conformer, n, ca, cb, hb)

            if (chi1 - hb_torsion) % 360 > (cg2_torsion - hb_torsion) % 360:
                new_conformer_indices[cg1] = cg2
                new_conformer_indices[atom_dict[resnum]["HG11"]] = atom_dict[resnum]["HG21"]
                new_conformer_indices[atom_dict[resnum]["HG12"]] = atom_dict[resnum]["HG22"]
                new_conformer_indices[atom_dict[resnum]["HG13"]] = atom_dict[resnum]["HG23"]
                new_conformer_indices[cg2] = cg1
                new_conformer_indices[atom_dict[resnum]["HG21"]] = atom_dict[resnum]["HG11"]
                new_conformer_indices[atom_dict[resnum]["HG22"]] = atom_dict[resnum]["HG12"]
                new_conformer_indices[atom_dict[resnum]["HG23"]] = atom_dict[resnum]["HG13"]

        # Gamma isopropyl CD1 and CD2
        if residue.residue_name == "LEU":
            ca = atom_dict[resnum]["CA"]
            cb = atom_dict[resnum]["CB"]
            cg = atom_dict[resnum]["CG"]
            hg = atom_dict[resnum]["HG"]
            cd1 = atom_dict[resnum]["CD1"]
            cd2 = atom_dict[resnum]["CD2"]

            chi2 = measure_dihedral(conformer, ca, cb, cg, cd1)
            cd2_torsion = measure_dihedral(conformer, ca, cb, cg, cd22)
            hg_torsion = measure_dihedral(conformer, ca, cb, cg, hg)

            if (chi2 - hg_torsion) % 360 > (cd2_torsion - hg_torsion) % 360:
                new_conformer_indices[cd1] = cd2
                new_conformer_indices[atom_dict[resnum]["HD11"]] = atom_dict[resnum]["HD21"]
                new_conformer_indices[atom_dict[resnum]["HD12"]] = atom_dict[resnum]["HD22"]
                new_conformer_indices[atom_dict[resnum]["HD13"]] = atom_dict[resnum]["HD23"]
                new_conformer_indices[cd2] = cd1
                new_conformer_indices[atom_dict[resnum]["HD21"]] = atom_dict[resnum]["HD11"]
                new_conformer_indices[atom_dict[resnum]["HD22"]] = atom_dict[resnum]["HD12"]
                new_conformer_indices[atom_dict[resnum]["HD23"]] = atom_dict[resnum]["HD13"]

    offmol.conformers[conformer_index] = conformer[new_conformer_indices]


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


def measure_dihedral(conformer: numpy.typing.ArrayLike, i: int, j: int, k: int, l: int) -> float:
    """Measure a dihedral for atoms with indices (i, j, k, l) in degrees."""

    r_ij = conformer[j] - conformer[i]
    r_jk = conformer[k] - conformer[j]
    r_kl = conformer[l] - conformer[k]

    cross_jkl = numpy.cross(r_jk, r_kl)
    dihedral = numpy.arctan2(
        numpy.dot(numpy.linalg.norm(r_jk) * r_ij, cross_jkl),
        numpy.dot(numpy.cross(r_ij, r_jk), cross_jkl),
    )

    return numpy.rad2deg(dihedral).m

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
