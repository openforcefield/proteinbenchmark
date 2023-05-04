"""Script to get initial PDBs for benchmark targets."""
from pathlib import Path

import requests

rcsb_url = "https://files.rcsb.org/view/"

# Hen egg white lysozyme (HEWL)
pdb_id = "1P7E"
initial_pdb = "gb3-1P7E.pdb"

if not Path(f"{pdb_id}.pdb").exists():
    with open(initial_pdb, "w") as pdb_file:
        # Get PDB text from RCSB website
        pdb_text = requests.get(f"{rcsb_url}{pdb_id}.pdb").text

        # Extract ATOM and TER lines
        for line in pdb_text.split("\n"):
            if line[:4] in ["ATOM", "TER "]:
                pdb_file.write(f"{line}\n")


# Hen egg white lysozyme (HEWL)
pdb_id = "1E8L"
initial_pdb = "hewl-1E8L-model-1.pdb"

if not Path(f"{pdb_id}.pdb").exists():
    with open(initial_pdb, "w") as pdb_file:
        # Get PDB text from RCSB website
        pdb_text = requests.get(f"{rcsb_url}{pdb_id}.pdb").text

        # Extract first model
        print_lines = False
        for line in pdb_text.split("\n"):
            if line[:6] == "ENDMDL":
                print_lines = False

            if print_lines:
                pdb_file.write(f"{line}\n")

            if line[:14] == "MODEL        1":
                print_lines = True
