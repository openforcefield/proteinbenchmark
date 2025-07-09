"""Script to get initial PDBs for benchmark targets."""
from pathlib import Path

import requests

rcsb_url = "https://files.rcsb.org/view/"


def get_rcsb_file(pdb_id, pdb_path, model=None):
    if pdb_path.exists():
        return

    with open(pdb_path, "w") as pdb_file:
        # Get PDB text from RCSB website
        pdb_text = requests.get(f"{rcsb_url}{pdb_id}.pdb").text

        if model is None:
            # Extract ATOM and TER lines for alternate location A
            for line in pdb_text.split("\n"):
                if line[:4] == "ATOM" and line[16] in {" ", "A"}:
                    # Replace deuterium with hydrogen
                    if line[12:16].lstrip()[0] == "D":
                        pdb_file.write(
                            f"{line[:12]}{line[12:16].replace('D', 'H', 1)}"
                            f"{line[16:77]}{line[77].replace('D', 'H')}{line[78:]}\n"
                        )
                    else:
                        pdb_file.write(f"{line}\n")

                elif line[:3] == "TER":
                    pdb_file.write(f"{line}\n")

        else:
            # Extract specified model
            print_lines = False
            for line in pdb_text.split("\n"):
                if line[:6] == "ENDMDL":
                    print_lines = False

                if print_lines:
                    pdb_file.write(f"{line}\n")

                if line[:14] == "MODEL        1":
                    print_lines = True


# Bovine pancreatic trypsin inhibitor (BPTI)
pdb_id = "1PIT"
nmr_model = 1
initial_pdb_path = Path(f"bpti-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

pdb_id = "5PTI"
initial_pdb_path = Path(f"bpti-{pdb_id}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path)

# NESG target CcR55
pdb_id = "2JQN"
nmr_model = 1
initial_pdb_path = Path(f"ccr55-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target CtR148A
pdb_id = "2KO1"
nmr_model = 1
initial_pdb_path = Path(f"ctr148a-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target DhR29B
pdb_id = "2KPU"
nmr_model = 1
initial_pdb_path = Path(f"dhr29b-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target DhR8C
pdb_id = "2KYI"
nmr_model = 1
initial_pdb_path = Path(f"dhr8c-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target ER382A
pdb_id = "2JN0"
nmr_model = 1
initial_pdb_path = Path(f"er382a-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# Streptococcal protein G third immunoglobulin-binding domain (GB3)
pdb_id = "1P7E"
initial_pdb_path = Path(f"gb3-{pdb_id}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path)

pdb_id = "1IGD"
initial_pdb_path = Path(f"gb3-{pdb_id}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path)

# Hen egg white lysozyme (HEWL)
pdb_id = "1E8L"
nmr_model = 1
initial_pdb_path = Path(f"hewl-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

pdb_id = "193L"
initial_pdb_path = Path(f"hewl-{pdb_id}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path)

# NESG target MrR110B
pdb_id = "2K5V"
nmr_model = 1
initial_pdb_path = Path(f"mrr110b-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target PsR293
pdb_id = "2KFP"
nmr_model = 1
initial_pdb_path = Path(f"psr293-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target SR478
pdb_id = "2JS1"
nmr_model = 1
initial_pdb_path = Path(f"sr478-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target SrR115C
pdb_id = "2KCV"
nmr_model = 1
initial_pdb_path = Path(f"srr115c-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# NESG target StR65
pdb_id = "2JN8"
nmr_model = 1
initial_pdb_path = Path(f"str65-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

# Human ubiquitin (Ubq)
pdb_id = "1D3Z"
nmr_model = 1
initial_pdb_path = Path(f"ubq-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)

pdb_id = "1UBQ"
initial_pdb_path = Path(f"ubq-{pdb_id}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path)

# NESG target XcR50
pdb_id = "1XPV"
nmr_model = 1
initial_pdb_path = Path(f"xcr50-{pdb_id}-model-{nmr_model}.pdb")
get_rcsb_file(pdb_id, initial_pdb_path, model=nmr_model)
