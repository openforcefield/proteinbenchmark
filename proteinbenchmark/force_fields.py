"""List of force fields and water models."""

from pathlib import Path

from proteinbenchmark.utilities import package_data_directory

ff_directory = Path(package_data_directory, "force-fields")

# List of force fields with force field XML file, water model, and water model
# XML file
force_fields = {
    "espaloma-0.3.2-opc": {
        "force_field_file": Path(ff_directory, "espaloma-0.3.2.pt"),
        "water_model": "opc",
    },
    "espaloma-0.3.2-opc3": {
        "force_field_file": Path(ff_directory, "espaloma-0.3.2.pt"),
        "water_model": "opc3",
    },
    "espaloma-0.3.2-tip3p": {
        "force_field_file": Path(ff_directory, "espaloma-0.3.2.pt"),
        "water_model": "tip3p",
    },
    "espaloma-0.3.2-tip3p-fb": {
        "force_field_file": Path(ff_directory, "espaloma-0.3.2.pt"),
        "water_model": "tip3p-fb",
    },
    "espaloma-0.3.2-tip4p-fb": {
        "force_field_file": Path(ff_directory, "espaloma-0.3.2.pt"),
        "water_model": "tip4p-fb",
    },
    "ff14sb-opc": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sb_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "opc",
    },
    "ff14sb-opc3": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sb_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "opc3",
    },
    "ff14sb-tian-opc": {
        "force_field_file": Path(ff_directory, "tian_ff14sb_c0ala.xml"),
        "water_model": "opc",
    },
    "ff14sb-tian-opc3": {
        "force_field_file": Path(ff_directory, "tian_ff14sb_c0ala.xml"),
        "water_model": "opc3",
    },
    "ff14sb-tian-tip3p": {
        "force_field_file": Path(ff_directory, "tian_ff14sb_c0ala.xml"),
        "water_model": "tip3p",
    },
    "ff14sb-tian-tip3p-fb": {
        "force_field_file": Path(ff_directory, "tian_ff14sb_c0ala.xml"),
        "water_model": "tip3p-fb",
    },
    "ff14sb-tian-tip4p-fb": {
        "force_field_file": Path(ff_directory, "tian_ff14sb_c0ala.xml"),
        "water_model": "tip4p-fb",
    },
    "ff14sb-tip3p": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sb_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "tip3p",
    },
    "ff14sb-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sb_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "tip3p-fb",
    },
    "ff14sb-tip4p-fb": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sb_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "tip4p-fb",
    },
    "ff14sbonlysc-opc": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sbonlysc_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "opc",
    },
    "ff14sbonlysc-opc3": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sbonlysc_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "opc3",
    },
    "ff14sbonlysc-tip3p": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sbonlysc_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "tip3p",
    },
    "ff14sbonlysc-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sbonlysc_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "tip3p-fb",
    },
    "ff14sbonlysc-tip4p-fb": {
        "force_field_file": Path(
            ff_directory, "nerenberg_ff14sbonlysc_c0ala_c0gly_c0val.xml"
        ),
        "water_model": "tip4p-fb",
    },
    "null-0.0.1-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.1.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.2-nbamber-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NBAmber.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.2-nbamber-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NBAmber.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.2-nmr-ala5-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-Ala5-OPC3.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.2-nmr-ala5-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.2-Ala5-TIP3P-FB.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.2-nmr-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-OPC3.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.2-nmr-2-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-OPC3-2.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.2-nmr-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-TIP3P-FB.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.2-nmr-2-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-TIP3P-FB-2.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.2-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NH2.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.2-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NH2.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.2-qamber-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-QAmber.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.2-qamber-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-QAmber.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.2-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NH2.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.2-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NH2.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.2-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.2-NH2.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-ai-dw-opc": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-abinitio-default-weights.offxml"
        ),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-ai-dw-opc3": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-abinitio-default-weights.offxml"
        ),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-ai-dw-tip3p": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-abinitio-default-weights.offxml"
        ),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-ai-dw-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-abinitio-default-weights.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-ai-dw-tip4p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-abinitio-default-weights.offxml"
        ),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-ai-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-abinitio.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-ai-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-abinitio.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-ai-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-abinitio.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-ai-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-abinitio.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-ai-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-abinitio.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-dw-opc": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-default-weights.offxml"
        ),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-dw-opc3": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-default-weights.offxml"
        ),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-dw-tip3p": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-default-weights.offxml"
        ),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-dw-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-default-weights.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-dw-tip4p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Null-0.0.3-default-weights.offxml"
        ),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-nagl-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NAGL.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-nagl-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-nagl-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NAGL.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-nagl-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NAGL.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-nagl-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NAGL.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-nbamber-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NBAmber.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-nbamber-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NBAmber.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-nbamber-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NBAmber.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-nbamber-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NBAmber.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-nbamber-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-NBAmber.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-qamber-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-QAmber.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-qamber-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-QAmber.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-qamber-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-QAmber.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-qamber-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-QAmber.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-qamber-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-QAmber.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-pair-general-nmr-1e2-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-general-nmr-1e2-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-general-nmr-1e3-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-general-nmr-1e3-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-general-nmr-1e4-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-general-nmr-1e4-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-general-nmr-1e5-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-general-nmr-1e5-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e2-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e2-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e3-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e3-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-2-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-2-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-3-0.7-1e3-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-3-0.7-1e3-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-3-0.7-1e4-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-3-0.7-1e4-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-3-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-3-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-4-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-4-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-5-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-5-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e4-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e4-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-0.8-1e5-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-0.8-1e5-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e2-2-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e2-2-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e2-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e2-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e3-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e3-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e4-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e4-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e5-2-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e5-2-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e5-3-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e5-3-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-nmr-1e5-opc3": {
        "force_field_file": Path(ff_directory, "null-0.0.3-pair-nmr-1e5-opc3-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-opc": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "null-0.0.3-pair-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "null-0.0.3-pair-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-pair-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-pair-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "null-0.0.3-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "null-0.0.3-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "null-0.0.3-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Null-0.0.3.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "sage-2.1.0-nagl-opc": {
        "force_field_file": Path(ff_directory, "Sage-2.1.0-NAGL.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "sage-2.1.0-nagl-opc3": {
        "force_field_file": Path(ff_directory, "Sage-2.1.0-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "sage-2.1.0-nagl-tip3p": {
        "force_field_file": Path(ff_directory, "Sage-2.1.0-NAGL.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "sage-2.1.0-nagl-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Sage-2.1.0-NAGL.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "sage-2.1.0-nagl-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Sage-2.1.0-NAGL.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "specific-0.0.1-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.1.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.2-nbamber-opc": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NBAmber.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.2-nbamber-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NBAmber.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.2-nmr-ala5-opc3": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.2-Ala5-OPC3.offxml"
        ),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.2-nmr-ala5-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.2-Ala5-TIP3P-FB.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.2-nmr-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-OPC3.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.2-nmr-2-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-OPC3-2.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.2-nmr-2-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.2-TIP3P-FB-2.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.2-nmr-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.2-TIP3P-FB.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.2-opc": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NH2.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.2-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NH2.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.2-qamber-opc": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-QAmber.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.2-qamber-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-QAmber.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.2-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NH2.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.2-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NH2.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.2-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.2-NH2.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "specific-0.0.3-ai-opc": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-abinitio.offxml"
        ),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.3-ai-opc3": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-abinitio.offxml"
        ),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.3-ai-tip3p": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-abinitio.offxml"
        ),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.3-ai-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-abinitio.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.3-ai-tip4p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-abinitio.offxml"
        ),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "specific-0.0.3-dw-opc": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-default-weights.offxml"
        ),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.3-dw-opc3": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-default-weights.offxml"
        ),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.3-dw-tip3p": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-default-weights.offxml"
        ),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.3-dw-tip3p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-default-weights.offxml"
        ),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.3-dw-tip4p-fb": {
        "force_field_file": Path(
            ff_directory, "Protein-Specific-0.0.3-default-weights.offxml"
        ),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "specific-0.0.3-opc": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.3-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.3-pair-opc": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.3-pair-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.3-pair-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.3-pair-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.3-pair-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Pairwise-NAGL.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "specific-0.0.3-sage-pair-opc": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Sage-Pairwise-NAGL.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "specific-0.0.3-sage-pair-opc3": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Sage-Pairwise-NAGL.offxml"),
        "water_model": "opc3",
        "water_model_file": Path(ff_directory, "opc3-1.0.0.offxml"),
    },
    "specific-0.0.3-sage-pair-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Sage-Pairwise-NAGL.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.3-sage-pair-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Sage-Pairwise-NAGL.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.3-sage-pair-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3-Sage-Pairwise-NAGL.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "specific-0.0.3-tip3p": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3.offxml"),
        "water_model": "tip3p",
        "water_model_file": None,
    },
    "specific-0.0.3-tip3p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3.offxml"),
        "water_model": "tip3p-fb",
        "water_model_file": Path(ff_directory, "tip3p_fb-1.1.0.offxml"),
    },
    "specific-0.0.3-tip4p-fb": {
        "force_field_file": Path(ff_directory, "Protein-Specific-0.0.3.offxml"),
        "water_model": "tip4p-fb",
        "water_model_file": Path(ff_directory, "tip4p_fb-1.0.0.offxml"),
    },
    "amber-ff-default-weight-tip3p": {
        "force_field_file": Path(ff_directory, "amber-ff-default-weight.offxml"),
        "water_model": "tip3p",
        "water_model_file": "tip3p-1.0.1.offxml",
    },
    "amber-ff-default-weight-opc": {
        "force_field_file": Path(ff_directory, "amber-ff-default-weight.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
    "amber-ff-equal-weights-tip3p": {
        "force_field_file": Path(ff_directory, "amber-ff-equal-weights.offxml"),
        "water_model": "tip3p",
        "water_model_file": "tip3p-1.0.1.offxml",
    },
    "amber-ff-equal-weights-opc": {
        "force_field_file": Path(ff_directory, "amber-ff-equal-weights.offxml"),
        "water_model": "opc",
        "water_model_file": "opc-1.0.0.offxml",
    },
}

# Add implicit water model files
water_model_files = {
    "tip3p": "amber/tip3p_standard.xml",
    "opc3": Path(ff_directory, "openmm-opc3.xml"),
    "opc": "amber/opc_standard.xml",
    "tip3p-fb": Path(ff_directory, "openmm-tip3p-fb.xml"),
    "tip4p-fb": Path(ff_directory, "openmm-tip4p-fb.xml"),
}

for force_field_name, ff_parameters in force_fields.items():
    if "water_model_file" not in ff_parameters:
        water_model_file = water_model_files[ff_parameters["water_model"]]
        ff_parameters["water_model_file"] = water_model_file
