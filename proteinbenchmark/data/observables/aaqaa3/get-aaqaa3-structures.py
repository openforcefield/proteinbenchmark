import click
from proteinbenchmark.benchmark_targets import benchmark_targets
from proteinbenchmark.system_setup import build_initial_coordinates


@click.command()
@click.option(
    "-o",
    "--output_prefix",
    default="aaqaa3",
    show_default=True,
    type=click.STRING,
    help="Prefix for path to write PDB files."
)
def main(output_prefix):

    aaqaa3_parameters = benchmark_targets["aaqaa3"]
    for build_method in ["alpha", "beta", "PII", "extended"]:
        build_initial_coordinates(
            build_method=build_method,
            ph=aaqaa3_parameters["ph"],
            initial_pdb=f"{output_prefix}-{build_method}-initial.pdb",
            protonated_pdb=f"{output_prefix}-{build_method}-protonated.pdb",
            aa_sequence=aaqaa3_parameters["aa_sequence"],
            nterm_cap=aaqaa3_parameters["nterm_cap"],
            cterm_cap=aaqaa3_parameters["cterm_cap"],
        )


if __name__ == "__main__":
    main()
