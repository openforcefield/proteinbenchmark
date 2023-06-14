import click
import numpy
import pandas
from openmm import unit

BOLTZMANN_CONSTANT = 0.001987204259 * unit.kilocalorie_per_mole / unit.kelvin


@click.command()
@click.option(
    "-o",
    "--output_path",
    default="cln025-fraction-folded-by-temperature.dat",
    show_default=True,
    type=click.STRING,
    help="Path to write fraction folded as a function of temperature.",
)
def main(output_path):
    # Fraction helix f as a function of temperature T is given by
    # -k_B T log(1 / f - 1) = \Delta G
    #     = \Delta H - T \Delta S - \Delta C_p (T_0 - T (1 - log(T / T_0)))
    # where \Delta H and \Delta S are functions of the reference temperature T_0
    # and the heat capacity \Delta C_p is independent of temperature

    # Reference values assuming delta_Cp = 0 from Table S1 of
    # Honda S, Akiba T, Kato YS, Sawada Y, Sekijima M, Ishimura M, Ooishi A,
    #     Watanabe H, Odahara T, Harata K. (2008). J. Am. Chem. Soc. 130,
    #     15327-15331.
    ref_Tm = 342.8 * unit.kelvin
    ref_delta_H = 47.7 * unit.kilojoule_per_mole
    ref_delta_S = ref_delta_H / ref_Tm
    ref_delta_Cp = 0.0 * unit.kilocalorie_per_mole / unit.kelvin

    # Array of temperatures between 0 deg C and 100 deg C in Kelvin
    T = numpy.arange(273.15, 373.16, 1.0) * unit.kelvin

    # \Delta G (T) for reference values at 100 deg C
    delta_G = ref_delta_H - T * ref_delta_S

    # Fraction helix as a function of temperature f(T)
    fraction_folded = 1 / (1 + numpy.exp(-delta_G / BOLTZMANN_CONSTANT / T))

    pandas.DataFrame(
        {
            "Temperature (Kelvin)": T,
            "Fraction Folded": fraction_folded,
        }
    ).to_csv(output_path)


if __name__ == "__main__":
    main()
