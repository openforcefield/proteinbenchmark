import click
import numpy
import pandas
from openmm import unit

BOLTZMANN_CONSTANT = 0.001987204259 * unit.kilocalorie_per_mole / unit.kelvin


@click.command()
@click.option(
    "-o",
    "--output_path",
    default="aaqaa-fraction-helix-by-temperature.dat",
    show_default=True,
    type=click.STRING,
    help="Path to write fraction helix as a function of temperature.",
)
def main(output_path):
    # Fraction helix f as a function of temperature T is given by
    # -k_B T log(1 / f - 1) = \Delta G
    #     = \Delta H - T \Delta S - \Delta C_p (T_0 - T (1 - log(T / T_0)))
    # where \Delta H and \Delta S are functions of the reference temperature T_0
    # and the heat capacity \Delta C_p is independent of temperature

    # Reference values at T_0 = 100 deg C multiplied by number of residues (15)
    # from Table 4 of
    # Shalongo W, Dugad L, Stellwagen E. (1994). J. Am. Chem. Soc. 116,
    #     8288-8293.
    ref_T0 = 373.15 * unit.kelvin
    ref_delta_H = 0.80 * 15 * unit.kilocalorie_per_mole
    ref_delta_S = 0.00286 * 15 * unit.kilocalorie_per_mole / unit.kelvin
    ref_delta_Cp = 0.00245 * 15 * unit.kilocalorie_per_mole / unit.kelvin

    # Array of temperatures between 0 deg C and 100 deg C in Kelvin
    T = numpy.arange(273.15, 373.16, 1.0) * unit.kelvin

    # \Delta G (T) for reference values at 100 deg C
    delta_G = (
        ref_delta_H
        - T * ref_delta_S
        - ref_delta_Cp * (ref_T0 - T * (1 - numpy.log(T / ref_T0)))
    )

    # Fraction helix as a function of temperature f(T)
    fraction_helix = 1 / (1 + numpy.exp(-delta_G / BOLTZMANN_CONSTANT / T))

    pandas.DataFrame(
        {
            "Temperature (Kelvin)": T,
            "Fraction Helix": fraction_helix,
        }
    ).to_csv(output_path)

    # Linear regression for \Delta H and \Delta S for reference temperature of
    # 25 deg C
    # \Delta G + \Delta C_p (T_0 + T log T)
    #     = (\Delta C_p (1 + log T_0) - \Delta S) * T + \Delta H
    T0 = 298.15 * unit.kelvin
    regression_matrix = numpy.vstack([T, numpy.ones(T.size)]).T
    y_values = delta_G + ref_delta_Cp * (T0 + T * numpy.log(T / unit.kelvin))

    fit_result = numpy.linalg.lstsq(regression_matrix, y_values, rcond=None)

    delta_H = fit_result[0][1] * unit.kilocalorie_per_mole
    delta_S = (
        ref_delta_Cp * (1 + numpy.log(T0 / unit.kelvin))
        - fit_result[0][0] * unit.kilocalorie_per_mole / unit.kelvin
    )


if __name__ == "__main__":
    main()
