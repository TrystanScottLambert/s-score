"""
Script to run nessie and create a test case for the compare R function.
"""

import polars as pl
import numpy as np
from nessie import FlatCosmology, RedshiftCatalog
from nessie.helper_funcs import create_density_function

INFILE = "/Users/00115372/Desktop/masking_mock_cat/galaxies_shark.parquet"
FRACTIONAL_AREA = 0.00145392285


def create_redshift_catalog(data_frame: pl.DataFrame) -> RedshiftCatalog:
    ras = data_frame["ra"].to_numpy()
    decs = data_frame["dec"].to_numpy()
    redshifts = data_frame["redshift_observed"].to_numpy()
    mock_ids = data_frame["id_fof"].to_numpy()
    cosmo = FlatCosmology(h=0.7, omega_matter=0.3)
    running_density = create_density_function(
        redshifts, len(redshifts), FRACTIONAL_AREA, cosmo
    )
    red_cat = RedshiftCatalog(ras, decs, redshifts, running_density, cosmo)
    red_cat.set_completeness()
    red_cat.mock_group_ids = mock_ids
    return red_cat


def main():
    df_waves = pl.read_parquet(INFILE)

    df_deep = df_waves.filter(
        (pl.col("ra") > 339)
        & (pl.col("ra") < 351)
        & (pl.col("dec") < -30)
        & (pl.col("dec") > -35)
        & (pl.col("mag_Z_VISTA") < 21.25)
    )

    redcat_deep = create_redshift_catalog(df_deep)
    redcat_deep.run_fof(0.05, 30)
    print("nessie s-score:  ", redcat_deep.compare_to_mock(5))
    group_catalog_deep = pl.DataFrame(
        redcat_deep.calculate_group_table(
            velocity_errors=np.ones(len(df_deep["ra"])) * 50,
            absolute_magnitudes=df_deep["mag_abs_r_SDSS"].to_numpy(),
        )
    )
    df_mock = pl.DataFrame(
        {
            "id_galaxy_sky": df_deep["id_galaxy_sky"].to_numpy(),
            "id_fof": df_deep["id_fof"].to_numpy(),
        }
    )

    df_deep = df_deep.with_columns(
        pl.Series(name="id_group", values=redcat_deep.group_ids)
    )
    df_deep = df_deep.rename({"id_galaxy_sky": "id_galaxy"})
    df_deep = df_deep.rename({"id_fof": "id_fof_backup"})
    df_deep.write_parquet("galaxies.parquet")
    group_catalog_deep = group_catalog_deep.rename({"group_id": "id_group"})
    group_catalog_deep.write_parquet("groups.parquet")

    df_mock.write_parquet("mock.parquet")


if __name__ == "__main__":
    main()
