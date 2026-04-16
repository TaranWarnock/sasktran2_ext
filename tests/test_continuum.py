from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import sasktran2 as sk

from sasktran2_ext.continuum import MTCKDContinuum


@pytest.fixture()
def test_data_dir():
    """Returns the path to the test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture()
def zero_absorption_file(test_data_dir):
    """Returns the path to the zero absorption file."""
    return test_data_dir / "zero_absorption.nc"


def _test_scenarios(absorption_file):
    config = sk.Config()
    config.emission_source = sk.EmissionSource.Standard
    config.single_scatter_source = sk.SingleScatterSource.NoSource

    altitude_grid = np.arange(0, 65001, 1000.0)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    viewing_geo.add_ray(sk.GroundViewingSolar(0.6, 0, 0.8, 200000))

    wavenum = np.arange(9174.3, 9259.3, 1)

    atmosphere = sk.Atmosphere(geometry, config, wavenumber_cminv=wavenum)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    for species in ["H2O", "O3", "CO2"]:
        opt_prop = sk.optical.database.OpticalDatabaseGenericAbsorber(absorption_file)
        constit = sk.climatology.mipas.constituent(species, opt_prop)
        atmosphere[species] = sk.constituent.VMRAltitudeAbsorber(
            opt_prop,
            altitude_grid,
            np.interp(altitude_grid, constit.altitudes_m, constit.vmr),
            "zero",
        )
    atmosphere["emission"] = sk.constituent.ThermalEmission()
    atmosphere["surface_emission"] = sk.constituent.SurfaceThermalEmission(300, 0.9)
    atmosphere["continuum"] = MTCKDContinuum(
        numeric_wf_fractional_change=1e-2,
        numeric_wf_central_difference=True,
    )

    scen = []

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": atmosphere,
        }
    )

    return scen


def test_continuum_wf(zero_absorption_file):
    scens = _test_scenarios(zero_absorption_file)

    for species in ["H2O", "O3", "CO2"]:
        for scen in scens:
            atmosphere = scen["atmosphere"]

            engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

            radiance = sk.test_util.wf.numeric_wf(
                atmosphere[species].vmr, 0.01, engine, atmosphere, f"wf_{species}_vmr"
            )

            radi = radiance.isel(los=0, stokes=0)
            radi["altitude"] = ("altitude", atmosphere.model_geometry.altitudes())
            radi["CO2_altitude"] = ("CO2_altitude", atmosphere["CO2"].altitudes_m)
            radi["H2O_altitude"] = ("H2O_altitude", atmosphere["H2O"].altitudes_m)
            radi["O3_altitude"] = ("O3_altitude", atmosphere["O3"].altitudes_m)

            sk.test_util.wf.validate_wf(
                radi[f"wf_continuum_{species}_vmr"],
                radi[f"wf_{species}_vmr_numeric"].interp(
                    {f"{species}_altitude": radi.altitude}
                ),
                wf_dim="altitude",
                decimal=5,
            )
