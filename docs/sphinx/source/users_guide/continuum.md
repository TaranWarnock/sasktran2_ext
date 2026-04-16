---
file_format: mystnb
---

(_users_continuum)=
# Continuum

The [MT_CKD continuum](https://github.com/AER-RC/MT_CKD) is a fortran program which calculates the continuum
absorption due to water vapour, nitrogen, oxygen, carbon dioxide, and ozone over the wavenumber range 0-20,000 cm{sup}`-1`.

## Basic Usage

The MT_CKD continuum can be added to a model atmosphere as a constituent.
See the example in [Quick Start](_quickstart).

The atmosphere must contain constituents that define the H2O, CO2, and O3 VMR profiles. By default, it is assumed that
these species have been added to the atmosphere with the keys "H2O", "CO2", and "O3". Different key names can be used
by setting `h2o_name`, `co2_name`, and `o3_name` in the initialization of {py:class}`sasktran2_ext.continuum.MTCKDContinuum`.
Temperature and pressure profiles must also be set in the atmosphere.

```{note}
The continuum calculated by MT_CKD accounts for contributions that are further than 25 cm{sup}`-1` from each line.
When using the continuum it is recommended to leave the parameter `line_contribution_width` in the optical properties
that do line-by-line calculations, such as {py:class}`sasktran2.optical.hitran.HITRANAbsorber` and {py:class}`sasktran2.optical.hitran.AERLineAbsorber`, at its default value of 25.
```

## Weighting Functions

Weighting functions are calculated with respect to pressure and temperature as well as H2O, CO2, and O3 VMR. In the output
the weighting functions have `continuum` in their variable names to indicate they are portion of the species derivative
resulting from the continuum. For example, if we perform a calculation with the continuum

```{code-cell}
:tags: ["remove-cell"]
import numpy as np
import sasktran2 as sk
from sasktran2_ext.continuum import MTCKDContinuum

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

wavenum = np.arange(530, 550, 0.01)

atmosphere = sk.Atmosphere(geometry, config, wavenumber_cminv=wavenum)
sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
for species in ["H2O", "O3", "CO2"]:
    hitran_db = sk.optical.database.OpticalDatabaseGenericAbsorber(sk.database.StandardDatabase().path(f"hitran/{species}/sasktran2/6e786d8451a300f48627c3239314a9c9214b44bb.nc"))
    atmosphere[species] = sk.climatology.mipas.constituent(species, hitran_db)
atmosphere["emission"] = sk.constituent.ThermalEmission()
atmosphere["surface_emission"] = sk.constituent.SurfaceThermalEmission(300, 0.9)
atmosphere["continuum"] = MTCKDContinuum()

engine = sk.Engine(config, geometry, viewing_geo)
output = engine.calculate_radiance(atmosphere)
```

```{code}
import numpy as np
import sasktran2 as sk
from sasktran2_ext.continuum import MTCKDContinuum

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

wavenum = np.arange(530, 550, 0.01)

atmosphere = sk.Atmosphere(geometry, config, wavenumber_cminv=wavenum)
sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
for species in ["H2O", "O3", "CO2"]:
    hitran_db = sk.database.HITRANDatabase(
        molecule=species,
        start_wavenumber=530,
        end_wavenumber=550,
        wavenumber_resolution=0.01,
        reduction_factor=1,
        backend="sasktran2",
        profile="voigt"
    )
    atmosphere[species] = sk.climatology.mipas.constituent(species, hitran_db)
atmosphere["emission"] = sk.constituent.ThermalEmission()
atmosphere["surface_emission"] = sk.constituent.SurfaceThermalEmission(300, 0.9)
atmosphere["continuum"] = MTCKDContinuum()

engine = sk.Engine(config, geometry, viewing_geo)
output = engine.calculate_radiance(atmosphere)
```

the output will be

```{code-cell}
print(output)
```

Notice there are two weighting functions related to ozone, `wf_O3_vmr` and `wf_continuum_O3_vmr`.
The first contains the ozone weighting function contribution due to line-by-line absorption and
the second is the contribution from continuum absorption.
