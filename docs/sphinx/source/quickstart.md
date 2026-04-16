---
file_format: mystnb
---

(_quickstart)=
# Quick Start

SASKTRAN2 extensions are designed to be used with the base {py:mod}`sasktran2`.
The MT_CKD continuum is implemented as a SASKTRAN2 constituent in {py:class}`sasktran2_ext.continuum.MTCKDContinuum`,

The continuum constituent calls the MT_CKD fortran module, which requires pressure, temperature, H2O, CO2, and O3 as
input parameters. Therefore, to add the continuum to our atmosphere, we must also include these species.

We use the base SASKTRAN2 package to set up our model for a nadir-viewing infrared calculation,

```{code-cell}
:tags: ["remove-cell"]
import numpy as np
import sasktran2 as sk

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
```

```{code}
import numpy as np
import sasktran2 as sk

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
```

For comparison, calculate the radiance before adding the continuum,
```{code-cell}
engine = sk.Engine(config, geometry, viewing_geo)
output_no_continuum = engine.calculate_radiance(atmosphere)
```

Now add the continuum constituent from SASKTRAN2 extensions and redo the calculation,

```{code-cell}
from sasktran2_ext.continuum import MTCKDContinuum

atmosphere["continuum"] = MTCKDContinuum()

output_with_continuum = engine.calculate_radiance(atmosphere)
```

Plot the result,
```{code-cell}
:tags: ["remove-input","remove-stdout","remove-stderr"]
import matplotlib.pyplot as plt
output_no_continuum["radiance"].isel(los=0, stokes=0).plot(label="No continuum")
output_with_continuum["radiance"].isel(los=0, stokes=0).plot(label="With continuum")
plt.legend()
```
