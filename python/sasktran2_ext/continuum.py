from __future__ import annotations

import numpy as np
import sasktran2 as sk
from sasktran2.constituent.base import Constituent
from sasktran2.util.interpolation import linear_interpolating_matrix

from sasktran2_ext import mt_ckd


class MTCKDContinuum(Constituent):
    def __init__(
        self,
        h2o_name: str = "H2O",
        co2_name: str = "CO2",
        o3_name: str = "O3",
        numeric_wf_fractional_change=1e-5,
        numeric_wf_central_difference=True,
    ):
        """
        The MT-CKD continuum absorption model.

        Parameters
        ----------
        h2o_name : str, optional
            The name of the H2O constituent in the atmosphere., by default "H2O"
        co2_name : str, optional
            The name of the CO2 constituent in the atmosphere., by default "CO2"
        o3_name : str, optional
            The name of the O3 constituent in the atmosphere., by default "O3"
        """
        self._h2o_name = h2o_name
        self._co2_name = co2_name
        self._o3_name = o3_name
        self._mtckd_wavenumbers = np.arange(
            -10, 19910, 10
        )  # same as the wavenumber grid hardcoded in the MTCKD fortran code
        self._fractional_change = numeric_wf_fractional_change
        self._central_difference = numeric_wf_central_difference

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the continuum constituent"
            raise ValueError(msg)

        if atmo.pressure_pa is None:
            msg = "It is required to set the pressure_pa property in the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo[self._h2o_name] is None:
            msg = f"It is required to add an {self._h2o_name} constituent to the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo[self._co2_name] is None:
            msg = f"It is required to add an {self._co2_name} constituent to the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo[self._o3_name] is None:
            msg = f"It is required to add an {self._o3_name} constituent to the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        wavenum_interp_matrix = linear_interpolating_matrix(
            self._mtckd_wavenumbers,
            atmo.wavenumbers_cminv,
            "zero",
        )
        h2o_alt_interp_matrix = linear_interpolating_matrix(
            atmo[self._h2o_name].altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        co2_alt_interp_matrix = linear_interpolating_matrix(
            atmo[self._co2_name].altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        o3_alt_interp_matrix = linear_interpolating_matrix(
            atmo[self._o3_name].altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )

        # interpolate vmrs to altitude grid
        h2o_vmr = h2o_alt_interp_matrix @ atmo[self._h2o_name].vmr
        co2_vmr = co2_alt_interp_matrix @ atmo[self._co2_name].vmr
        o3_vmr = o3_alt_interp_matrix @ atmo[self._o3_name].vmr

        # mt_ckd returns optical depth. set path length to 1.0 m (100.0 cm) so return value is equivalent to absorption in m^-1
        continuum_absorption = mt_ckd(
            atmo.pressure_pa,
            atmo.temperature_k,
            h2o_vmr,
            co2_vmr,
            o3_vmr,
            100.0,  # path length in cm
        )

        # remove unused portion of array
        continuum_absorption = continuum_absorption[:, 0 : len(self._mtckd_wavenumbers)]

        # interpolate continuum to atmosphere grid
        atmo.storage.total_extinction[:] += (
            np.nan_to_num(continuum_absorption) @ wavenum_interp_matrix.T
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        wavenum_interp_matrix = linear_interpolating_matrix(
            self._mtckd_wavenumbers,
            atmo.wavenumbers_cminv,
            "zero",
        )
        h2o_alt_interp_matrix = linear_interpolating_matrix(
            atmo[self._h2o_name].altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        co2_alt_interp_matrix = linear_interpolating_matrix(
            atmo[self._co2_name].altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        o3_alt_interp_matrix = linear_interpolating_matrix(
            atmo[self._o3_name].altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )

        p_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_pressure_pa")
        t_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_temperature_k")
        h2o_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_{self._h2o_name}_vmr"
        )
        co2_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_{self._co2_name}_vmr"
        )
        o3_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_{self._o3_name}_vmr"
        )

        pressure_pa = np.copy(atmo.pressure_pa)

        temperature_k = np.copy(atmo.temperature_k)

        alts = atmo.model_geometry.altitudes()

        h2o_vmr = h2o_alt_interp_matrix @ atmo[self._h2o_name].vmr
        co2_vmr = co2_alt_interp_matrix @ atmo[self._co2_name].vmr
        o3_vmr = o3_alt_interp_matrix @ atmo[self._o3_name].vmr

        base_cont = mt_ckd(
            pressure_pa,
            temperature_k,
            h2o_vmr,
            co2_vmr,
            o3_vmr,
            100.0,
        )[:, 0 : len(self._mtckd_wavenumbers)]

        for input_var_base, d_mapping in zip(
            [pressure_pa, temperature_k, h2o_vmr, co2_vmr, o3_vmr],
            [p_mapping, t_mapping, h2o_mapping, co2_mapping, o3_mapping],
            strict=True,
        ):
            input_var = input_var_base
            dx = input_var * self._fractional_change

            input_var += dx
            cont_above = mt_ckd(
                pressure_pa,
                temperature_k,
                h2o_vmr,
                co2_vmr,
                o3_vmr,
                100.0,
            )[:, 0 : len(self._mtckd_wavenumbers)]

            if self._central_difference:
                # central diff
                input_var -= 2 * dx
                cont_below = mt_ckd(
                    pressure_pa,
                    temperature_k,
                    h2o_vmr,
                    co2_vmr,
                    o3_vmr,
                    100.0,
                )[:, 0 : len(self._mtckd_wavenumbers)]
                input_var += dx

                central_diff_wf = (cont_above - cont_below) / (2 * dx[:, np.newaxis])
            else:
                # forward diff
                central_diff_wf = (cont_above - base_cont) / (dx[:, np.newaxis])

                input_var -= dx

            d_mapping.d_extinction[:] += (
                np.nan_to_num(central_diff_wf) @ wavenum_interp_matrix.T
            )

            d_ssa = (
                np.nan_to_num(central_diff_wf)
                @ wavenum_interp_matrix.T
                * (-atmo.storage.ssa)
                / atmo.storage.total_extinction
            )
            d_mapping.d_ssa[:] += d_ssa

            d_mapping.interp_dim = "altitude"
            d_mapping.interpolator = np.eye(len(alts))
