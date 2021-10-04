import numpy as np
from astropy.io import fits

"""
Given a measured SDHOD, this library provides the 
Ncen(log halo mass, redshift) 
and 
Nsat(log halo mass, redshift)
functions, together with functions lcen and lsat
that give the luminosity of a central or satellite galaxy.

"""


def divide(a, b):
    filtro = b > 0
    res = np.zeros_like(a)
    res[filtro] = a[filtro] / b[filtro]
    return res


# Fraction Number of Centrals
def Ncen(sdhod, logM, z):

    mass_index = np.digitize(logM, sdhod["M_bins"][0]) - 1
    z_index = np.digitize(z, sdhod["z_bins"][0]) - 1

    # Assuming the indexes are arrays
    try:

        # Reassignement of the out of bound indexes
        mass_index[mass_index == sdhod["M_bins"][0].size - 1] = (
            sdhod["M_bins"][0].size - 2
        )
        z_index[z_index == sdhod["z_bins"][0].size - 1] = sdhod["z_bins"][0].size - 2

        return divide(
            sdhod["n_cen_gaus"][0, 0, z_index, mass_index],
            sdhod["n_hal_gaus"][0, z_index, mass_index],
        )

    # Raised exception in case the object is not indexable
    except TypeError:

        # Reassignement of the out of bound indexes
        if z_index == sdhod["z_bins"][0].size - 1:

            z_index -= 1

        elif mass_index == sdhod["M_bins"][0].size - 1:

            mass_index -= 1

        return divide(
            sdhod["n_cen_gaus"][0, 0, z_index, mass_index],
            sdhod["n_hal_gaus"][0, z_index, mass_index],
        )


# log10 of the flux of the given central
def lfcen(sdhod, logM, z):

    mass_index = np.digitize(logM, sdhod["M_bins"][0]) - 1
    z_index = np.digitize(z, sdhod["z_bins"][0]) - 1

    try:
        # Reassignement of the out of bound indexes
        mass_index[mass_index == sdhod["M_bins"][0].size - 1] = (
            sdhod["M_bins"][0].size - 2
        )
        z_index[z_index == sdhod["z_bins"][0].size - 1] = sdhod["z_bins"][0].size - 2
        # Getting the corresponding Flux inside the Inverse CDF
        u = np.random.random(logM.size)
        f = [
            np.interp(
                ui,
                1
                - sdhod["n_cen_gaus"][0, :, zi, mi] / sdhod["n_cen_gaus"][0, 0, zi, mi],
                sdhod["f_bins"][0],
            )
            for ui, zi, mi in zip(u, z_index, mass_index)
        ]
        f = np.array(f)

    except TypeError:

        # Reassignement of the out of bound indexes
        if z_index == sdhod["z_bins"][0].size - 1:

            z_index -= 1

        elif mass_index == sdhod["M_bins"][0].size - 1:

            mass_index -= 1

        u = np.random.random()
        f = np.interp(
            u,
            1
            - sdhod["n_cen_gaus"][0, :, z_index, mass_index]
            / sdhod["n_cen_gaus"][0, 0, z_index, mass_index],
            sdhod["f_bins"][0],
        )

    return f


# Number of Satelites
def Nsat(sdhod, logM, z):

    mass_index = np.digitize(logM, sdhod["M_bins"][0]) - 1
    z_index = np.digitize(z, sdhod["z_bins"][0]) - 1

    # Assuming the indexes are arrays
    try:

        # Reassignement of the out of bound indexes
        mass_index[mass_index == sdhod["M_bins"][0].size - 1] = (
            sdhod["M_bins"][0].size - 2
        )
        z_index[z_index == sdhod["z_bins"][0].size - 1] = sdhod["z_bins"][0].size - 2

        return divide(
            sdhod["n_sat_gaus"][0, 0, z_index, mass_index],
            sdhod["n_hal_gaus"][0, z_index, mass_index],
        )

    except TypeError:

        # Reassignement of the out of bound indexes
        if z_index == sdhod["z_bins"][0].size - 1:

            z_index -= 1

        elif mass_index == sdhod["M_bins"][0].size - 1:

            mass_index -= 1

        return divide(
            sdhod["n_sat_gaus"][0, 0, z_index, mass_index],
            sdhod["n_hal_gaus"][0, z_index, mass_index],
        )


# log10 of the flux of the satelites
def lfsat(sdhod, logM, z):

    mass_index = np.digitize(logM, sdhod["M_bins"][0]) - 1
    z_index = np.digitize(z, sdhod["z_bins"][0]) - 1

    try:

        # Reassignement of the out of bound indexes
        mass_index[mass_index == sdhod["M_bins"][0].size - 1] = (
            sdhod["M_bins"][0].size - 2
        )
        z_index[z_index == sdhod["z_bins"][0].size - 1] = sdhod["z_bins"][0].size - 2
        # Getting the corresponding Flux inside the Inverse CDF
        u = np.random.random(logM.size)
        f = [
            np.interp(
                ui,
                1
                - sdhod["n_sat_gaus"][0, :, zi, mi] / sdhod["n_sat_gaus"][0, 0, zi, mi],
                sdhod["f_bins"][0],
            )
            for ui, zi, mi in zip(u, z_index, mass_index)
        ]

    except TypeError:

        # Reassignement of the out of bound indexes
        if z_index == sdhod["z_bins"][0].size - 1:

            z_index -= 1

        elif mass_index == sdhod["M_bins"][0].size - 1:

            mass_index -= 1

        u = np.random.random()
        f = np.interp(
            u,
            1.0
            - sdhod["n_sat_gaus"][0, :, z_index, mass_index]
            / sdhod["n_sat_gaus"][0, 0, z_index, mass_index],
            sdhod["f_bins"][0],
        )

    return f
