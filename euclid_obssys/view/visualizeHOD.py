################################################################################
### Authors: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                  ###
###                                                                          ###
################################################################################
from . import register_view_tool
from ..config import readConfig


@register_view_tool
def visualizeHOD(config: str, ext: str = "png"):
    import numpy as np
    import matplotlib
    import os

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ..disk import DefaultCatalogRead
    from .. import filenames
    import sys

    input = readConfig(config)

    def divide(a, b):
        filtro = b > 0
        res = np.zeros_like(a)
        res[filtro] = a[filtro] / b[filtro]
        return res

    print(f"# Running visualizeHOD.py with {config}")

    c = np.array(
        [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
    )

    with DefaultCatalogRead(filenames.SDHOD(input)) as store:
        sdhod = store["sdhod"]

    izbins = np.arange(12) * 10
    Mx = 10.0 ** (0.5 * (sdhod["M_bins"][0, 1:] + sdhod["M_bins"][0, :-1]))
    zx = 0.5 * (sdhod["z_bins"][0, 1:] + sdhod["z_bins"][0, :-1])

    plt.figure()
    plt.title("SD-HOD, centrals")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([0.8e11, 1.0e15])
    plt.ylim([1.0e-5, 100])

    for ic in np.arange(10):
        z = zx[izbins[ic]]
        y = divide(
            sdhod["n_cen_gaus"][0, 0, izbins[ic]], sdhod["n_hal_gaus"][0, izbins[ic]]
        )
        ff = y > 0
        plt.plot(Mx[ff], y[ff] * z ** 8, label="%3.1f" % z, c=c[ic], linewidth=3)

        y = divide(sdhod["n_cen"][0, 0, izbins[ic]], sdhod["n_halos"][0, izbins[ic]])
        ff = y > 0
        plt.plot(Mx[ff], y[ff] * z ** 8, c=c[ic])

    plt.xlabel("halo mass (Msun/h)")
    plt.ylabel(r"$N_{\rm cen} (M_h|z)$")
    plt.legend()
    plt.savefig(filenames.plot_hod(input, type="centrals"))

    plt.figure()
    plt.title("SD-HOD, satellites")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim([1.0e-5, 1e3])

    for ic in np.arange(10):
        z = zx[izbins[ic]]
        y = divide(
            sdhod["n_sat_gaus"][0, 0, izbins[ic]], sdhod["n_hal_gaus"][0, izbins[ic]]
        )
        ff = y > 0
        plt.plot(Mx[ff], y[ff] * z ** 8, label="%3.1f" % z, c=c[ic], linewidth=3)

        y = divide(sdhod["n_sat"][0, 0, izbins[ic]], sdhod["n_halos"][0, izbins[ic]])
        ff = y > 0
        plt.plot(Mx[ff], y[ff] * z ** 8, c=c[ic])

    plt.xlabel("halo mass (Msun/h)")
    plt.ylabel(r"$N_{\rm cen} (M_h|z)$")
    plt.legend()
    plt.savefig(filenames.plot_hod(input, type="sats"))
