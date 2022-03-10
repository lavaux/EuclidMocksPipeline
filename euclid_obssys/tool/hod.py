################################################################################
### Authors: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig


@register_tool
def createSmoothHOD(config: str) -> None:
    """Create a HOD curves to prepare as SDHOD catalog.

    Based on a master file, it measures the HOD curves and smooths them in redshift.

    Args:
        config (str): Pipeline configuration file.
    """
    from .. import sdhod, NFW, utils, filenames
    import numpy as np
    from colossus.cosmology import cosmology
    from euclid_obssys.disk import DefaultCatalogRead, DefaultCatalogWrite
    from scipy.ndimage import gaussian_filter
    from scipy.stats import poisson
    import healpy as hp

    input = readConfig(config)

    print(f"# Running createSmoothHOD.py with {config}")

    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
    print(f"# The footprint covers {100*sky_fraction}% of the sky")

    # Import the Raw Catalog
    fname = filenames.indices(input)
    print("# Reading indices from file {}...".format(fname))
    with DefaultCatalogRead(fname) as cat_file:
        rawcat_indices = cat_file["indices"]["indices"]

    master_fname = filenames.master(input)
    print(f"# Reading Catalog {master_fname}...")
    with DefaultCatalogRead(master_fname) as cat_file:
        rawcat = cat_file["catalog"]

    print("# Extracting colums from the catalog...")
    good = np.array(rawcat_indices) > 0
    zs = rawcat["true_redshift_gal"][rawcat_indices][good]

    del rawcat_indices

    flux = rawcat[input.flux_key][good]
    lm = rawcat["halo_lm"][good]
    kind = rawcat["kind"][good]

    del rawcat
    del good

    print("# Binning galaxies in flux...")

    # Flux cuts from the minimum flux defined in input.py to
    # two magnitude higher
    Ncut = 400
    DeltaFlux = 1.6
    fcut = np.linspace(input.logflux_limit, input.logflux_limit + DeltaFlux, Ncut)

    # Selecting the galaxies inside the fcut bins
    gal_fluxbin = np.digitize(flux, fcut)

    # Redshift bins
    z1 = footprint_zrange[0]
    z2 = footprint_zrange[1]
    zbins = np.linspace(z1, z2, round((z2 - z1) / input.deltazbin + 1))

    # Volume bins
    Vbins = (
        (4 * np.pi / 3)
        * (
            input.cosmo.comovingDistance(0.0, zbins[1:]) ** 3
            - input.cosmo.comovingDistance(0.0, zbins[:-1]) ** 3
        )
        * sky_fraction
    )

    # Halo mass Bins
    Mbins = np.linspace(lm.min(), lm.max())

    # Number of Centrals and Satelites inside each bin of redshift and flux
    Ncen = np.zeros((fcut.size, zbins.size - 1, Mbins.size - 1))
    Nsat = np.zeros((fcut.size, zbins.size - 1, Mbins.size - 1))
    Nhalos = np.zeros((zbins.size - 1, Mbins.size - 1))
    Ncen_g = np.zeros((fcut.size, zbins.size - 1, Mbins.size - 1))
    Nsat_g = np.zeros((fcut.size, zbins.size - 1, Mbins.size - 1))
    Nhal_g = np.zeros((zbins.size - 1, Mbins.size - 1))

    print("# Selecting centrals and satellites...")

    cen = kind == 0

    print("# Binning galaxies in mass...")
    Nhalos, bins1, bins2 = np.histogram2d(zs[cen], lm[cen], bins=[zbins, Mbins])
    Nhalos = (Nhalos.T / Vbins).T

    print("# Constructing the SDHOD...")

    # Working on the flux limited sample
    for i in range(1, fcut.size + 1):

        print("    Working on flux cut {} of {}...".format(i, Ncut))

        # Reading kind array
        idsi = gal_fluxbin == i
        this_kind = kind[idsi]
        this_zs = zs[idsi]
        this_lm = lm[idsi]
        this_cen = this_kind == 0
        this_sat = ~this_cen

        Ncen[i - 1], bins1, bins2 = np.histogram2d(
            this_zs[this_cen], this_lm[this_cen], bins=[zbins, Mbins]
        )
        Nsat[i - 1], bins1, bins2 = np.histogram2d(
            this_zs[this_sat], this_lm[this_sat], bins=[zbins, Mbins]
        )

        Ncen[i - 1] = ((Ncen[i - 1].T) / Vbins).T
        Nsat[i - 1] = ((Nsat[i - 1].T) / Vbins).T

        # Smoothing the dn/dz
        for j in range(Mbins.size - 1):

            Ncen_g[i - 1, :, j] = gaussian_filter(
                Ncen[i - 1, :, j], input.smoothing_length
            )
            Nsat_g[i - 1, :, j] = gaussian_filter(
                Nsat[i - 1, :, j], input.smoothing_length
            )

    print("# Smoothing the halo SDHOD...")

    # Smoothing the halo dn/dz
    for j in range(Mbins.size - 1):
        Nhal_g[:, j] = gaussian_filter(Nhalos[:, j], input.smoothing_length)

    print("# Cumulating the number of galaxies...")

    # Cumulative Summation on the fluxes
    Ncen = Ncen.sum(axis=0) - Ncen.cumsum(axis=0) + Ncen
    Nsat = Nsat.sum(axis=0) - Nsat.cumsum(axis=0) + Nsat
    Ncen_g = Ncen_g.sum(axis=0) - Ncen_g.cumsum(axis=0) + Ncen_g
    Nsat_g = Nsat_g.sum(axis=0) - Nsat_g.cumsum(axis=0) + Nsat_g

    print("# Writing file {}...".format(filenames.SDHOD(input)))

    with DefaultCatalogWrite(filenames.SDHOD(input)) as out_file:

        cat = out_file.new_array(
            "sdhod",
            shape=(1,),
            dtype=[
                ("n_halos", float, (zbins.size - 1, Mbins.size - 1)),
                ("n_cen", float, (fcut.size, zbins.size - 1, Mbins.size - 1)),
                ("n_sat", float, (fcut.size, zbins.size - 1, Mbins.size - 1)),
                ("n_hal_gaus", float, (zbins.size - 1, Mbins.size - 1)),
                ("n_cen_gaus", float, (fcut.size, zbins.size - 1, Mbins.size - 1)),
                ("n_sat_gaus", float, (fcut.size, zbins.size - 1, Mbins.size - 1)),
                ("z_bins", float, zbins.size),
                ("V_bin", float, Vbins.size),
                ("M_bins", float, Mbins.size),
                ("f_bins", float, fcut.size),
            ],
        )

        cat["n_halos"] = Nhalos
        cat["n_cen"] = Ncen
        cat["n_sat"] = Nsat
        cat["n_hal_gaus"] = Nhal_g
        cat["n_cen_gaus"] = Ncen_g
        cat["n_sat_gaus"] = Nsat_g
        cat["z_bins"] = zbins
        cat["V_bin"] = Vbins
        cat["M_bins"] = Mbins
        cat["f_bins"] = fcut

    print("# done!")


###############################################################################


@register_tool
def createSDHOD_Catalog(config: str) -> None:
    """Creates a galaxy catalog by
    populating the halos of the master catalog with galaxies, using the
    SDHOD.

    Args:
        config (str): [description]

    """
    from .. import sdhod, NFW, utils, filenames
    from colossus.halo import concentration
    from colossus.cosmology import cosmology
    from astropy.io import fits
    from scipy.ndimage import gaussian_filter
    from scipy.stats import poisson
    import healpy as hp
    import numpy as np
    import numba as nb
    from euclid_obssys.disk import DefaultCatalogRead, DefaultCatalogWrite
    from dask import delayed
    from dask.distributed import Client

    input = readConfig(config)

    print("# Running createSDHOD_Catalog.py with {}".format(config))

    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
    del footprint

    # this is needed to calibrate the 1-halo term in redshift space
    SCALE_VELOCITIES = 0.7

    # seed for random numbers
    np.random.seed(seed=input.SEED_hod)

    print("# Reading the HOD table {}...".format(filenames.SDHOD(input)))
    with DefaultCatalogRead(filenames.SDHOD(input)) as in_file:
        hodtable = in_file["sdhod"]

    print("# Starting Dask...")
    client = Client()

    print("# Reading the halo catalog from {}...".format(filenames.master(input)))
    with DefaultCatalogRead(filenames.master(input)) as store:
        rawcat = store["catalog"]

        c_kind = rawcat["kind"]
        c_z = rawcat["true_redshift_gal"]
        c_logm = rawcat["halo_lm"]
        c_ra_gal = rawcat["ra_gal"]
        c_dec_gal = rawcat["dec_gal"]
        c_x_gal = rawcat["x_gal"]
        c_y_gal = rawcat["y_gal"]
        c_z_gal = rawcat["z_gal"]
        c_vrad_gal = rawcat["vrad_gal"]
        c_halo_id = rawcat["halo_id"]

        del rawcat

    print("# Filtering the catalog for main halos...")
    filt = c_kind == 0
    del c_kind

    c_z = c_z[filt]
    c_logm = c_logm[filt]
    c_ra_gal = c_ra_gal[filt]
    c_dec_gal = c_dec_gal[filt]
    c_x_gal = c_x_gal[filt]
    c_y_gal = c_y_gal[filt]
    c_z_gal = c_z_gal[filt]
    c_vrad_gal = c_vrad_gal[filt]
    c_halo_id = c_halo_id[filt]

    Nhalos = c_z.size

    print(f"# Found {Nhalos} main halos")

    print("# Random sampling the Centrals...")
    NCen = sdhod.Ncen(hodtable, c_logm, c_z)
    # The halo hosts an Halpha emitter if a random variable is less than NCen
    hostCentral = np.random.rand(Nhalos) <= NCen
    totCentrals = hostCentral.sum()

    print(f"# There will be {totCentrals} central galaxies")

    print("# Random sampling the Satellites...")
    # Getting the number of Satellites according to a poisson distribution
    Nsat = poisson.rvs(sdhod.Nsat(hodtable, c_logm, c_z))
    hostSat = Nsat > 0
    totSat = Nsat.sum()

    print(f"# There will be {totSat} satellite galaxies")

    # Total number of Galaxies Sat+Cen
    Ngal = int(totSat + totCentrals)

    Mass = 10.0 ** c_logm

    print("# Creating the galaxy variables...")
    ra_gal = np.empty(Ngal, dtype=np.float)
    dec_gal = np.empty(Ngal, dtype=np.float)
    xgal = np.empty(Ngal, dtype=np.float)
    ygal = np.empty(Ngal, dtype=np.float)
    zgal = np.empty(Ngal, dtype=np.float)
    rgal = np.empty(Ngal, dtype=np.float)
    vlosgal = np.empty(Ngal, dtype=np.float)
    true_zgal = np.empty(Ngal, dtype=np.float)
    obs_zgal = np.empty(Ngal, dtype=np.float)
    halo_m = np.empty(Ngal, dtype=np.float)
    log10f = np.empty(Ngal, dtype=np.float)
    kind = np.empty(Ngal, dtype=np.int)
    haloid = np.empty(Ngal, dtype=np.int64)

    cid = np.arange(Ngal)

    print("# Working on the {} centrals...".format(totCentrals))

    ra_gal[:totCentrals] = c_ra_gal[hostCentral]
    dec_gal[:totCentrals] = c_dec_gal[hostCentral]
    xgal[:totCentrals] = c_x_gal[hostCentral]
    ygal[:totCentrals] = c_y_gal[hostCentral]
    zgal[:totCentrals] = c_z_gal[hostCentral]
    rgal[:totCentrals] = np.sqrt(
        xgal[:totCentrals] ** 2 + ygal[:totCentrals] ** 2 + zgal[:totCentrals] ** 2
    )
    vlosgal[:totCentrals] = c_vrad_gal[hostCentral]
    true_zgal[:totCentrals] = c_z[hostCentral]
    halo_m[:totCentrals] = Mass[hostCentral]
    kind[:totCentrals] = 0
    haloid[:totCentrals] = c_halo_id[hostCentral]

    print("# Assigning fluxes to the centrals...")

    log10f[:totCentrals] = sdhod.lfcen(hodtable, c_logm[hostCentral], c_z[hostCentral])

    print("# Working on the {} Satellites".format(totSat))

    print("# Calculating concentrations according to {}...".format(input.cmrelation))

    kind[totCentrals:] = 1

    MDelta = Mass[hostSat].astype(np.float64)
    z2 = c_z[hostSat].astype(np.float64)
    model_name = input.cmrelation

    #@delayed
    def process_concentrations(Mdelta, z2):
       return concentration.concentration(Mdelta, "200c", z=z2, model=model_name)

    MDelta_da = da.array(MDelta).rechunk()
    z2_da = da.array(z2).rechunk()

    concentrations = da.map_blocks(process_concentrations, Mdelta_da, z2_da, dtype=float).compute()

    RDelta = (3.0 * MDelta / 4.0 / np.pi / 200.0 / input.cosmo.rho_c(z2)) ** (1.0 / 3.0)
    RDelta /= 1e3  # To Mpc/h

    print("# Distributing satellite galaxies in halos...")

    NwithSat = np.count_nonzero(hostSat)
    Nsat = Nsat[hostSat]

    myhalo = np.zeros(totSat, dtype=np.int)
    cc = 0
    for i in range(NwithSat):
        myhalo[cc : cc + Nsat[i]] = i
        cc += Nsat[i]

    conv = np.pi / 180.0

    randr = np.array(
        [
            NFW.getr(c, u)[0]
            for c, u in zip(concentrations[myhalo], np.random.rand(totSat))
        ]
    )
    randt, randp = utils.randomSpherePoint(totSat)
    randv = np.array(
        [utils.circularVelocity(c, r) for c, r in zip(concentrations[myhalo], randr)]
    ) * np.sqrt(MDelta[myhalo] / RDelta[myhalo])

    print("# random numbers done")
    print("# updating spatial positions...")

    dist = np.sqrt(
        c_x_gal[hostSat] ** 2 + c_y_gal[hostSat] ** 2 + c_z_gal[hostSat] ** 2
    )
    sint = np.sin(randt)
    cosd = np.cos(c_dec_gal[hostSat][myhalo] * conv)
    xgal[totCentrals:] = dist[myhalo] * cosd * np.cos(
        c_ra_gal[hostSat][myhalo] * conv
    ) + RDelta[myhalo] * randr * sint * np.cos(randp)
    ygal[totCentrals:] = dist[myhalo] * cosd * np.sin(
        c_ra_gal[hostSat][myhalo] * conv
    ) + RDelta[myhalo] * randr * sint * np.sin(randp)
    zgal[totCentrals:] = dist[myhalo] * np.sin(
        c_dec_gal[hostSat][myhalo] * conv
    ) + RDelta[myhalo] * randr * np.cos(randt)

    rgal[totCentrals:] = np.sqrt(
        xgal[totCentrals:] ** 2 + ygal[totCentrals:] ** 2 + zgal[totCentrals:] ** 2
    )

    print("# updating angular positions...")
    dec_gal[totCentrals:] = np.arcsin(zgal[totCentrals:] / rgal[totCentrals:])
    ra_gal[totCentrals:] = np.arctan2(ygal[totCentrals:], xgal[totCentrals:])

    # v = np.array( [ randv * np.random.rand(totSat), randv * np.random.rand(totSat), randv * np.random.rand(totSat) ] )
    v = SCALE_VELOCITIES * np.array(
        [
            randv * np.random.normal(size=totSat),
            randv * np.random.normal(size=totSat),
            randv * np.random.normal(size=totSat),
        ]
    )

    vlosgal[totCentrals:] = (
        c_vrad_gal[hostSat][myhalo]
        + (
            v[0, :] * xgal[totCentrals:]
            + v[1, :] * ygal[totCentrals:]
            + v[2, :] * zgal[totCentrals:]
        )
        / rgal[totCentrals:]
    )

    haloid[totCentrals:] = c_halo_id[hostSat][myhalo]
    halo_m[totCentrals:] = Mass[hostSat][myhalo]
    true_zgal[totCentrals:] = c_z[hostSat][myhalo]
    log10f[totCentrals:] = sdhod.lfsat(
        hodtable, c_logm[hostSat][myhalo], true_zgal[totCentrals:]
    )

    obs_zgal = true_zgal + vlosgal / input.SPEEDOFLIGHT * (1.0 + true_zgal)

    ra_gal *= 180.0 / np.pi
    ra_gal[ra_gal < 0] += 360.0
    dec_gal *= 180.0 / np.pi

    print("# Shuffling galaxy fluxes...")
    shuffled_log10f = np.copy(log10f)

    zbins = np.linspace(
        footprint_zrange[0],
        footprint_zrange[1],
        round((footprint_zrange[1] - footprint_zrange[0]) / input.deltazbin + 1),
    )
    zindex = np.digitize(true_zgal, zbins)
    for iz in np.arange(zbins.size - 1):
        ff = zindex == iz
        extract = shuffled_log10f[ff]
        np.random.shuffle(extract)
        shuffled_log10f[ff] = extract

    print("# Saving the catalog to file {}".format(filenames.hodcat(input)))
    with DefaultCatalogWrite(filenames.hodcat(input)) as output:
        catalog = output.new_array(
            "catalog",
            shape=(Ngal,),
            dtype=[
                ("x_gal", float),
                ("y_gal", float),
                ("z_gal", float),
                ("ra_gal", float),
                ("dec_gal", float),
                ("kind", int),
                ("true_redshift_gal", float),
                ("observed_redshift_gal", float),
                ("halo_lm", float),
                ("id", int),
                ("halo_id", int),
                (input.flux_key, float),
                ("sh_" + input.flux_key, float),
            ],
        )

        catalog["x_gal"] = xgal
        catalog["y_gal"] = ygal
        catalog["z_gal"] = zgal
        catalog["true_redshift_gal"] = true_zgal
        catalog["observed_redshift_gal"] = obs_zgal
        catalog["ra_gal"] = np.arctan2(ygal, xgal) * 180.0 / np.pi
        catalog["dec_gal"] = (
            90.0
            - np.arccos(zgal / (np.sqrt(xgal ** 2 + ygal ** 2 + zgal ** 2)))
            * 180.0
            / np.pi
        )
        catalog["halo_lm"] = np.log10(halo_m)
        catalog["kind"] = kind
        catalog[input.flux_key] = log10f
        catalog["sh_" + input.flux_key] = shuffled_log10f
        catalog["id"] = cid
        catalog["halo_id"] = haloid

    print("# !!DONE!!")
