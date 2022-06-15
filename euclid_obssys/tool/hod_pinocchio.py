################################################################################
### Authors: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                  ###
###                                                                          ###
################################################################################
from . import register_tool, need_dask
from ..config import readConfig


chunks_lfsat=1000000
chunks_c=1000000

def process_concentrations(Mdelta, z2, cmrelation, cosmo):
    from colossus.halo import concentration
    from colossus.cosmology import cosmology

    cosmology.setCurrent(cosmo)
    return concentration.concentration(Mdelta, "200c", z=z2, model=cmrelation)


@register_tool
@need_dask
def createHODFromPinocchio(config: str, starting_run: int, last_run: int) -> None:
    """Create a HOD catalog from a set of Pinocchio runs.

    Args:
        config (str): Pipeline configuration file.
        starting_run (int): ID of the first Pinocchio run to consider
        last_run (int): ID of the last Pinocchio run to consider
    """

    from astropy.io import fits
    import numpy as np
    from .. import NFW, utils, sdhod, filenames, time
    from .. import match_pinocchio_masses as match
    from colossus.halo import concentration
    from colossus.cosmology import cosmology
    from scipy.stats import poisson
    from .. import pinocchio as rp
    from ..disk import DefaultCatalogWrite
    import healpy as hp
    import sys
    import os.path
    import numexpr as ne
    from dask import delayed
    import dask.array as da
    from dask.distributed import Client

    input = readConfig(config)

    print(f"# Running createSDHODfromPinocchio.py from Pinocchio with {config}")

    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()

    # this is needed to calibrate the 1-halo term in redshift space
    SCALE_VELOCITIES = 0.7

    # seed for random numbers
    # NB: this implies that to have the same result you need to process the same pinocchio mocks
    # to have a more flexible approach one should create an array of seeds, one for each run
    np.random.seed(seed=input.SEED_hod)

    print("# Reading the HOD table {}...".format(filenames.SDHOD(input)))
    hodtable = fits.getdata(filenames.SDHOD(input))
    hodtable_da = da.array(hodtable)

    for myrun in np.arange(starting_run, last_run + 1):

        # reading the PLC

        print("\n# -- starting with catalog %d --" % myrun)

        pinfname = filenames.pinplc(input, myrun)
        if not os.path.isfile(pinfname) and not os.path.isfile(pinfname + '.0'):
            print(f"ERROR, file {pinfname} not found, skipping this run")
            continue
        print("# Reading plc from file {}".format(pinfname))
        pincat = rp.plc(pinfname)
        print("Reading done")

        # filtering halos to reduce processing time

        print(
            "# Filtering the catalog for halos with true redshift in the range [{}-{})".format(
                footprint_zrange[0], footprint_zrange[1]
            )
        )
        z_filter = (pincat.redshift < footprint_zrange[1]) & (
            pincat.redshift >= footprint_zrange[0]
        )

        print("# Rotating coordinates to match the footprint")
        conv = np.pi / 180.0
        dec, ra = input.pinocchio_rotator(
            np.pi / 2.0 - pincat.theta[z_filter] * conv, pincat.phi[z_filter] * conv
        )

        z = pincat.redshift[z_filter]
        obsz = pincat.obsz[z_filter]
        Mass = pincat.Mass[z_filter]
        dist = np.sqrt(
            pincat.pos[z_filter, 0] ** 2
            + pincat.pos[z_filter, 1] ** 2
            + pincat.pos[z_filter, 2] ** 2
        )
        vlos = pincat.vlos[z_filter]
        name = pincat.name[z_filter]

        del pincat

        print("# Filtering the catalog for halos in the footprint")
        pix = hp.ang2pix(footprint_res, dec, ra)
        foot_filter = footprint[pix]==1
        dec = dec[foot_filter]
        ra = ra[foot_filter]
        z = z[foot_filter]
        Mass = Mass[foot_filter]
        dist = dist[foot_filter]
        vlos = vlos[foot_filter]
        name = name[foot_filter]

        dec = np.pi / 2.0 - dec

        Nhalos = z.size

        print("I will process {} halos".format(Nhalos))

        print("# Applying the HOD")
        # this implements the calibration based on clustering-matched or abundance-matched halos
        with time.check_time("HOD"):
            logM = match.logFSmasses(Mass, z)
            if input.MASS_SHIFT is not None:

                # number of central and satellite galaxies with clustering-matched masses
                mass_shift = input.MASS_SHIFT[0] + (1.0 - z) * input.MASS_SHIFT[1]
                factor = match.correct_number_density(mass_shift, z)
                Ncen = sdhod.Ncen(hodtable, logM - mass_shift, z) * factor
                Nsat_a = sdhod.Nsat(hodtable, logM - mass_shift, z) * factor
                print(f"# Nsat_a total = {Nsat_a.sum()}, factor = {factor}, mass_shift = {mass_shift}")

            else:

                # number of central and satellite galaxies with abundance-matched masses
                Ncen = sdhod.Ncen(hodtable, logM, z)
                Nsat_a = sdhod.Nsat(hodtable, logM, z)
                print(f"# (MASS_SHIFT none) Nsat_a total = {Nsat_a.sum()}")

        print("# Random sampling the Centrals")
        with time.check_time("Random sampling Centrals"):
            hostCentral = np.random.rand(Nhalos) <= Ncen
            totCentrals = hostCentral.sum()

        print("There will be {} central galaxies".format(totCentrals))

        print("# Random sampling the Satellites")
        with time.check_time("Random sampling satellites"):
            Nsat = poisson.rvs(Nsat_a)
            totSat = Nsat.sum()
            hostSat = Nsat > 0
            del Nsat_a

        print("# There will be {} satellite galaxies".format(totSat))

        # Total number of Galaxies Sat+Cen
        Ngal = int(totSat + totCentrals)

        print("# Creating the galaxy variables")
        ra_gal = np.empty(Ngal, dtype=np.float32)
        dec_gal = np.empty(Ngal, dtype=np.float32)
        xgal = np.empty(Ngal, dtype=np.float32)
        ygal = np.empty(Ngal, dtype=np.float32)
        zgal = np.empty(Ngal, dtype=np.float32)
        rgal = np.empty(Ngal, dtype=np.float32)
        vlosgal = np.empty(Ngal, dtype=np.float32)
        true_zgal = np.empty(Ngal, dtype=np.float32)
        obs_zgal = np.empty(Ngal, dtype=np.float32)
        halo_m = np.empty(Ngal, dtype=np.float32)
        log10f = np.empty(Ngal, dtype=np.float32)
        kind = np.empty(Ngal, dtype=np.int32)
        haloid = np.empty(Ngal, dtype=np.int64)

        cid = np.arange(Ngal)

        print("# Working on the {} Centrals".format(totCentrals))
        with time.check_time("Centrals"):
            with time.check_time("RA"):
                ra_gal[:totCentrals] = ra[hostCentral]
            with time.check_time("DEC"):
                dec_gal[:totCentrals] = dec[hostCentral]
            with time.check_time("dist"):
                rgal[:totCentrals] = dist[hostCentral]
            with time.check_time("cos_dec"):
                cos_dec = np.cos(dec_gal[:totCentrals])
            with time.check_time("xyzgal"):
                xgal[:totCentrals] = (
                    dist[hostCentral] * cos_dec * np.cos(ra_gal[:totCentrals])
                )
                ygal[:totCentrals] = (
                    dist[hostCentral] * cos_dec * np.sin(ra_gal[:totCentrals])
                )
                zgal[:totCentrals] = dist[hostCentral] * np.sin(dec_gal[:totCentrals])

            with time.check_time("leftovers"):
                vlosgal[:totCentrals] = vlos[hostCentral]
                true_zgal[:totCentrals] = z[hostCentral]
                halo_m[:totCentrals] = Mass[hostCentral]

            with time.check_time("mass shift"):
                if input.MASS_SHIFT is not None:
                    log10f[:totCentrals] = sdhod.lfcen(
                        hodtable,
                        logM[hostCentral] - mass_shift[hostCentral],
                        true_zgal[:totCentrals],
                    )
                else:
                    log10f[:totCentrals] = sdhod.lfcen(
                        hodtable, logM[hostCentral], true_zgal[:totCentrals]
                    )

        kind[:totCentrals] = 0
        haloid[:totCentrals] = name[hostCentral]

        print("# Working on the {} Satellites".format(totSat))

        print(
            "# Calculating concentrations according to {}...".format(input.cmrelation)
        )

        kind[totCentrals:] = 1

        with time.check_time("Satellite"):
            MDelta = Mass[hostSat]
            zhost = z[hostSat]

            Mdelta_da = da.array(MDelta).rechunk(chunks=chunks_c)
            zhost_da = da.array(zhost).rechunk(chunks=chunks_c)

            with time.check_time("Concentrations"):
                concentrations = da.map_blocks(
                    process_concentrations,
                    Mdelta_da,
                    zhost_da,
                    input.cmrelation,
                    cosmology.getCurrent(),
                    dtype=float,
                )

        RDelta = (3.0 * MDelta / 4.0 / np.pi / 200.0 / input.cosmo.rho_c(zhost)) ** (
            1.0 / 3.0
        )
        RDelta /= 1e3  # To Mpc/h

        print("# Distributing satellite galaxies in halos...")

        NwithSat = np.count_nonzero(hostSat)
        Nsat = Nsat[hostSat]

        myhalo = np.zeros(totSat, dtype=np.int32)
        cc = 0
        for i in range(NwithSat):
            myhalo[cc : cc + Nsat[i]] = i
            cc += Nsat[i]

        conv = np.pi / 180.0

        concentrations_da = concentrations[myhalo].rechunk(chunks=chunks_c)
        randr_da = da.map_blocks(
            lambda c, u: np.array(list(NFW.getr(c0, u0)[0] for c0, u0 in zip(c, u))),
            concentrations_da.rechunk(chunks=chunks_c),
            da.random.random(totSat).rechunk(chunks=chunks_c),
        )
        with time.check_time("randt"):
            randt, randp = utils.randomSpherePoint(totSat)
        with time.check_time("randv"):
            print(randr_da.shape, concentrations_da.shape)
            randv_da = da.map_blocks(
                lambda c, r: np.array(
                    list(utils.circularVelocity(c0, r0) for c0, r0 in zip(c, r))
                ),
                concentrations_da,
                randr_da,
            )
            print(randv_da.shape)
            randv_da = randv_da * da.sqrt(da.array(MDelta[myhalo] / RDelta[myhalo]).rechunk(chunks=chunks_c))
            print(randv_da.shape)

        with time.check_time("randr+randv"):
            randr = randr_da.compute()
            randv = randv_da.compute()

        print("# random numbers done")
        print("# updating spatial positions...")

        sint = np.sin(randt)
        cosd = np.cos(dec[hostSat][myhalo])
        xgal[totCentrals:] = dist[hostSat][myhalo] * cosd * np.cos(
            ra[hostSat][myhalo]
        ) + RDelta[myhalo] * randr * sint * np.cos(randp)
        ygal[totCentrals:] = dist[hostSat][myhalo] * cosd * np.sin(
            ra[hostSat][myhalo]
        ) + RDelta[myhalo] * randr * sint * np.sin(randp)
        zgal[totCentrals:] = dist[hostSat][myhalo] * np.sin(
            dec[hostSat][myhalo]
        ) + RDelta[myhalo] * randr * np.cos(randt)
        print("# updated spatial positions done")

        rgal[totCentrals:] = np.sqrt(
            xgal[totCentrals:] ** 2 + ygal[totCentrals:] ** 2 + zgal[totCentrals:] ** 2
        )

        print("# updating angular positions...")
        dec_gal[totCentrals:] = np.arcsin(zgal[totCentrals:] / rgal[totCentrals:])
        ra_gal[totCentrals:] = np.arctan2(ygal[totCentrals:], xgal[totCentrals:])
        print("# updated angular positions done")

        v = SCALE_VELOCITIES * np.array(
            [
                randv * np.random.normal(size=totSat),
                randv * np.random.normal(size=totSat),
                randv * np.random.normal(size=totSat),
            ]
        )

        vlosgal[totCentrals:] = (
            vlos[hostSat][myhalo]
            + (
                v[0, :] * xgal[totCentrals:]
                + v[1, :] * ygal[totCentrals:]
                + v[2, :] * zgal[totCentrals:]
            )
            / rgal[totCentrals:]
        )

        print("# updated velocities and redshift done")

        haloid[totCentrals:] = name[hostSat][myhalo]
        halo_m[totCentrals:] = Mass[hostSat][myhalo]
        true_zgal[totCentrals:] = z[hostSat][myhalo]

        with time.check_time("log10f"):
            if input.MASS_SHIFT is not None:
                log10f[totCentrals:] = da.map_blocks(
                    lambda lM, ms, zg: (sdhod.lfsat(hodtable_da.compute(), lM - ms, zg)),
                    da.array(logM[hostSat][myhalo]).rechunk(chunks=chunks_lfsat),
                    da.array(mass_shift[hostSat][myhalo]).rechunk(chunks=chunks_lfsat),
                    da.array(true_zgal[totCentrals:]).rechunk(chunks=chunks_lfsat),
                    dtype=np.float,
                ).compute()
            else:
                log10f[totCentrals:] = da.map_blocks(
                    lambda lM, zg: (sdhod.lfsat(hodtable_da.compute(), lM, zg)),
                    da.array(logM[hostSat][myhalo]).rechunk(chunks=chunks_lfsat),
                    da.array(true_zgal[totCentrals:]).rechunk(chunks=chunks_lfsat),
                    dtype=np.float,
                ).compute()

        halo_m = np.log10(halo_m)
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

        fname = filenames.pincat(input, myrun)
        print("# Saving the catalog to file {}".format(fname))

        with DefaultCatalogWrite(fname) as store:
            catalog = store.new_array(
                "catalog",
                shape=(Ngal,),
                dtype=[
                    ("x_gal", np.float32),
                    ("y_gal", np.float32),
                    ("z_gal", np.float32),
                    ("ra_gal", np.float32),
                    ("dec_gal", np.float32),
                    ("kind", np.int32),
                    ("true_redshift_gal", np.float32),
                    ("observed_redshift_gal", np.float32),
                    ("halo_lm", np.float32),
                    ("id", np.int32),
                    ("halo_id", np.int32),
                    (input.flux_key, np.float32),
                    ("sh_" + input.flux_key, np.float32),
                ],
            )

            catalog["x_gal"] = xgal
            catalog["y_gal"] = ygal
            catalog["z_gal"] = zgal
            catalog["true_redshift_gal"] = true_zgal
            catalog["observed_redshift_gal"] = obs_zgal
            catalog["ra_gal"] = ra_gal
            catalog["dec_gal"] = dec_gal
            catalog["halo_lm"] = halo_m
            catalog["kind"] = kind
            catalog[input.flux_key] = log10f
            catalog["sh_" + input.flux_key] = shuffled_log10f
            catalog["id"] = cid
            catalog["halo_id"] = haloid

        print(f"# -- done with catalog {myrun} --")

    print("# !!DONE!!")
