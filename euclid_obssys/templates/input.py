################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology
import astropy.units as u
from scipy.integrate import quad
import importlib
import healpy as hp
import euclid_obssys.pozzettiLF as pozzettiLF
from astropy.io import fits
from os import path
import sys

# directory where products are stored, absolute path
repo = "REPODIR"
project = "PROJECTDIR"

# master catalog queried from Cosmohub
query = "8614"


# additional label in case of specific tests
otherlab = None

# tags in this list will be skipped in file names
skip_tags = ["alpha", "sm", "lf_model", "cm", "rsd", "ngrid", "Lbox"]

# FOOTPRINT
# Survey footprint, including an angular mask and a global redshift range
# it is read by the read_footprint() function
# None indicates to read the full Flagship octant
# footprint_tag = '100sqdeg'
# footprint_tag = None
footprint_tag = FOOTTAG #"BenSC8"


# PINOCCHIO
# this should give the path to the repository where a set of pinocchio runs is stored
# WITHOUT the final /
# None if pinocchio runs are not used
pinocchio_repo = "PINOCCHIOREPO"
# this should give the kernel for the run names, without numbers
pinocchio_kernel = "GeppettoFC"
pinocchio_rotator = hp.Rotator(
    rot=[-81.9608336, 22.91156899, 33.37698195]
)  # rot=[0,-26,-53] in equatorial coordinates
# Here we specify the range of pinocchio runs to be used by dndz.py and createRandom.py
# in case pinocchio_last_run is None, dndz and random are done as for the data catalog
# otherwise pinocchio_last_run MUST BE larger than pinocchio_first_run
# createSDHODfromPinocchio does not use these variables but requires the first and last
#   run to be provided on the command line, for greater flexibility
pinocchio_first_run = 1
pinocchio_last_run = 22
MASS_SHIFT = [0.2, 0.2]

#  GALAXY CATALOG SECTION

# LF model - may be '1' for Pozzetti model 1 Halpha galaxies
#                   '3' for Pozzetti model 3 Halpha galaxies
#                   '9' for photometric catalog limited to H<=24
#                   '0' for DM halo mass
lf_model = "3"
logflux_limit = np.log10(2.0e-16)  # This is the nominal flux limit
smoothing_length = 5  # Gaussian smoothing for the SDHOD and for dn/dz
cat_type = "CATTYPE"  # can be 'flagship', 'sdhod', 'pinocchio' or 'box'
shuffled_fluxes = SHUFFLE

deltazbin = 0.01  # this is the small binning used to compute number counts and measure the SDHOD
SEED_hod = (
    54321  # Seed for random number in the generation of galaxy catalog with SDHOD
)

# RANDOM CATALOG SECTION
SEED_random = 12345
alpha = 50
smooth_dndz_in_random = (
    False  # if True the random follows the smoothed dn/dz TO BE IMPLEMENTED
)
random_realization_number = None


# SELECTION SECTION
selection_data_tag = SELDATA
selection_random_tag = SELRAND 
apply_dataselection_to_random = True


# MEASUREMENT SECTION
# Redshift Space Distorsions - if True use observed_redshift_gal as redshift,
# else use true_redshift_gal
RSD = RSDFLAG


# redshift binning
# finalCatZShell = [(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3),(1.3,1.4),(1.4,1.5),(1.5,1.6),(1.6,1.7),(1.7,1.8)]
# finalCatZShell = [(1.0,1.1)]
finalCatZShell = [(0.9, 1.2)]  # ,(1.2,1.5),(1.5,1.8)]

# plotting and showing
PLOT = False
SHOW = False

# rotation applied to mollweide view plots - not a healpy rotator
rotmoll = [40, 55, -90]


# COSMOLOGY SECTION
LCDMmodel = FlatLambdaCDM(Om0=0.319, H0=67.0)
LCDMmodelLF = FlatLambdaCDM(Om0=0.300, H0=70.0)
Mpc_h = u.def_unit("Mpc_h", u.Mpc / LCDMmodel.h)
Mpc_h_LF = u.def_unit("Mpc_h_LF", u.Mpc / LCDMmodel.h)
cm_LF = u.def_unit("cm_LF", u.cm)
l_unit = Mpc_h
params = {
    "flat": True,
    "H0": 67.0,
    "Om0": 0.319,
    "Ob0": 0.049,
    "sigma8": 0.83,
    "ns": 0.96,
    "relspecies": False,
}
cosmology.addCosmology("FS", params)
cosmo = cosmology.setCosmology("FS")
cmrelation = "diemer19"


# LE3 CATALOGS SECTION
cat4le3_format = "fits"
WriteLE3Random = True
max_PKs_in_script = 50
max_2PCFs_in_script = 10
# PK parameters; with Lbox=None compute_Lbox is set to 'true'
Lbox = 2500.0
ngrid = 512


###############################################
# What is written below should not be changed #
###############################################

# check that repo and project end with /
if repo[-1] != "/":
    repo += "/"
if project[-1] != "/":
    project += "/"

# applying data selection to the construction of random
if (apply_dataselection_to_random) & (selection_data_tag is not None):
    selection_random_tag = "DS" + selection_data_tag

# sets all pinocchio variables to None if pinocchio runs are not used
if cat_type != "pinocchio":
    pinocchio_repo = None
    pinocchio_kernel = None
    pinocchio_rotator = None
    pinocchio_first_run = None
    pinocchio_last_run = None
    MASS_SHIFT = None

# sets the redshift key in the catalog
if RSD:
    redshift_key = "observed_redshift_gal"
else:
    redshift_key = "true_redshift_gal"


# this is for reading a footprint file
def read_footprint():

    if footprint_tag is None:
        footprint_fname = repo + "Footprints/FullOctant.fits"
    else:
        footprint_fname = repo + "Footprints/{}.fits".format(footprint_tag)

    print("# Loading footprint from {}...".format(footprint_fname))

    if not path.exists(footprint_fname):
        print("ERROR: footprint %s does not exist, please create it" % footprint_fname)
        return None

    ff = fits.open(footprint_fname)
    footprint_res = ff[1].header["RES"]
    footprint_zrange = [ff[1].header["MINZ"], ff[1].header["MAXZ"]]
    sky_fraction = ff[1].header["SKYFRAC"]

    return footprint_res, footprint_zrange, sky_fraction, ff[1].data["FOOTPRINT_G"]


SPEEDOFLIGHT = 299792.0  # Km/s
sqdegonthesky = (180 / np.pi) ** 2 * 4 * np.pi
flux_limit = 10.0 ** logflux_limit


def dummy(L, z):
    return 0


# Pozzetti model
dphi_dL = {
    "0": dummy,
    "1": pozzettiLF.model1Fast,
    "3": pozzettiLF.model3Fast,
    "9": dummy,
}[lf_model]
flux_key = {
    "0": "halo_lm",
    "1": "logf_halpha_model1_ext",
    "3": "logf_halpha_model3_ext",
    "9": "euclid_nisp_H",
}[lf_model]

if shuffled_fluxes:
    flux_key = "sh_" + flux_key

# Conversion from Flux to Luminosity
def fluxToLuminosity(z, flux=None):

    if flux is None:
        flux = flux_limit
    return (
        flux * 4.0 * np.pi * (LCDMmodelLF.luminosity_distance(z).to(u.cm)) ** 2
    ).value


# Integrated number of Galaxies per sq degree
def Pozzetti_dndz(z1, z2, minflux=None, maxflux=None):

    if maxflux is None:
        maxflux = 1.0e-13

    integrandL = lambda L, z: dphi_dL(L, z)
    integrandz = (
        lambda z: 4.0
        * np.pi
        / sqdegonthesky
        * LCDMmodelLF.differential_comoving_volume(z).value
        * quad(
            integrandL,
            fluxToLuminosity(z, minflux),
            fluxToLuminosity(z, maxflux),
            args=z,
        )[0]
    )

    return quad(integrandz, z1, z2)[0]
