################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology
import astropy.units as u
#from healpy.rotator import Rotator
from scipy.integrate import quad
import importlib
import healpy as hp
import euclid_obssys.pozzettiLF as pozzettiLF
from astropy.io import fits
from os import path

# directory where products are stored, absolute path
outdir = "OUTDIR"

# master catalog queried from Cosmohub
query = '8614'

# additional label in case of specific tests
otherlab=None

# FOOTPRINT
# Survey footprint, including an angular mask and a global redshift range
# it is read by the read_footprint() function
# None indicates to read the full Flagship octant
#footprint_tag = '100sqdeg'
footprint_tag = FOOTTAG

# GALAXY CATALOG SECTION

# LF model - may be '1' for Pozzetti model 1
#                   '3' for Pozzetti model 3
#                   '9' for photometric catalog limited to H<=24
#                   '0' for halo mass
lf_model='LFMODEL'
logflux_limit = np.log10(2.e-16) # This is the nominal flux limit
smoothing_length = 5             # Gaussian smoothing for the SDHOD and for dn/dz
cat_type = 'CATTYPE'               # can be 'flagship', 'sdhod', 'pinocchio' or 'box'
shuffled_fluxes=SHUFFLE

deltazbin = 0.01        # this is the small binning used to compute number counts and measure the SDHOD
SEED_hod=54321          # Seed for random number in the generation of galaxy catalog with SDHOD

# SELECTION SECTION
selection_data_tag = SELDATA #'MWext'
selection_random_tag = SELRAND #'MWext'

# RANDOM CATALOG SECTION
SEED_random=12345
alpha=50
smooth_dndz_in_random = False   # if True the random follows the smoothed dn/dz

# MEASUREMENT SECTION
# Redshift Space Distorsions - if True use observed_redshift_gal as redshift, 
# else use true_redshift_gal
RSD = RSDFLAG

# redshift binning
finalCatZShell = [(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3),(1.3,1.4),(1.4,1.5),(1.5,1.6),(1.6,1.7),(1.7,1.8)]
#finalCatZShell = [(1.0,1.1),(1.6,1.7)]

# plotting and showing
PLOT = False
SHOW = False
# rotation applied to mollweide view plots - not a healpy rotator
rotmoll=[0,0,0]

# COSMOLOGY SECTION
LCDMmodel    = FlatLambdaCDM(Om0 = 0.319, H0 = 67.0)
LCDMmodelLF  = FlatLambdaCDM(Om0 = 0.300, H0 = 70.0)
Mpc_h    = u.def_unit('Mpc_h', u.Mpc/LCDMmodel.h)
Mpc_h_LF = u.def_unit('Mpc_h_LF', u.Mpc/LCDMmodel.h)
cm_LF    = u.def_unit('cm_LF', u.cm)
l_unit   = Mpc_h
params = {'flat': True, 'H0': 67.0, 'Om0': 0.319, 'Ob0': 0.049, 'sigma8': 0.83, 'ns': 0.96, 'relspecies':False}
cosmology.addCosmology('FS', params)
cosmo = cosmology.setCosmology('FS')
cmrelation = "diemer19"

# LE3 CATALOGS SECTION
cat4le3_format="fits"
WriteLE3Random=True
max_PKs_in_script=15
# PK parameters; with Lbox=None compute_Lbox is set to 'true'
Lbox = None
ngrid = 1024

###############################################
# What is written below should not be changed #
###############################################


# this is for reading a footprint file
def read_footprint():
    if footprint_tag is None:
        footprint_fname = path.join(outdir,'Footprints','FullOctant.fits')
    else:
        footprint_fname = path.join(outdir,'Footprints','{}.fits'.format(footprint_tag))

    print("# Loading footprint from {}...".format(footprint_fname))

    if not path.exists(footprint_fname):
        print("ERROR: footprint %s does not exist, please create it"%footprint_fname)
        return None

    ff=fits.open(footprint_fname)
    footprint_res    = ff[0].header['RES']
    footprint_zrange = [ ff[0].header['MINZ'], ff[0].header['MAXZ'] ]
    sky_fraction     = ff[0].header['SKYFRAC']

    print(ff)
    return footprint_res, footprint_zrange, sky_fraction, ff[1].data['FOOTPRINT_G']

SPEEDOFLIGHT = 299792.0 # Km/s
sqdegonthesky = (180/np.pi)**2*4*np.pi
flux_limit = 10.**logflux_limit

def dummy(L,z):
    return 0

# Pozzetti model
dphi_dL  = {'0':dummy, 
            '1':pozzettiLF.model1Fast, 
            '3':pozzettiLF.model3Fast, 
            '9':dummy}[lf_model]
flux_key = {'0':'halo_lm', 
            '1':'logf_halpha_model1_ext', 
            '3':'logf_halpha_model3_ext', 
            '9':'euclid_nisp_H'}[lf_model]

if shuffled_fluxes:
    my_flux_key = 'sh_'+flux_key
else:
    my_flux_key = flux_key

# Conversion from Flux to Luminosity
def fluxToLuminosity (z, flux=None):

    if flux is None:
        flux=flux_limit
    return (flux * 4.0 * np.pi * (LCDMmodelLF.luminosity_distance(z).to(u.cm)) ** 2).value

# Integrated number of Galaxies per sq degree
def Pozzetti_dndz (z1, z2, minflux=None, maxflux=None):
    
    if maxflux is None:
        maxflux=1.e-13

    integrandL = lambda L,z: dphi_dL( L, z )
    integrandz = lambda z: 4.0 * np.pi / sqdegonthesky \
                 * LCDMmodelLF.differential_comoving_volume(z).value \
                 * quad(integrandL, fluxToLuminosity(z,minflux), fluxToLuminosity(z,maxflux), args=z)[0]

    return quad(integrandz, z1, z2)[0]


def build_fname(thisdir,tags,ext='.fits'):

    if thisdir[-1] != '/':
        thisdir += '/'

    basepath = path.join(outdir, thisdir)
    first = True

    fname = "catalogs_"
    for tag in tags:
        if not tag is None:
            if not first:
                if tag[0] != '_':
                    tag = '_' + tag
            first = False
            fname += tag
    fname += ext
    return path.join(basepath, fname)


def exclude_dir(fname):
    slashes=[pos for pos, char in enumerate(fname) if char == '/']
    return fname[slashes[-1]+1:]


# tags to be applied to file names

def lf_model_tag():
    if shuffled_fluxes:
        return 'm{}sh'.format(lf_model)
    else:
        return 'm{}'.format(lf_model)

def sm_tag():
    return 'smL{}'.format(smoothing_length)

def random_sm_tag():
    if smooth_dndz_in_random:
        return sm_tag()
    else:
        return 'nosm'

if RSD:
    redshift_key = 'observed_redshift_gal'
else:
    redshift_key = 'true_redshift_gal'

def rsd_tag():
    if RSD:
        return 'obsz'
    else:
        return 'truez'


def cm_tag():
    return 'cM{}'.format(cmrelation)
def alpha_tag():
    return '{}x'.format(alpha)
def ngrid_tag():
    return 'gr{}'.format(ngrid)
def Lbox_tag():
    return 'box{}'.format(Lbox)

# file names

def master_fname():
    return build_fname('RawCatalogs',[query,footprint_tag], ext='.zarr')

def indices_fname():
    return build_fname('RawCatalogs',[query,footprint_tag,'indices'], ext='.zarr')

def boxcat_fname():
    return build_fname('GalaxyCatalogs',['box',      otherlab,lf_model_tag(),cm_tag()])

def flagcat_fname():
    return build_fname('GalaxyCatalogs',['flagship', query,footprint_tag,lf_model_tag()])

def hodcat_fname():
    return build_fname('GalaxyCatalogs',['hodcat',   query,footprint_tag,otherlab,lf_model_tag(),cm_tag()])

def pincat_fname(r):
    return build_fname('GalaxyCatalogs',['pinocchio',r,  footprint_tag,otherlab,lf_model_tag(),cm_tag()])

def SDHOD_fname():
    return build_fname('SDHOD',['SDHOD',query,lf_model_tag(),sm_tag(),cm_tag()])

def galcat_fname(tag=cat_type):

    if tag == 'box':
        return boxcat_fname()
    elif tag == 'flagship':
        return flagcat_fname()
    elif tag == 'sdhod':
        return hodcat_fname()
    elif tag == 'pinocchio':
        return pincat_fname()
    else:
        print("ERROR in input.galcat_fname: unrecognised cat_type")
        raise ValueError()

def galcat_kernel(tag=cat_type):
    return exclude_dir(galcat_fname(tag))[:-5]

def selection_data_fname(sel_tag=selection_data_tag):
    if not sel_tag is None:
        return build_fname('Selections',['data',sel_tag,galcat_kernel()])
    else:
        return ''
def selection_random_fname(sel_tag=selection_random_tag):
    if not sel_tag is None:
        return build_fname('Selections',['random',sel_tag,galcat_kernel()])
    else:
        return ''

def selection_tag():
    if (selection_data_tag is None) & (selection_random_tag is None):
        return None
    else:
        return 'ds{}_rs{}'.format(selection_data_tag,selection_random_tag)

def dndz_fname(sel_tag=selection_data_tag):
    return build_fname('NumberCounts',['dndz',sel_tag,galcat_kernel(),sm_tag(),rsd_tag()])
def random_fname():
    return build_fname('RandomCatalogs',['random',galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag()])

def numbercounts_fname( z1, z2, sel_tag=selection_data_tag): 
    return build_fname('NumberCounts',['numbercounts',sel_tag,galcat_kernel(),rsd_tag(),'z{}-{}'.format(z1,z2)])
def cls_fname( z1, z2):
    return build_fname('Cls',['Cls',selection_tag(),galcat_kernel(),sm_tag(),rsd_tag(),'z{}-{}'.format(z1,z2)])

def LE3_data_fname( z1, z2):
    return build_fname('Catalogs4LE3',['data',selection_data_tag,galcat_kernel(),rsd_tag(),'z{}-{}'.format(z1,z2)])
def LE3_random_fname( z1, z2): 
    return build_fname('Catalogs4LE3',['random',selection_random_tag,galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag(),'z{}-{}'.format(z1,z2)])

def delta_data_fname( z1, z2):
    return build_fname('Cls',['delta_data',selection_data_tag,galcat_kernel(),rsd_tag(),'z{}-{}'.format(z1,z2)])
def delta_random_fname( z1, z2): 
    return build_fname('Cls',['delta_random',selection_random_tag,galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag(),'z{}-{}'.format(z1,z2)])

def params_fname( z1, z2): 
    return build_fname('PkParams',['params',selection_tag(),galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag(),ngrid_tag(),Lbox_tag(),'z{}-{}'.format(z1,z2)],ext='.ini')
def pk_fname( z1, z2): 
    return build_fname('Pks',['pk',selection_tag(),galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag(),ngrid_tag(),Lbox_tag(),'z{}-{}'.format(z1,z2)],ext='')
def script_fname( num): 
    return build_fname('PkScripts',['script',selection_tag(),galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag(),ngrid_tag(),Lbox_tag(),'{}'.format(num)],ext='.sh')

def plot_dndz_fname(sel_tag=selection_data_tag):
    return build_fname('Plots',['dndz',sel_tag,galcat_kernel(),sm_tag(),rsd_tag()],ext='.png')
def plot_numbercounts_fname( z1, z2, sel_tag=selection_data_tag):
    return build_fname('Plots',['numbercounts',sel_tag,galcat_kernel(),rsd_tag(),'z{}-{}'.format(z1,z2)],ext='.png')
def plot_cls_fname( z1, z2): 
    return build_fname('Plots',['cls',selection_tag(),galcat_kernel(),sm_tag(),rsd_tag(),'z{}-{}'.format(z1,z2)],ext='.png')
def plot_pk_fname( z1, z2): 
    return build_fname('Plots',['pk',selection_tag(),galcat_kernel(),random_sm_tag(),rsd_tag(),alpha_tag(),ngrid_tag(),Lbox_tag(),'z{}-{}'.format(z1,z2)],ext='.png')

