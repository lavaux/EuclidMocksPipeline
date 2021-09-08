################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from astropy.io import fits
import numpy as np
import sdhod
import NFW
import utils
from colossus.halo import concentration
from scipy.stats import poisson
import sys

if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running createSDHODCatalog.py with {}".format(sys.argv[1]))

footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
del footprint

# this is needed to calibrate the 1-halo term in redshift space
SCALE_VELOCITIES = 0.7

# seed for random numbers
np.random.seed(seed=input.SEED_hod)

print("# Reading the HOD table {}...".format(input.SDHOD_fname()))
hodtable=fits.getdata(input.SDHOD_fname())

print('# Reading the halo catalog from {}...'.format(input.master_fname()))
rawcat = fits.getdata(input.master_fname())

c_kind     = rawcat['kind']
c_z        = rawcat['true_redshift_gal']
c_logm     = rawcat['halo_lm']
c_ra_gal   = rawcat['ra_gal']
c_dec_gal  = rawcat['dec_gal']
c_x_gal    = rawcat['x_gal']
c_y_gal    = rawcat['y_gal']
c_z_gal    = rawcat['z_gal']
c_vrad_gal = rawcat['vrad_gal']
c_halo_id  = rawcat['halo_id']

del rawcat

print('# Filtering the catalog for main halos...')
filt = c_kind==0
del c_kind

c_z       = c_z       [filt]
c_logm    = c_logm    [filt]
c_ra_gal  = c_ra_gal  [filt]
c_dec_gal = c_dec_gal [filt]
c_x_gal   = c_x_gal   [filt]
c_y_gal   = c_y_gal   [filt]
c_z_gal   = c_z_gal   [filt]
c_vrad_gal= c_vrad_gal[filt]
c_halo_id = c_halo_id [filt]

Nhalos  = c_z.size

print('# Found {} main halos'.format(Nhalos))

print("# Random sampling the Centrals...")
NCen        = sdhod.Ncen(hodtable, c_logm, c_z)
# The halo hosts an Halpha emitter if a random variable is less than NCen
hostCentral = np.random.rand(Nhalos) <= NCen
totCentrals = hostCentral.sum()

print("# There will be {} central galaxies".format(totCentrals))

print("# Random sampling the Satellites...")
# Getting the number of Satellites according to a poisson distribution
Nsat    = poisson.rvs( sdhod.Nsat(hodtable,c_logm, c_z) )
hostSat = Nsat>0
totSat  = Nsat.sum()

print("# There will be {} satellite galaxies".format(totSat))

# Total number of Galaxies Sat+Cen
Ngal = int(totSat + totCentrals)

Mass= 10.**c_logm

print('# Creating the galaxy variables...')
ra_gal     = np.empty(Ngal, dtype=np.float)
dec_gal    = np.empty(Ngal, dtype=np.float)
xgal       = np.empty(Ngal, dtype=np.float)
ygal       = np.empty(Ngal, dtype=np.float)
zgal       = np.empty(Ngal, dtype=np.float)
rgal       = np.empty(Ngal, dtype=np.float)
vlosgal    = np.empty(Ngal, dtype=np.float)
true_zgal  = np.empty(Ngal, dtype=np.float)
obs_zgal   = np.empty(Ngal, dtype=np.float)
halo_m     = np.empty(Ngal, dtype=np.float)
log10f     = np.empty(Ngal, dtype=np.float)
kind       = np.empty(Ngal, dtype=np.int)
haloid     = np.empty(Ngal, dtype=np.int64)

cid        = np.arange(Ngal)

print('# Working on the {} centrals...'.format(totCentrals))

ra_gal[:totCentrals]    = c_ra_gal[hostCentral]
dec_gal[:totCentrals]   = c_dec_gal[hostCentral]
xgal[:totCentrals]      = c_x_gal[hostCentral]
ygal[:totCentrals]      = c_y_gal[hostCentral]
zgal[:totCentrals]      = c_z_gal[hostCentral]
rgal[:totCentrals]      = np.sqrt(xgal[:totCentrals]**2 + ygal[:totCentrals]**2 + zgal[:totCentrals]**2)
vlosgal[:totCentrals]   = c_vrad_gal[hostCentral]
true_zgal[:totCentrals] = c_z[hostCentral]
halo_m[:totCentrals]    = Mass[hostCentral]
kind[:totCentrals]      = 0
haloid[:totCentrals]    = c_halo_id[hostCentral]

print("# Assigning fluxes to the centrals...")

log10f[:totCentrals]    = sdhod.lfcen(hodtable,c_logm[hostCentral], c_z[hostCentral])

print('# Working on the {} Satellites'.format(totSat))

print('# Calculating concentrations according to {}...'.format(input.cmrelation))

kind[totCentrals:] = 1

MDelta         = Mass[hostSat]
z2             = c_z[hostSat]
concentrations = np.array( [concentration.concentration(MM, '200c', z=zz, model = input.cmrelation ) for MM, zz in zip(MDelta, z2)] )
RDelta         = ( 3.0*MDelta/4.0/np.pi/200.0/input.cosmo.rho_c(z2) )**(1.0/3.0)
RDelta        /= 1e3 # To Mpc/h

print('# Distributing satellite galaxies in halos...')

NwithSat = np.count_nonzero(hostSat)
Nsat=Nsat[hostSat]

myhalo = np.zeros(totSat,dtype=np.int)
cc=0
for i in range(NwithSat):
    myhalo[cc:cc+Nsat[i]]=i
    cc+=Nsat[i]

conv=np.pi/180.

randr = np.array( [ NFW.getr(c,u)[0] for c,u in zip(concentrations[myhalo],np.random.rand(totSat))] )
randt, randp = utils.randomSpherePoint(totSat)
randv = np.array( [ utils.circularVelocity(c,r) for c,r in zip(concentrations[myhalo],randr) ] ) * np.sqrt(MDelta[myhalo]/RDelta[myhalo])

print("# random numbers done")
print("# updating spatial positions...")

dist = np.sqrt(c_x_gal[hostSat]**2+c_y_gal[hostSat]**2+c_z_gal[hostSat]**2)
sint = np.sin(randt)
cosd = np.cos(c_dec_gal[hostSat][myhalo]*conv)
xgal[totCentrals:] = dist[myhalo] * cosd * np.cos(c_ra_gal[hostSat][myhalo]*conv) + RDelta[myhalo] * randr * sint * np.cos(randp)
ygal[totCentrals:] = dist[myhalo] * cosd * np.sin(c_ra_gal[hostSat][myhalo]*conv) + RDelta[myhalo] * randr * sint * np.sin(randp)
zgal[totCentrals:] = dist[myhalo] * np.sin(c_dec_gal[hostSat][myhalo]*conv)       + RDelta[myhalo] * randr * np.cos(randt)

rgal[totCentrals:] = np.sqrt(xgal[totCentrals:]**2 + ygal[totCentrals:]**2 + zgal[totCentrals:]**2)

print("# updating angular positions...")
dec_gal[totCentrals:] = np.arcsin(zgal[totCentrals:]/rgal[totCentrals:])
ra_gal[totCentrals:] = np.arctan2(ygal[totCentrals:],xgal[totCentrals:])

v = SCALE_VELOCITIES * np.array( [ randv * np.random.normal(size=totSat), randv * np.random.normal(size=totSat), randv * np.random.normal(size=totSat) ] )

vlosgal[totCentrals:] = c_vrad_gal[hostSat][myhalo] + (v[0,:]*xgal[totCentrals:]+v[1,:]*ygal[totCentrals:]+v[2,:]*zgal[totCentrals:])/rgal[totCentrals:]

haloid[totCentrals:] = c_halo_id[hostSat][myhalo]
halo_m[totCentrals:] = Mass[hostSat][myhalo]
true_zgal[totCentrals:] = c_z[hostSat][myhalo]
log10f[totCentrals:] = sdhod.lfsat(hodtable,c_logm[hostSat][myhalo],true_zgal[totCentrals:])

obs_zgal = true_zgal + vlosgal / input.SPEEDOFLIGHT * (1.0 + true_zgal)

ra_gal *= 180./np.pi
ra_gal[ra_gal<0] += 360.
dec_gal *= 180./np.pi

print("# Shuffling galaxy fluxes...")
shuffled_log10f = np.copy(log10f)

zbins = np.linspace(footprint_zrange[0],footprint_zrange[1],round( (footprint_zrange[1]-footprint_zrange[0] )/input.deltazbin + 1 ) )
zindex = np.digitize( true_zgal, zbins )
for iz in np.arange(zbins.size-1):
    ff = zindex==iz
    extract = shuffled_log10f[ff]
    np.random.shuffle(extract)
    shuffled_log10f[ff] = extract

fname=input.hodcat_fname()
print("# Saving the catalog to file {}".format(fname))
catalog = np.empty(Ngal,
              dtype=[('x_gal', np.float), ('y_gal', np.float), ('z_gal', np.float),
                     ('ra_gal', np.float), ('dec_gal', np.float), ('kind', np.int), 
                     ('true_redshift_gal', np.float), ('observed_redshift_gal', np.float), 
                     ('halo_lm', np.float), ('id', np.int), ('halo_id', np.int), 
                     (input.flux_key, np.float), ('sh_'+input.flux_key, np.float)])

catalog['x_gal']                  = xgal
catalog['y_gal']                  = ygal
catalog['z_gal']                  = zgal
catalog['true_redshift_gal']      = true_zgal
catalog['observed_redshift_gal']  = obs_zgal
catalog['ra_gal']                 = np.arctan2(ygal, xgal) * 180.0/np.pi
catalog['dec_gal']                = 90.0 - np.arccos( zgal/(np.sqrt( xgal**2 + ygal**2 + zgal**2 )) ) * 180.0/np.pi
catalog['halo_lm']                = np.log10(halo_m)
catalog['kind']                   = kind
catalog[input.flux_key]           = log10f
catalog['sh_'+input.flux_key]     = shuffled_log10f
catalog['id']                     = cid
catalog['halo_id']                = haloid

fits.writeto(fname, catalog, overwrite=True)
print("# !!DONE!!")
