################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.io import fits
from fitsio import FITS, FITSHDR
import healpy as hp
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


print("Sorry, this is still work in progress!")
sys.exit(0)


def delta_map (ra, dec, footprint, radians=False):
# this function builds a density contrast from a list of galaxies
# pix: healpix sky pixel where the galaxy falls
# npix: number of pixels
# footprint: survey footprint

    print("## building map...")

    if radians:
        conv=1.
    else:
        conv=np.pi/180.

    if footprint is None:
        footprint = np.ones(hp.nside2npix(input.footprint_res),dtype=bool)

    pix = hp.ang2pix(input.footprint_res, np.pi/2. - dec*conv, ra*conv)
    delta=np.zeros(hp.nside2npix(input.footprint_res),dtype=np.float)
    for this in pix:
        delta[this]+=1.0
    #delta  = np.histogram(pix, bins=np.linspace(0, npix, npix+1))[0].astype(np.float)
    av = np.mean(delta[footprint])
    if av == 0:
        print("Error in delta_map: zero average")
        return None

    delta = (delta-av)/av
    delta[~footprint] = hp.UNSEEN

    return delta

def write_delta(d,fname):
    towrite=np.empty(delta.size, dtype=[('DELTA', np.float)])
    towrite['DELTA']=d
    fits.writeto(fname,towrite,overwrite=True)
    del towrite

print("# Running maps_and_cls.py with {}".format(sys.argv[1]))


# reads the survey footprint in equatorial coordinates
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
print("# the footprint covers a sky fraction of %f"%sky_fraction)

conv=np.pi/180.

ell=np.arange(3*footprint_res)

for zmin,zmax in input.finalCatZShell:

    fname = input.LE3_random_fname(zmin, zmax)
    print("## reading file %s"%fname)
    data = fits.getdata(fname)
    print("## computing density map on %d galaxies"%len(data))
    delta = delta_map(data['RA'], data['DEC'],footprint)

    if delta is not None:
        write_delta(delta,input.delta_random_fname(zmin,zmax))
        print("## computing Cls")
        cl_rand = hp.anafast(delta)/sky_fraction
    else:
        cl_rand = np.zeros_like(3*footprint_res)

    fname = input.LE3_data_fname(zmin,zmax, run=input.pinocchio_first_run)
    print("## reading file %s"%fname)
    data = fits.getdata(fname)
    print("## computing density map on %d galaxies"%len(data))
    delta = delta_map(data['RA'], data['DEC'],footprint)

    if delta is not None:
        write_delta(delta,input.delta_data_fname(zmin,zmax, run=input.pinocchio_first_run))
        print("## computing Cls")
        cl_data = hp.anafast(delta)/sky_fraction
    else:
        cl_data = np.zeros_like(3*footprint_res)

    columns = [['ell', ell,'f8'],
               ['cl_data', cl_data,'f8'],
               ['cl_rand', cl_rand,'f8']]

    header_keywords = {}

    types = []
    # keep just wanted columns ...bad but..
    for c in columns :
        types.append((c[0], c[2]))

    hdr = FITSHDR()
    print("+ Add keywords")
    for k in header_keywords :
        hdr[k] = header_keywords[k]

    keep_table = {}
    for c in columns :
        # keep column with  requested position in the input file
        print(str("+  -Col '%s'" % (c[0])))
        keep_table[c[0]] = c[1] #.astype(np.float64)

    fullsize = len(keep_table)*len(c[1])* 8 / 1024 / 1024.
    print(str("+   ~%.2f MB in memory" % fullsize))

    fname = input.cls_fname(zmin,zmax, run=input.pinocchio_first_run)
    print("## Writing FITS: %s" % fname)
    fitsfile = FITS(fname,'rw')
    fitsfile.write_table(data=keep_table, header=hdr)
    fitsfile.close()

print("# DONE!")
