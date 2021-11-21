################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.io import fits
import healpy as hp
from scipy.ndimage import gaussian_filter1d
import sys
from os import path
import filenames


if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running dndz.py with {}".format(sys.argv[1]))

# loads the survey footprint in equatorial coordinates
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()

# redshift binning for computing dn/dz
print("# setting binnings...")

ztab = np.linspace(footprint_zrange[0], footprint_zrange[1], int( np.round( (footprint_zrange[1]-footprint_zrange[0])/input.deltazbin) + 1) )

if np.abs((ztab[1]-ztab[0])/input.deltazbin-1) > 1e-5:
    print("!The chosen redshift binning and redshift range are not compatible!")
    sys.exit(-1)

# comoving distances for bins
bin_edges = np.asarray([ input.LCDMmodel.comoving_distance(z).to(input.l_unit).value for z in ztab ])
# bin volume and center
bin_volume = 4. * np.pi/3. * (bin_edges[1:]**3 - bin_edges[:-1]**3) * sky_fraction
bin_center = 0.5 * (bin_edges[:-1] + bin_edges[1:])


if (input.cat_type is 'pinocchio') & (input.pinocchio_last_run is not None):

    Ngal=None

    nruns=input.pinocchio_last_run-input.pinocchio_first_run +1
    for myrun in np.arange(input.pinocchio_first_run,input.pinocchio_last_run+1):

        myfname=filenames.galcat(input,myrun)
        print("# loading catalog {}...".format(myfname))
        
        if not path.exists(myfname):
            print("ERROR: galaxy catalog {} does not exist".format(myfname))
            sys.exit(-1)

        cat = fits.getdata(myfname)

        if input.selection_data_tag is not None:
            myfname=filenames.selection_data(input,myrun)
            print("# loading selection {}...".format(myfname))
            mysel = fits.getdata(myfname)['SELECTION']
        else:
            mysel = np.ones(len(cat),dtype=bool)

        print("# Processing catalog...")

        # Histogram
        isCen    = cat['kind'][mysel]==0
        if Ngal is None:
            Ngal  = (np.histogram(cat[input.redshift_key][mysel], bins=ztab)[0]).astype(np.float64)
            Ncen  = (np.histogram(cat[input.redshift_key][mysel][isCen], bins=ztab)[0]).astype(np.float64)
        else:
            Ngal += (np.histogram(cat[input.redshift_key][mysel], bins=ztab)[0]).astype(np.float64)
            Ncen += (np.histogram(cat[input.redshift_key][mysel][isCen], bins=ztab)[0]).astype(np.float64)

    Ngal = (Ngal/float(nruns)).astype(float)
    Ncen = (Ncen/float(nruns)).astype(float)

else:

    myfname=filenames.galcat(input,input.pinocchio_first_run)
    print("# loading catalog {}...".format(myfname))

    if not path.exists(myfname):
        print("ERROR: galaxy catalog {} does not exist".format(myfname))
        sys.exit(-1)

    cat = fits.getdata(myfname)

    # selection
    if input.selection_data_tag is not None:
        myfname=filenames.selection_data(input,input.pinocchio_first_run)
        print("# loading selection {}...".format(myfname))
        mysel = fits.getdata(myfname)['SELECTION']
    else:
        mysel = np.ones(len(cat),dtype=bool)


    print("# Processing catalog...")

    # Histogram

    Ngal     = (np.histogram(cat[input.redshift_key][mysel], bins=ztab)[0]).astype(float)
    isCen    = cat['kind'][mysel]==0
    Ncen     = (np.histogram(cat[input.redshift_key][mysel][isCen], bins=ztab)[0]).astype(float)

## Writes on file
myfname = filenames.dndz(input)
print("# Writing results on file {}...".format(myfname))

dndz = np.empty(Ngal.size,
                         dtype=[('N_gal', float), ('N_gal_gaus', float), 
                                ('N_cen', float), ('N_cen_gaus', float),
                                ('z_center', float), ('z_lower', float), ('z_upper', float),
                                ('bin_center', float), ('bin_lower', float), 
                                ('bin_upper', float), ('bin_volume', float)])

dndz['N_gal'] = Ngal;
dndz['N_cen'] = Ncen;
dndz['bin_center'] = bin_center;
dndz['bin_lower'] = bin_edges[:-1]; dndz['bin_upper'] = bin_edges[1:];
dndz['bin_volume'] = bin_volume;
dndz['z_center'] = (ztab[1:] + ztab[:-1])/2; dndz['z_lower'] = ztab[:-1]; dndz['z_upper'] = ztab[1:];

dndz['N_gal_gaus'] = gaussian_filter1d(Ngal, input.smoothing_length);
dndz['N_cen_gaus'] = gaussian_filter1d(Ncen, input.smoothing_length);

fits.writeto(myfname, dndz, overwrite=True)

print("# done!")
