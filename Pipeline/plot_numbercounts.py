################################################################################
### Author: Tiago Castro, Pierluigi Monaco                                   ###
###                                                                          ###
################################################################################
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
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

print("# Running plot_numbercounts.py with {}".format(sys.argv[1]))

if not input.PLOT:
    print("Plots are not required, exiting...")
    sys.exit(0)

myrun=None
if len(sys.argv)>=3:
    if sys.argv[2].isdigit and input.cat_type is 'pinocchio':
        myrun=int(sys.argv[2])
        print("# I will process run number {}".format(myrun))
    else:
        print("# WARNING: unrecognised command-line option {}".format(sys.argv[2]))
else:
    if input.cat_type is 'pinocchio':
        myrun=input.pinocchio_first_run
        print("# I will process run number {}".format(myrun))

# survey footprint in equatorial coordinate
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
sky_coverage=input.sqdegonthesky * sky_fraction
print("This survey covers {} sq deg".format(sky_coverage))
del footprint

nzbins=len(input.finalCatZShell)
zbins=np.asarray(input.finalCatZShell)

for i in range(nzbins):

    z1=zbins[i,0]
    z2=zbins[i,1]

    fname = input.numbercounts_fname(z1,z2,run=myrun)
    print("# Reading number counts from file {}".format(fname))
    LF = fits.getdata(fname)

    fig=plt.figure(figsize=(8, 8))
    plt.suptitle(input.exclude_dir(fname))

    gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1], hspace=0)

    panel1 = plt.subplot(gs[0])
    panel2 = plt.subplot(gs[1])

    panel1.set_xlim([1.5e-16, 1.e-14])
    panel2.set_xlim([1.5e-16, 1.e-14])
    panel1.set_ylim([1., 1.5e4])
    panel1.set_xscale('log')
    panel1.set_yscale('log')
    panel2.set_xscale('log')
    panel2.set_yscale('linear')
    panel2.set_ylim([0.7,1.3])
    panel2.set_xlabel(r'flux (erg cm$^{-2}$ s$^{-1}$)')
    panel1.set_ylabel(r'number counts (sq deg$^{-1}$ redshift$^{-1}$ mag$^{-1}$)')
    panel2.set_ylabel(r'residuals vs model')
    plt.tight_layout(pad=1.5)

    Nmodel=[]
    for f1,f2 in zip(LF['f_lower'],LF['f_upper']):
        DeltaF=2.5*(np.log10(f2)-np.log10(f1))
        Ng=input.Pozzetti_dndz(z1,z2,f1,f2)
        Nmodel.append(Ng/(z2-z1)/DeltaF)
    Nmodel=np.asarray(Nmodel)
    pos=Nmodel>0

    panel1.plot(LF['f_center'],Nmodel,label='model',c='k')
    panel2.plot(LF['f_center'],np.ones_like(Nmodel),c='k')

    panel1.plot(LF['f_center'],LF['LF'],label='data',c='red')
    panel2.plot(LF['f_center'][pos],LF['LF'][pos]/Nmodel[pos],c='red')
    panel1.plot(LF['f_center'],LF['LF_cen'],'--',label='data, centrals',c='red')
    panel1.plot(LF['f_center'],LF['LF_sat'],':',label='data, satellites',c='red')

    panel1.legend()

    plt.savefig(input.plot_numbercounts_fname(z1,z2,myrun))
    print("## written image in file {}".format(input.plot_numbercounts_fname(z1,z2)))



if input.SHOW:
    plt.show()



