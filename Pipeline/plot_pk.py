################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from glob import glob
import matplotlib.gridspec as gridspec
import sys
import importlib
import filenames

if len(sys.argv)<2:
    print("Usage: python {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file {} not found".format(sys.argv[1]))
    print("Usage: python {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running plot_pk.py with {}".format(sys.argv[1]))

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

type1='pinocchio'
type2='pinocchio'
lab1='base'
lab2='extinction'

for z1,z2 in input.finalCatZShell:

    importlib.reload(input)
    input.cat_type=type1
    fname = glob(filenames.estimator_measure(input,'PK',z1,z2,myrun)+"/EUC_LE3_GCL_PK_*0Z_1D.fits")[0]
    print("Reading PK from file {}".format(fname))
    data1 = fits.getdata( fname )

    if type2 is not None:
        old_sdt = input.selection_data_tag
        old_srt = input.selection_random_tag
        input.selection_data_tag='lutns'
        input.selection_random_tag='lutns'
        #input.apply_dataselection_to_random=True
        fname = glob(filenames.estimator_measure(input,'PK',z1,z2,myrun)+"/EUC_LE3_GCL_PK_*0Z_1D.fits")[0]
        print("Reading PK from file {}".format(fname))
        data2 = fits.getdata( fname )

        input.selection_data_tag=old_sdt
        input.selection_random_tag=old_srt
        #input.apply_dataselection_to_random=False
    

    for order in [0,2,4]:

        title = r"z=[{},{}], PK{}".format(z1,z2,order)
        fig  = plt.figure(figsize=(6, 8))
        plt.suptitle(title)
        gs   = gridspec.GridSpec(2,1,height_ratios=[1.5,1], hspace=0)
        plot = fig.add_subplot(gs[0])
        res  = fig.add_subplot(gs[1], sharex=plot)
        plt.xscale('log')
        #plt.yscale('log')
        plot.set_yscale('log')
        string=r'P_{}'.format(order)
        plot.set_ylabel(r'$'+string+'(k)$ ($h^{-3}$ Mpc)')
        res.set_xlabel(r'$k$ ($h$/Mpc)')
        res.set_ylabel('data2/data1 -1')
        res.set_ylim([-0.15,0.15])

        plot.plot(data1['K'], data1['PK{}'.format(order)], label=lab1,c='r')
        plot.plot(data1['K'],-data1['PK{}'.format(order)], ':', c='r')

        if type2 is not None:

            plot.plot(data2['K'], data2['PK{}'.format(order)], label=lab2,c='b')
            plot.plot(data2['K'],-data2['PK{}'.format(order)], ':', c='b')

            res.plot(data1['K'], data2['PK{}'.format(order)]/data1['PK{}'.format(order)]-1., c='b')

        res.plot(data1['K'],np.zeros_like(data1['K']),c='r')
        plot.legend()

        fname=filenames.plot_estimator(input,'PK',order,z1,z2,myrun)
        fig.savefig(fname)
        print("## written image in file "+fname)



if input.SHOW:
    plt.show()

print("# DONE!")
