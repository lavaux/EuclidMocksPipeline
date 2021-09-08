################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
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


def divide(a,b):
    filtro=b>0
    res=np.zeros_like(a)
    res[filtro]=a[filtro]/b[filtro]
    return res

print("# Running visualizeHOD.py with {}".format(sys.argv[1]))

c=np.array(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

sdhod = fits.getdata(input.SDHOD_fname())

izbins=np.arange(12)*10
Mx=10.**(0.5*(sdhod['M_bins'][0,1:]+sdhod['M_bins'][0,:-1]))
zx=0.5*(sdhod['z_bins'][0,1:]+sdhod['z_bins'][0,:-1])

plt.figure()
plt.title('SD-HOD, centrals')
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.8e11,1.e15])
plt.ylim([1.e-5,100])

for ic in np.arange(10):
    z=zx[izbins[ic]]
    y=divide(sdhod['n_cen_gaus'][0, 0, izbins[ic]],sdhod['n_hal_gaus'][0, izbins[ic]])
    ff=y>0
    plt.plot(Mx[ff],y[ff]*z**8,label="%3.1f"%z,c=c[ic],linewidth=3)

    y=divide(sdhod['n_cen'][0, 0, izbins[ic]],sdhod['n_halos'][0, izbins[ic]])
    ff=y>0
    plt.plot(Mx[ff],y[ff]*z**8,c=c[ic])

plt.xlabel('halo mass (Msun/h)')
plt.ylabel(r'$N_{\rm cen} (M_h|z)$')
plt.legend()
plt.show()
plt.savefig(input.outdir + 'Plots/SDHOD_centrals.png')

plt.figure()
plt.title('SD-HOD, satellites')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1.e-5,1e3])

for ic in np.arange(10):
    z=zx[izbins[ic]]
    y=divide(sdhod['n_sat_gaus'][0, 0, izbins[ic]],sdhod['n_hal_gaus'][0, izbins[ic]])
    ff=y>0
    plt.plot(Mx[ff],y[ff]*z**8,label="%3.1f"%z,c=c[ic],linewidth=3)

    y=divide(sdhod['n_sat'][0, 0, izbins[ic]],sdhod['n_halos'][0, izbins[ic]])
    ff=y>0
    plt.plot(Mx[ff],y[ff]*z**8,c=c[ic])

plt.xlabel('halo mass (Msun/h)')
plt.ylabel(r'$N_{\rm cen} (M_h|z)$')
plt.legend()
plt.show()
plt.savefig(input.outdir + 'Plots/SDHOD_centrals.png')

