import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

from scipy.interpolate import RegularGridInterpolator


lut = fits.getdata('lookup_table_15by15.fits')

nred=14
nHal=14
nexp=4
nbkg=13

ndata = nred * nHal * nexp * nbkg
if ndata != len(lut):
    print("Error: table length is %d while I expected %d")
    sys.exit(0)

print("Table dimension: {} * {} * {} * {} = {}".format(nred,nHal,nexp,nbkg,ndata))

xred = lut['z'].reshape(nbkg,nexp,nHal,nred)
xHal = np.log10(lut['flux_Ha']).reshape(nbkg,nexp,nHal,nred)
xexp = lut['exp_time'].reshape(nbkg,nexp,nHal,nred)
xbkg = lut['sky_bg'].reshape(nbkg,nexp,nHal,nred)

prob = lut['galaxy_frac'].reshape(nbkg,nexp,nHal,nred)

lim=np.log10(2.e-16)

plt.figure()

iexp=3
ibkg=0
plt.title('background={}, exposure time={}'.format(xbkg[ibkg,0,0,0],xexp[0,iexp,0,0]))
for ired in np.arange(nred):
    plt.plot(xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired],label='z={:4.2f}'.format(xred[0,0,0,ired]))
plt.legend()
plt.plot([lim,lim],[-0.2,1.2],':',c='k')
plt.ylim([-0.05,1.1])
plt.xlabel('Log Halpha flux')
plt.ylabel('detection probability')

plt.savefig('prob_redshift.png')

plt.figure()

ired=5
ibkg=0
plt.title('background={}, z={:4.2f}'.format(xbkg[ibkg,0,0,0],xred[0,0,0,ired]))
for iexp in np.arange(nexp):
    plt.plot(xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired],label='exp={}'.format(xexp[0,iexp,0,0]))
plt.plot([lim,lim],[-0.2,1.2],':',c='k')
plt.legend()
plt.ylim([-0.05,1.1])
plt.xlabel('Log Halpha flux')
plt.ylabel('detection probability')

plt.savefig('prob_exptime.png')

plt.figure()

ired=5
iexp=3
plt.title('exposure time={}, z={:4.2f}'.format(xexp[0,iexp,0,0],xred[0,0,0,ired]))
for ibkg in np.arange(nbkg):
    plt.plot(xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired],label='backgr={}'.format(xbkg[ibkg,0,0,0]))
plt.plot([lim,lim],[-0.2,1.2],':',c='k')
plt.legend()
plt.ylim([-0.05,1.1])
plt.xlabel('Log Halpha flux')
plt.ylabel('detection probability')

plt.savefig('prob_background.png')

iexp=3
ibkg=0
plt.figure()
plt.title('background={}, exposure time={}'.format(xbkg[ibkg,0,0,0],xexp[0,iexp,0,0]))
plt.xlabel('redshift')
plt.ylabel('completeness at flux limit')
plt.ylim([-0.05,1.1])

for iflux in np.arange(2,7):
    flux = lim + np.log10(iflux/2.)
    compl=[]
    for ired in np.arange(nred):
        compl.append(np.interp(flux,xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired]))
    plt.plot(xred[0,0,0],np.asarray(compl),label='limit={}e-16'.format(iflux))
plt.legend()

plt.savefig('compl_flux.png')

ibkg=0
plt.figure()
plt.title('limit=2e-16, background={}'.format(xbkg[ibkg,0,0,0]))
plt.xlabel('redshift')
plt.ylabel('completeness at flux limit')
plt.ylim([-0.05,1.1])

for iexp in np.arange(nexp):
    flux = lim
    compl=[]
    for ired in np.arange(nred):
        compl.append(np.interp(flux,xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired]))
    plt.plot(xred[0,0,0],np.asarray(compl),label='exp={}'.format(xexp[0,iexp,0,0]))
plt.legend()

plt.savefig('compl_exptime.png')


iexp=3
plt.figure()
plt.title('limit=2e-16, exposure time={}'.format(xexp[0,iexp,0,0]))
plt.xlabel('redshift')
plt.ylabel('completeness at flux limit')
plt.ylim([-0.05,1.1])

for ibkg in np.arange(nbkg):
    flux = lim
    compl=[]
    for ired in np.arange(nred):
        compl.append(np.interp(flux,xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired]))
    plt.plot(xred[0,0,0],np.asarray(compl),label='backgr={}'.format(xbkg[ibkg,0,0,0]))
plt.legend()

plt.savefig('compl_background.png')

plt.show()


my_interp = RegularGridInterpolator( ( xbkg[:,0,0,0],xexp[0,:,0,0],xHal[0,0,:,0],xred[0,0,0,:] ),
                                  prob, method='linear' )


ibkg=0
iexp=3
plt.figure()
plt.title('limit=2e-16, background={}'.format(xbkg[ibkg,0,0,0]))
plt.xlabel('redshift')
plt.ylabel('completeness at flux limit')
plt.ylim([-0.05,1.1])

compl=[]
for ired in np.arange(nred):
    compl.append(np.interp(lim,xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired]))

plt.plot(xred[0,0,0],np.asarray(compl))

iexp=2
compl=[]
for ired in np.arange(nred):
    compl.append(np.interp(lim,xHal[ibkg,iexp,:,ired],prob[ibkg,iexp,:,ired]))

plt.plot(xred[0,0,0],np.asarray(compl))

xx=np.linspace(0.9,1.8,num=44)
ones=np.ones_like(xx)
for ee in np.linspace(xexp[0,2,0,0],xexp[0,3,0,0],num=7):
    points = np.array([ xbkg[ibkg,0,0,0]*ones, ee*ones, flux*ones, xx ]).transpose()
    plt.plot(xx,my_interp(points))

plt.legend()

