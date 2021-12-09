import numpy as np
from . import sdhod
from astropy.io import fits
from scipy import interpolate as interp
from colossus.halo import concentration


# match = fits.getdata('../Products/Pinocchio/match_FS_pin_masses.fits')
# FSmass = interp.RectBivariateSpline(match['pinmasses'],match['redshift'][0],match['FSmasses'][0])

cin = np.array([-0.00018763, -0.04794911, 0.12291895])
csl = np.array([0.01885588, -0.0737355, 0.8910409])
ccu = np.array([-0.00316799, 0.01848572, -0.00340885])
CMFc = np.array([-0.05002674, -0.10366184, -0.9558882])

patch = np.array([0.561958, -2.96723156, 5.62666982, -4.52795846, 2.32903134])


def logFSmasses(M, z):

    inte = cin[0] * z ** 2 + cin[1] * z + cin[2]
    slop = csl[0] * z ** 2 + csl[1] * z + csl[2]
    curv = ccu[0] * z ** 2 + ccu[1] * z + ccu[2]

    x = np.log10(M) - 12.0

    return 12.0 + inte + slop * x + curv * x ** 2


def correct_number_density(mass_shift, z):
    # this gives a correction to galaxy number density to compensate for mass shift
    return 10.0 ** (-mass_shift * np.poly1d(CMFc)(z)) * (
        patch[0] * z ** 4
        + patch[1] * z ** 3
        + patch[2] * z ** 2
        + patch[3] * z
        + patch[4]
    )
