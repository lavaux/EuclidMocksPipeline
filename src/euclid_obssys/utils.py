import numpy as np
from astropy.constants import G
from .NFW import normalizedNFWMass

c = 299792.0 # Km/s
G = G.to("Mpc Msun^-1 km^2 s^-2").value # G Newton in Mpc/Msun x (km/s)^2

def randomSpherePoint (n=1, theta_max=np.pi, phi_max=2.0*np.pi):
    '''
    Returns random distributes theta (0,pi) and phi (0, 2pi)
    '''
    u = 1.0 - np.random.rand(n) * theta_max/np.pi
    v = np.random.rand(n)

    return np.arccos(2*u-1), phi_max *v

def getDeltaZ (x, y, z, vx, vy, vz):
    '''
    Returns the redshift due to Pec. Vel.

    x,   y,  z: Cartesian coordinates [Not Spec. Units]
    vx, vy, vz: Cartesian velocities  [It must be in km/s]
    '''
    return (x*vx + y*vy + z*vz)/np.sqrt( x*x + y*y + z*z )/c

def circularVelocity(cx, r):
    '''
    Returns the circular velocity at r=R/Rx divided by np.sqrt(MDelta/RDelta)

    cx is the concentration of the Delta_x SO halo
    '''
    return np.sqrt(G * normalizedNFWMass(cx, r)/r)

def randomPerpendicularDirection (x, y, z):
    '''
    Returns random vector perpendicular to x, y, z
    '''

    norm = np.sqrt(x*x + y*y +z*z)

    vn = np.array([x/norm, y/norm, z/norm])

    # Creating a vector perpendicular to [x,y,z]
    perp = np.cross(vn, [1.0, 0.0, 0.0])
    perp = perp/np.sqrt( np.dot(perp, perp) )

    # Applying Rodrigues formula for random phi
    phi = np.random.rand() * 2 * np.pi
    vrot = vn * np.cos(phi) + np.perp(perp, vn) * np.sin(phi) + perp * np.dot(perp, vn) * (1-np.cos(phi))
