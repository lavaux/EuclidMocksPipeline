import numpy as np
import astropy.units as u

# Model 1
alpha_1       = -1.35
L_star0       = 10 ** 41.5 * u.erg/u.s
L_star0_u     = L_star0.value
phi_star0_1   = 10 ** -2.8 / (u.Mpc)**3
phi_star0_1_u = phi_star0_1.value
delta         = 2.0
epsilon       = 1.0
zbreak        = 1.3

# Model 3
alpha_3          = -1.587
beta             =  1.615
Delta            =  2.288
phi_star0_3      = 10 ** (-2.920) / (u.Mpc)**3
phi_star0_3_u    = phi_star0_3.value
log10_L_star20   = 42.557
log10_L_star05   = 41.733
log10_L_starIn   = 42.956

# Model 1
def model1(L, z):

    L_starz = L_star0 * (1.0 + z) ** delta

    if z < zbreak:

        phi_starz = phi_star0_1 * (1.0 + z)**epsilon

    else:

        phi_starz = phi_star0_1 * (1.0 + z) ** (-epsilon) * (1.0 + zbreak) ** (2 * epsilon)


    return  phi_starz * (L/L_starz) ** alpha_1 * np.exp(-L/L_starz) / L_starz

# Model 1 wout Units (FAST)
def model1Fast(L, z):

    L_starz = L_star0_u * (1.0 + z) ** delta

    if z < zbreak:

        phi_starz = phi_star0_1_u * (1.0 + z)**epsilon

    else:

        phi_starz = phi_star0_1_u * (1.0 + z) ** (-epsilon) * (1.0 + zbreak) ** (2 * epsilon)


    return  phi_starz * (L/L_starz) ** alpha_1 * np.exp(-L/L_starz) / L_starz

# Model 3
def model3(L, z):

    L_starz =  10 ** ( log10_L_starIn + ((1.5/(1.0+z)) ** beta) * (log10_L_star05-log10_L_starIn) )  * u.erg / u.s
    
    return phi_star0_3/L_starz * (L/L_starz) ** alpha_3 * 1.0/( 1.0 + (np.exp(1.0) - 1.0) * (L/L_starz) ** Delta )

# Model 3 wout Unit (FAST)
def model3Fast(L, z):

    L_starz =  10 ** ( log10_L_starIn + ((1.5/(1.0+z)) ** beta) * (log10_L_star05-log10_L_starIn) ) 

    return phi_star0_3_u/L_starz * (L/L_starz) ** alpha_3 * 1.0/( 1.0 + (np.exp(1.0) - 1.0) * (L/L_starz) ** Delta )
