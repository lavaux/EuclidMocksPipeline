################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np

# This file containts all the information needed to define a selection

selection_tag ='lut'         # this defines a name for a specific selection

# selection keys can contain 'extinction','visibilitymask','fluxcut','central','satellite'
# 'visibilitymask' and 'fluxcut' cannot be present together
# 'central' and 'satellite' cannot be present together
selection_keys=['lookup']   

# lookup table section
lookup_table_fname   = 'LookUpTable/lookup_table_15by15.fits'
lookup_exptime_fname = None # 'LookUpTable/sc8_expmap_G.fits'
lookup_noise_fname   = 'LookUpTable/sc8_varmap_G.fits'
lookup_Nside         = 4096

# flux limit
selection_logflux_limit = -100. # log of flux limit for 'fluxcut' and 'extinction'

# visibility mask
selection_VM_fname = None       # name of visibility mask file giving flux limit in every pixel
selection_VM_res = 2048         # needed resolution for the visibility mask

# extinction map, filename starts from input.outdir
extinctionmap_fname = 'ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
extinctionmap_field = 2         # this is fine for Planck map
extinctionmap_res = 2048        # this allows to resample the map
extinction_curve = 'standard'   # THIS PART SHOULD BE EXPANDED

