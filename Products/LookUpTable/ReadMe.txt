 Written By: Maxwell W Kuschel
 Updated: 20/04/2021
#
 Here is a quick description of the contents of this repository, including
 a description of each folder, as well as a description of how the binned
 lookup tables were constructed. 
#
#################
# unbinned_fits #
#################
#
 Contains the output of the pypelid_lookup.py simulation. The output is 
 organized into 53 fits files, each file corresponds to a run with a specific
 background value and number of exposures. The file names are written to show
 the background and number of exposures for each file. The files have the 
 following columns, the units are available in the header of each file, and do 
 not change.
#
# index	#   label   # description #
###################################
   0 	- id	    # The haloid
   1 	- z	    # redshift
   2 	- amp	    # SNR
   3 	- flux_Ha   # Halpha flux
   4 	- flux_N2   # [NII] flux
   5	- radius    # Galaxy radius
   6	- sky_bg    # Background
   7	- int_time  # Exposure time
#
# Naming Convention for Output
#
 run_e{extinction}_b{background}_n{nexp}_c{catalog}.csv

 **Example:**

 run_e0_b12_n4_full.fits would be

 extinction = 0
 background = 1.2
 nexp = 4
#
#################
# binned_tables #
#################
#
 Contains the binned lookup tables, created using the files in the 
 unbinned_fits folder. These files were created using the 
 outputBinning.ipynb notebook in the notebooks folder. The file names once
 again correspond to what the file contains, in this instance the numbers
 in the file name represent the two changeable bins. The data is automatically
 binned by number of exposures and background, as that is how it was simulated,
 so only redshift and Halpha bins can be canged by the user. The naming
 convention is explained using an example below
#
# Naming Convention for Output
#
 lookup_table_{N_z_bins}by{N_fHa_bins}.fits
 
 **Example**

 lookup_table_15by10.fits

 Number of z bins = 15
 Number of Halpha bins = 10
 Number of exposure bins = 4
 Number of background bins = 13
#
# How to Use the Table
#
 Each column corresponds to a specific variable we care about, 

# index	#    label     # description #
######################################
   0 	- z	       # redshift
   1 	- flux_Ha      # Halpha flux
   2 	- exp_time     # Exposure time
   3	- sky_bg       # Background
   4	- galaxy_frac  # Fraction of observable Galaxies (SNR>3.5)
   5	- num_detected # Number of galaxies detected
   6	- num_galaxies # Number of galaxies

 The values for z and flux_Ha are the middle of the bins used, for exp_time
 and sky_bg they are they exact values used in the simulation.

#
###################
# python_programs #
###################
#
 This contains the programs written to help facilitate the plotting of the 
 lookup table. As of this writing this folder is empty, as I am still 
 working on getting these programs created and organized.
#
#############
# notebooks #
#############
#
 This contains the notebooks either used to create the lookup table and
 notebooks that can be used instead of the python programs in the 
 python_programs folder. These also contain tutorials for using the 
 functions and processes defined in the python_programs folder. In 
 addition to all of this we also contain the raw notebooks used to
 create the plots in the written report, in case you wish to reuse or
 tweak those plots. Runtime warnings are available where they should apply,
 they are there because some of these processes can take 10+ minutes to 
 complete.
