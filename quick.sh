#!/bin/sh

set -e

PYTHON=$(which python3)
BASETOOL="$PYTHON -m euclid_obssys.tool"
VIEWTOOL="$PYTHON -m euclid_obssys.view" 

FOOTPRINT="createFootprint"
#FOOTPRINT="createFullOctantFootprint"
CONFIG=myconf.py

if false; then
$BASETOOL generateConfig $CONFIG --footprintTag=100sqdeg --repodir=$(pwd)/Repo --projectdir=TestProject
## these are typically not used in a production pipeline but are run only once
$BASETOOL $FOOTPRINT                # create catalog footprint

# these are the main scripts to be used in a pipeline
fi


$BASETOOL applyFootprintToMaster $CONFIG   # produces a smaller (rotated?) master
exit 0
$BASETOOL extractGalaxyCatalogFromMaster $CONFIG  # extracts a galaxy catalog from the master
$BASETOOL createIndicesForSats $CONFIG      # must be run on a master to measure the SDHOD
$BASETOOL createSmoothHOD $CONFIG     # measures the HOD from a flagship master


$BASETOOL createSDHOD_Catalog $CONFIG              # creates an SDHOD galaxy catalog from a halo catalog

$BASETOOL dN_dZ $CONFIG                            # measures the n(z) of a galaxy catalog
$BASETOOL createRandom $CONFIG                  # creates the random (before selection)

$BASETOOL createSelection --use-data=False $CONFIG       # creates a selection for random
$BASETOOL createSelection $CONFIG  # creates a selection for galaxy

$BASETOOL writeCatalogs4LE3 $CONFIG               # writes the catalog in redshift bins

$BASETOOL numbercounts $CONFIG                    # number counts for redshift bins

$VIEWTOOL angular_map Products/GalaxyCatalogs/catalogs_hodcat_8614_100sqdeg_m1_cMdiemer19.fits Products/Selections/catalogs_data_lutmw_catalogs_hodcat_8614_100sqdeg_m1_cMdiemer19.fits 0 1 --Nside 1024 --output angular_default.pdf

$VIEWTOOL visualizeHOD $CONFIG                    

#maps_and_cls.py input                    # measurement of angular clustering  WORK IN PROGRESS
$VIEWTOOL dn_dz $CONFIG

exit 0

#createSDHODfromPinocchio.py input        # creates a set of galaxy catalog from pinocchio lightcones

createPkScripts.py input                 # creates the script and param files for PK

# plotting scritps
plot_numbercounts.py input
plot_pk.py input
plot_cls.py input                   WORK IN PROGRESS
angular_map.py filename selection zmin zmax   # quick plot of an angular map

# tools
sdhod.py
matchPinocchioHaloMasses.py
NFW.py
pozzettiLF.py
ReadPinocchio.py
