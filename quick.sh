#!/bin/sh

PYTHON=$(which python3)
BASETOOL="$PYTHON -m euclid_obssys.tool"

FOOTPRINT="createFootprint"
#FOOTPRINT="createFullOctantFootprint"
CONFIG=myconf.py

if false; then
$BASETOOL generateConfig $CONFIG --footprintTag=100sqdeg --outdir=$(pwd)/Products
## these are typically not used in a production pipeline but are run only once
$BASETOOL $FOOTPRINT                # create catalog footprint

# these are the main scripts to be used in a pipeline


$BASETOOL applyFootprintToMaster $CONFIG   # produces a smaller (rotated?) master
$BASETOOL extractGalaxyCatalogFromMaster $CONFIG  # extracts a galaxy catalog from the master
$BASETOOL createIndicesForSats $CONFIG      # must be run on a master to measure the SDHOD
$BASETOOL createSmoothHOD $CONFIG     # measures the HOD from a flagship master


$BASETOOL createSDHOD_Catalog $CONFIG              # creates an SDHOD galaxy catalog from a halo catalog
$BASETOOL dN_dZ $CONFIG                            # measures the n(z) of a galaxy catalog

$BASETOOL createRandom $CONFIG                  # creates the random (before selection)

$BASETOOL createSelection $CONFIG  # creates a selection for galaxy
$BASETOOL createSelection --use-data=False $CONFIG       # creates a selection for random
fi


$BASETOOL writeCatalogs4LE3 $CONFIG               # writes the catalog in redshift bins

exit 0


createSDHODfromPinocchio.py input        # creates a set of galaxy catalog from pinocchio lightcones





numbercounts.py input                    # number counts for redshift bins

maps_and_cls.py input                    # measurement of angular clustering  WORK IN PROGRESS

createPkScripts.py input                 # creates the script and param files for PK

# plotting scritps
visualizeHOD.py input                    
plot_dndz.py input
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
