# Pipeline for processing dark matter halo and galaxy mocks in Euclid
# V0.3

## Starting point:
Query 8614 of the Flagship mock, obtained as follows:
SELECT `ra_gal`, `dec_gal`, `x_gal`, `y_gal`, `z_gal`, `vx_gal`, `vy_gal`, `vz_gal`, `true_redshift_gal`, `observed_redshift_gal`, `vrad_gal`, `logf_halpha_model3_ext`, `logf_halpha_model1_ext`, `euclid_vis`, `euclid_nisp_h`, `euclid_nisp_j`, `euclid_nisp_y`, `halo_lm`, `kind`, `halo_id`, `galaxy_id`
        FROM flagship_mock_1_8_4_s 
        WHERE `true_redshift_halo`>=0.8 
        AND `true_redshift_halo`<2.0
        AND `halo_lm` >=11
        AND `z_gal` >= 0
        AND `y_gal` >= 0
        AND `x_gal` >= 0

## Repository and project: 
Some files of general use are in the Repo/ directory, while files
produced by the pipeline are found in subdirectories of a project
directory, defined in the input file (see below). The script
createProjectDirectoryTree.py creates a project tree directory. To run
it:
> python createProjectDirectoryTree.py [project name]

## Directory Tree:
[root]/Pipeline                        # code
[root]/Repo
               /RawCatalogs            # Cosmohub query with various footprints applied
               /Footprints             # survey footprints (see below)
               /ExtinctionMaps         # extinction and reddening maps from various sources
               /SDHOD                  # Stripped-down HOD versions
	       /LookUpTable            # lookup table to bypass realistic galaxy selection
	       /SelectionInputs        # selection input files
[root]/[project]
               /GalaxyCatalogs         # complete galaxy catalogs
               /RandomCatalogs         # complete random catalogs
               /Selections             # selections to be applied to galaxy catalogs
               /Catalogs4LE3           # catalogs readable by LE3 PFs, in redshift slices
               /NumberCounts           # galaxy counts
               /Cls                    # angular maps and angular clustering measurements
               /Plots                  # various plots
for each estimator (2PCF, 3PCF, CM-2PCF, PK, BK, CM-PK):
               /[estimator]/Params     # parameter files (and product lists for CM-*)
               /[estimator]/Scripts    # scripts to run the estimator
               /[estimator]/Measures   # measured quantities

## Definitions:

* input file: most scripts require an input file that contains all
  needed information for ALL the scripts. The name of the input file
  (without .py) is provided to the scripts as a command line argument.

* Survey footprint: a fits file with the fields FOOTPRINT_G,
  containing boolean healpy mapx that are True/False in the sky pixels
  that are/are not observed, in galactic coordinates, and a redshift
  range of validity. It is defined by a footprint_tag in the input
  file. The original Flagship octant is defined by footprint_tag=None.
  The footprint is read by the input.read_footprint() function defined
  in the input file, that returns the healpy resolution of the map,
  the redshift range in which the footprint is defined, the fraction
  of the sky covered by the survey and the footprint itself.

  NOTE: it is assumed below that sky coordinates are ALWAYS galactic.

* Master catalog: a catalog containing information on all dark matter
  halos and their central/satellite gal

axies, coming from a CosmoHub
  query, or from applying a footprint on a Cosmohub query.x

* SDHOD: the pipeline can create a stripped-down HOD version of the
  Flagship HOD, that is used to populate sets of DM halos in the
  lightcone, both from Flagship and Pinocchio.

* Galaxy catalog: a catalog that contains only galaxies subject to a
  basic selection, like Halpha flux > 2e-16 erg/cm^2/s. The basic
  selection is based on a flux_key, that may be logf_halpha_model1_ext
  (model '1'), logf_halpha_model3_ext (model '3'), euclid_nisp_H
  (model '9') or halo_lm (model '0'). The threshold criterion is based
  on logflux_limit, specified in the input file. A galaxy catalog
  contains the following fields (with the same meaning as the Flagship
  catalog): x_gal, y_gal, z_gal, ra_gal, dec_gal, kind,
  true_redshift_gal, observed_redshift_gal, halo_lm, id, halo_id,
  [flux_key], sh_[flux_key]. This last column contains the same
  fluxes, shuffled among galaxies at fixed redshift to suppress
  luminosity-dependent bias.

* Selection: a boolear array of the same size of its reference galaxy
  or random catalog, that selects galaxies according to:
  - fluxcut: a simple cut in a reference variable (may be Halpha flux,
    H band magnitude, halo mass);
  - central / satellite: central or satellite galaxies;
  - extinction: galaxies are selected accorting to their extincted
    flux;
  - visibilitymask: galaxies are selected in flux using a healpy map
    giving the flux limit in each pixel;
  - lookup: a lookup table, giving the detection probability as a
    function of flux, redshift, exposure time, noise level and MW
    extinction, is used to determine what galaxies are detected.

  A selection is defined by a tag specified in the input file, and
  requires a separate sel_input_[tag].py file to be created. It can be
  applied to a galaxy or a random catalog. The input file specifies
  two selection tags, selection_data_tag for the galaxy catalog and
  selection_random_tag for the random. For systematics, applying the
  same selection to the data and to the corresponding random catalog
  will give a perfect mitigation of the systematics. Because selection
  changes the number density, the case of no mitigation will be
  obtained by creating a random based on the data catalog after
  applying the selection, not by using the parent random catalog.

* Random catalog: it is created by (randomly) replicating the data
  catalog (with or without applying a selection) alpha times in
  redshift and flux, and by distributing galaxies uniformly on the sky
  (within the footprint). Its number density thus depends only on
  redshift. A random with modulated number density is obtained by
  applying a selection to it. **TO BE IMPLEMENTED**: Instead of having
  exactly the same fluctuations in number density (measured in bins of
  0.01 width in redshift), the random can follow a smoothed version of
  the data number density.

* Catalogs for LE3: these are extracted from galaxy and random
  catalogs, with selections applied as needed, in a sequence of
  redshift bins as specified in the input file, and are readable by
  LE3 codes.

* Pinocchio light cones: the pipeline can process a set of pinocchio
  light cones. In this case it is possible to compute the galaxy
  number density as the average over a set of mocks, and a random will
  follow this averaged number density.

## Preliminary scripts:

* createProjectDirectoryTree.py: it creates the directory tree for a
  project, and copies the selection input files from the repository to
  the Selections/ directory where they are expected to be found;
  indeed, one may want to edit them.

* createIndicesForSats.py: must be run on a master catalog to create
  pointers of satellite galaxies to their main halo in that catalog.
  They are needed to create the SDHOD. This script must be edited to
  type in the name of the master catalog. Most likely, it will be run
  only once on the main master.

* createFullOctantFootprint.py: creates a survey footprint for the
  Flagship octant.

* createFootprint.py: this is an example script for creating a survey
  footprint and writing its file. It creates a small 100 square deg
  field in a position of intermediate galactic extinction, for testing
  purposes.

* applyFootprintToMasterCatalog.py [input file]: it applies a
  footprint to the master catalog, extracting another (smaller) master
  catalog that goes in the Repo/ directory.

## Galaxy catalogs from master:

* extractGalaxyCatalogFromMaster.py [input file]: extracts a galaxy
  catalog from a master catalog by applying the standard flux limit,
  according to the model. It must be run separately for each model.

## Creation of the Stripped-Down HOD from master:

* createSmoothHOD.py [input file]: based on a master file, it measures
  the HOD curves and smooths them in redshift. It should be run on the
  full octant, there is no reason to run it many times for a specified
  model.

## Galaxy catalog created from master DM halos using SDHOD

* createSDHODCatalog.py [input file]: it creates a galaxy catalog by
  populating the halos of the master catalog with galaxies, using the
  SDHOD, including standard recipes to distribute satellites in halos.

## Galaxy catalog created from pinocchio light cones using SDHOD

* createSDHODfromPinocchio.py [input file] [first run] [last run]: it
  creates a set of galaxy catalog by populating the halos of pinocchio
  light cones with galaxies, using the SDHOD, including standard
  recipes to distribute satellites in halos. A calibration of halo
  masses is applied, to match them to the Flagship catalog either in
  abundance or in clustering amplitude. Here the first and last runs
  are not taken from the input file but provided as an argument, to
  guaranteee more flexibility.

## Random catalog

* createRandom.py [input file]: it creates a random for the galaxy
  catalog (see the explanation above), with abundance alpha times the
  data catalog, alpha as specified in the input catalog. If
  apply_dataselection_to_random is set to True, it applies a selection
  to data catalog before replicating it, so as to obtain a random that
  follows the catalog after selection but has a mean density constant
  in time (resulting in no mitigation of the systematics). In case of
  a set of mocks, a sequence of integers is created, as long as the
  number of mocks, that sums up to the alpha value specified in the
  input file and says how many times a mock will be replicated. This
  implies that the number of random galaxy is not exactly the
  (average) number of data galaxies times alpha, but this is not a
  problem for the estimators. **TO BE IMPLEMENTED**: On request it can
  follow a smoothed version of the galaxy number density.

## Selection

* createSelection.py [input file] [run number / random]: creates a
  selection for a galaxy or a random catalog. The second option is set
  to the run number when a set of pinocchio mocks is processed
  (selections must be produced catalog by catalog) or to the string
  'random' to create a selection for the random. The parameters of the
  selection are specified in the sel_input_[input.selection_tag].py,
  where selection_tag is given in the input file.

## Measurements of galaxy catalogs:

* dndz.py [input file]: given a galaxy catalog (possibly with a
  selection), it measures the number density of its galaxies in
  redshift bins of delta z=0.01. It also provides number densities
  smoothed in redshift.

* writeCatalogs4LE3.py [input file]: it writes galaxy and random
  catalogs (with selections as requested) in a sequence redshift bins
  as specified in the input, in files that are readable by LE3
  processing functions.

* create[estimator]Scripts.py [input file]: creates the parameter
  files and the scripts needed to run the LE3 estimator code.
  Currently provided scripts are for PK, 2PCF, CM-PK.

* numbercounts.py [input file]: given a data catalog and a selection,
  it measures the differential number counts (number of galaxies per
  sq deg per redshift interval per magnitude as a function of flux) in
  redshift intervals.

* maps_and_cls.py [input file]: it builds the density contrast maps
  data and random galaxy catalogs, read from LE3 files, and measures
  their angular clustering with anafast.
  THIS IS STILL WORK IN PROGRESS

## Plots

All the scripts below should be considered as example script to
produce the wanted plots.

* visualizeHOD.py [input file]: plots the Ncen(Mh,z) and Nsat(Mh,z)
  curves for the SDHOD.

* plot_dndz.py [input file]: plots the number density of galaxies as a
  function of redshift, also separated into centrals and satellites,
  and compares it to the relative Pozzetti model if relevant. It is
  easy to edit the script to plot two sets of curves.

* plot_numbercounts.py [input file]: plots the differential number
  counts as a function of flux for specified redshift intervals, also
  separated into centrals and satellites, and compares it to Pozzetti
  model if relevant.

* plot_pk.py [input file]: plot the measured PKs (monopole, quadrupole,
  exadecapole) in redshift intervals. It is easy to edit the script to
  plot two sets of curves.

* plot_cls.py [input file]: plot the measured angular clustering in redshift
  intervals. It is easy to edit the script to plot two sets of curves.
  THIS IS STILL WORK IN PROGRESS


## Working Example

Please see the Example/ directory in Pipeline/ for a complete
description of a full working example.
