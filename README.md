# Pipeline for processing dark matter halo and galaxy mocks in Euclid

## Starting point:

Build and install python package:
```
python2 setup.py install
```

Query 8614 of the Flagship mock, obtained as follows:
```
SELECT `ra_gal`, `dec_gal`, `x_gal`, `y_gal`, `z_gal`, `vx_gal`, `vy_gal`, `vz_gal`, `true_redshift_gal`, `observed_redshift_gal`, `vrad_gal`, `logf_halpha_model3_ext`, `logf_halpha_model1_ext`, `euclid_vis`, `euclid_nisp_h`, `euclid_nisp_j`, `euclid_nisp_y`, `halo_lm`, `kind`, `halo_id`, `galaxy_id`
        FROM flagship_mock_1_8_4_s 
        WHERE `true_redshift_halo`>=0.8 
        AND `true_redshift_halo`<2.0
        AND `halo_lm` >=11
        AND `z_gal` >= 0
        AND `y_gal` >= 0
        AND `x_gal` >= 0
```

## Directory Tree:
* [root]/Pipeline                        # code
* [root]/Products
*               /RawCatalogs            # Cosmohub query with various footprints applied
*               /Footprints             # survey footprints (see below)
*               /ExtinctionMaps         # extinction and reddening maps from various sources
*               /SDHOD                  # Stripped-down HOD versions
*               /GalaxyCatalogs         # complete galaxy catalogs
*               /RandomCatalogs         # complete random catalogs
*               /Selections             # selections to be applied to galaxy catalogs
*               /Catalogs4LE3           # catalogs readable by LE3 PFs, in redshift slices
*               /NumberCounts           # galaxy counts
*               /Cls                    # angular maps and angular clustering measurements
*               /PkParams               # parameter files for PK code
*               /PkScripts              # scripts to run the PK code
*               /Pks                    # measured PKs
*               /Plots                  # various plots

## Definitions:

* input file: most scripts require an input file that contains all
  needed information for ALL the scripts. The name of the input file
  (without .py) is provided to the scripts as a command line argument.

* Survey footprint: a fits file with the fields FOOTPRINT_G,
  containing boolean healpy mapx that are True/False in the sky pixels
  that are/are not observed, in galactic coordinates, and a redshift
  range of validity. It is defined by a footprint_tag in the input
  file. The original Flagship octant is defined by footprint_tag=None.
  The footprint is read by the read_footprint() function defined in
  the input file, that returns the healpy resolution of the map, the
  redshift range in which the footprint is defined, the fraction of
  the sky covered by the survey and the footprint itself.

  NOTE: it is assumed below that sky coordinates are ALWAYS galactic.

* Master catalog: a catalog containing information on all dark matter
  halos and its central/satellite galaxies, coming from a CosmoHub
  query, or from applying a footprint on a Cosmohub query.

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
  - visibilitymask: a boolean visibility mask is applied to the catalog.
  A selection is defined by a tag specified in the input file, and
  requires a separate sel_input_[tag].py file to be created. It can be
  applied to a galaxy or a random catalog. The input file specified
  two selection tags, selection_data_tag for the galaxy catalog and
  selection_random_tag for the random.

* Random catalog: it is created by (randomly) replicating the data
  catalog alpha times in redshift and flux, and by distributing
  galaxies uniformly on the sky (within the footprint). It assumed to
  be unmodulated on the sky. It can follow a specific number density
  if required **TO BE IMPLEMENTED**.

* Catalogs for LE3: these are extracted from galaxy and random
  catalogs, give galaxies in specified redshift bins, and are readable
  by LE3 codes.

## Preliminary scripts:

* createIndicesForSats.py: must be run on a master catalog to create
  pointer of satellite galaxies to their main halo in that catalog.
  Needed to create the SDHOD. Must be edited to type in the name of
  the master catalog.

* createFullOctantFootprint.py: creates a survey footprint for the
  Flagship octant.

* createFootprint.py: this is an example script for creating a survey
  footprint and writing its file. It creates a small 100 square deg
  field in a position of intermediate galactic extinction, for testing
  purposes.

* applyFootprintToMasterCatalog.py [input file]: it applies a
  footprint to the master catalog, extracting another (smaller) master
  catalog.

## Galaxy catalogs from master:

* extractGalaxyCatalogFromMaster.py [input file]: extracts a galaxy
  catalog from a master catalog by applying the standard flux limit,
  according to the model. Run it separately for each model.

## Stripped-Down HOD from master:

* createSmoothHOD.py [input file]: based on a master file, it measures
  the HOD curves and smooths them in redshift.

* createSDHOD_Catalog.py [input file]: it creates a galaxy catalog by
  populating the halos of the master catalog with galaxies, using the
  SDHOD.

## Random catalog

* createRandom.py [input file]: it creates a random for the galaxy
  catalog, with abundace alpha times the data catalog, alpha as
  specified in the input catalog. TO BE IMPLEMENTED: on request it can
  follow a specified galaxy number density.

## Selection

* createSelection.py [input file] [random]: creates a selection for a
  galaxy or a random catalog; the second option is chose if "random"
  or "1" is provided when calling the script. The parameters of the
  selection are specified in the sel_input_[input.selection_tag].py,
  where selection_tag is given in the input file.


## Measurements of galaxy catalogs:

* dndz.py [input file]: given a data catalog and a selection, it
  measures the number density of galaxies (number of galaxies brighter
  than the flux limit per sq deg per redshift interval as a function
  of redshift) in the catalog, needed to create the random catalog. It
  also provides number densities smoothed in redshift.

* writeCatalogs4LE3.py [input file]: given a selection, it writes
  galaxy and random catalogs in redshift shells, in files that are
  readable by LE3 processing functions.

* createPkScripts.py [input file]: creates the parameter files and the
  scripts needed to run the LE3 PK code.

* numbercounts.py [input file]: given a data catalog and a selection,
  it measures the differential number counts (number of galaxies per
  sq deg per redshift interval per magnitude as a function of flux) in
  redshift intervals.

* maps_and_cls.py [input file]: it builds the density contrast maps
  data and random galaxy catalogs, read from LE3 files, and measures
  their angular clustering with anafast.

## Plots

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

* angular_map.py [filename] [zmin] [zmax]: reads the galaxy or random
  catalog specified in the command line argument and plots its angular
  map of density contrast (with Nside=256) in the specified redshift
  interval. Very useful for quick checks of what a catalog is like.


