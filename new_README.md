# Initial setup


## Using CONDA

Setup a typical CONDA installation with python3. Activate your environment first. You may install the requirements using:
```
python3 -m pip install -r requirements.txt
```

Then you may run the installation of the euclid observational systematics package (`euclid_obssys`) using
```
python3 setup.py install
```

Two basic tool are made available through `$PYTHON -m euclid_obssys.tool` and `$PYTHON -m euclid_obssys.view`. Each command comes with
an help message describing the input/output and the optional arguments. Please consult `quick.sh` for a complete example.

## Workflow

In all the text below, we will use `OUTPUT` to indicate the base output directory for the files. Please expand it with the appropriate directory of your
chosing.


### Download Flagship halos

For this, go to CosmoHub and use the following query on the Flagship mock database:
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

Ask for a FITS file output. Once it is ready download the file in `OUTPUT/RawCatalogs/8614.fits`. 


### Quick pipeline

A default script that runs an entire pipeline is available in `quick.sh`. You need to have downloaded the flagship halos in `Products/RawCatalogs/catalogs_8614.fits`.

### Setting up configuration

Next you want to create a configuration file for your mock catalog. The pipeline comes bundled with
a generator that can be invoked as:
```shell
python3 -m euclid_obssys.tool generateInput --type sdhod -o mytest.py  --outdir OUTPUT
```
(Remember: you have to replace `OUTPUT`  by that directory path)


### First step with the pipeline

Now we may use that configuration file with pipeline elements
```shell
python3 -m euclid_obssys.tool applyFootprintToMaster mytest.py 
```
