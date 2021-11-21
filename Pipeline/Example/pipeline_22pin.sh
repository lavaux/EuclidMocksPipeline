#!/bin/bash
#SBATCH -J Mocks
#SBATCH -p euclid
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=240000
#SBATCH -o batch-%j.out
#SBATCH -e batch-%j.err
#SBATCH --time=720

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/u/pmonaco/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/u/pmonaco/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/u/pmonaco/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/u/pmonaco/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate py3
cd /euclid_data/pmonaco/EuclidMocks/Pipeline

# creates all the needed input files
python -u createInput_22pin.py base     None  None  False
python -u createInput_22pin.py MWmit    MWext MWext False
python -u createInput_22pin.py MWunm    MWext None  True
python -u createInput_22pin.py lutbase  lut0  lut0  False
python -u createInput_22pin.py exptimit lutet lutet False
python -u createInput_22pin.py exptiunm lutet None  True
python -u createInput_22pin.py lutmwmit lutmw lutmw False
python -u createInput_22pin.py lutmwunm lutmw None  True
python -u createInput_22pin.py noisemit lutns lutns False
python -u createInput_22pin.py noiseunm lutns None  True
python -u createInput_22pin.py allutmit lut2  lut2  False
python -u createInput_22pin.py allutunm lut2  None  True

# process 22 pinocchio light cones to produce 22 galaxy catalogs
python -u createSDHODfromPinocchio.py input_22pinbase 1 22

# process the base configuration without lookup table
# (NB: the configuration file of 2PCF is still to be written)
python -u createRandom.py input_22pinbase
python -u dndz.py input_22pinbase
python -u writeCatalogs4LE3.py input_22pinbase
python -u createPKScripts.py input_22pinbase
#python -u create2PCFScripts.py input_22pinbase

# process the base configuration with lookup table
python -u createRandom.py input_22pinlutbase
for i in {1..22}
do
    python -u createSelection.py input_22pinlutbase $i
done
python -u createSelection.py input_22pinlutbase random
python -u dndz.py input_22pinlutbase
python -u writeCatalogs4LE3.py input_22pinlutbase
python -u createPKScripts.py input_22pinlutbase
# python -u create2PCFScripts.py input_22pinlutbase

# process all the other configurations, with realistic selections
# systematics applied to base catalogs
for tag in MW noise expti lutmw allut
do
    echo $tag

    # data selection is created only once with unmit
    for i in {1..22}
    do
	python -u createSelection.py input_22pin"$tag"unm $i
    done
    # random selection is created for mit
    python -u createSelection.py input_22pin"$tag"mit random
    # while unmitigated case requires to create a new random
    python -u createRandom.py input_22pin"$tag"unm
    
    python -u dndz.py input_22pin"$tag"mit
    
    for case in unm mit
    do
        python -u dndz.py input_22pin"$tag""$case"
        python -u writeCatalogs4LE3.py input_22pin"$tag""$case"
        python -u createPKScripts.py input_22pin"$tag""$case"
#        python -u create2PCFScripts.py input_22pin"$tag""$case"
    done
done    

