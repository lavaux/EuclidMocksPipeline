#!/bin/bash
#SBATCH -J Mocks
#SBATCH -p euclid
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=240000
#SBATCH -o batch-%j.out
#SBATCH -e batch-%j.err
#SBATCH --time=360

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

# python -u applyFootprintToMasterCatalog.py input_SC8
# python -u extractGalaxyCatalogFromMaster.py input_SC8
python -u dndz.py input_SC8
python -u createRandom.py input_SC8
python -u writeCatalogs4LE3.py input_SC8
python -u createPkScripts.py input_SC8

for tag in 0 ns et mw 1 2
do

    python -u createSelection.py input_SC8_lut$tag
    python -u dndz.py input_SC8_lut$tag
    python -u createSelection.py input_SC8_lut$tag 1
    python -u writeCatalogs4LE3.py input_SC8_lut$tag
    python -u createPkScripts.py input_SC8_lut$tag
    python -u createRandom.py input_SC8_lut${tag}_dsrand
    python -u writeCatalogs4LE3.py input_SC8_lut${tag}_dsrand
    python -u createPkScripts.py input_SC8_lut${tag}_dsrand

done
