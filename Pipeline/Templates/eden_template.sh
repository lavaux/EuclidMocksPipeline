#!/bin/bash
#SBATCH -J PK
#SBATCH -p euclid
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=240000
#SBATCH -o PK-%j.out
#SBATCH -e PK-%j.err
#SBATCH --time=6:00:00

# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$$(CONDA_REPORT_ERRORS=false '/u/pmonaco/anaconda3/bin/conda' shell.bash hook 2> /dev/null)"
if [ $$? -eq 0 ]; then
    \eval "$$__conda_setup"
else
    if [ -f "/u/pmonaco/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/u/pmonaco/anaconda3/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        \export PATH="/u/pmonaco/anaconda3/bin:$$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<

conda deactivate

singularity shell -s /bin/bash /euclid_data/lodeen/EDEN_2.0_cvmfs.img
conda deactivate
###source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.0/etc/profile.d/euclid.sh
source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.1/bin/activate
ELogin.sh


