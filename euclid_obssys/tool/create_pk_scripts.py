################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig
import pkg_resources


def batch_script(iterator, name_gen, batch):
    count = 0
    batch_count = 0
    scriptfname=None
    
    for i in iterator:
        if count == 0:
            scriptfname=name_gen(batch_count)
            script_template = pkg_resources.resource_string(
                "euclid_obssys", "templates/eden_template.sh"
            )
            with open(scriptfname, mode="wt") as new_script:
                new_script.write(script_template.decode('utf-8'))
            print(f"Writing {scriptfname}")
            batch_count += 1

        with open(scriptfname,"a") as script:
            yield script,i
            count = (count + 1) % batch
            if count == 0:
                script.write("echo '### {} DONE! ###\n'".format(scriptfname))

    with open(scriptfname,"a") as script:
        script.write(f"echo '### {scriptfname} DONE! ###'\n")


@register_tool
def createPkScripts(config: str) -> None:
    """Create scripts to execute LE3 power spectrum estimator on the mock catalogs

    Args:
        config (str): Pipeline configuration file
    """
    from euclid_obssys.disk import DefaultCatalogWrite
    from string import Template
    from shutil import copyfile
    import itertools
    import numpy as np
    from .. import filenames
    import sys

    input = readConfig(config)

    print(f"# Running createPKScript.py with {config}")

    # Open a file: params

    pk_template = pkg_resources.resource_string(
        "euclid_obssys", "templates/parameters_PK_template.ini"
    ).decode("utf-8")
    params = Template( pk_template )

    count=0
    nscript=0

    r1=input.pinocchio_first_run
    r2=input.pinocchio_last_run
    if input.cat_type is 'pinocchio':
        n1=r1
        if r2 is not None:
            n2=r2
        else:
            n2=n1
    else:
        n1=n2=0

    if input.cat_type is not 'pinocchio':
        toprocess=[None]
    else:
        toprocess=range(n1,n2+1)

    for script, (myrun,(z1,z2)) in batch_script(itertools.product(toprocess, input.finalCatZShell), name_gen=lambda n: filenames.estimator_script(input,'PK',n), batch=input.max_PKs_in_script):

        paramfname=filenames.estimator_params(input,'PK',z1,z2,myrun)
        print("writing %s"%paramfname)

        if input.Lbox is None:
            clbox='true'
            lbox=0.
        else:
            clbox='false'
            lbox=input.Lbox
        with open(paramfname, "w") as f:
            f.write(params.substitute(GRID   = input.ngrid,
                                    DATA   = filenames.exclude_dir(filenames.LE3_data(input, z1, z2, myrun)),
                                    RANDOM = filenames.exclude_dir(filenames.LE3_random(input, z1, z2)),
                                    OUTPUT = 'PK/Measures/'+filenames.exclude_dir(filenames.estimator_measure(input, 'PK', z1, z2, myrun)),
                                    CLBOX  = clbox, 
                                    LBOX   = lbox))

        script.write("\n")
        script.write("mkdir -p "+filenames.estimator_measure(input,'PK',z1,z2,myrun)+"\n")
        script.write("E-Run LE3_GC_PowerSpectrum  LE3_GC_ComputePowerSpectrum --log-level=DEBUG --parfile=PK/Params/"+filenames.exclude_dir(filenames.estimator_params(input,'PK',z1,z2,myrun))+" --workdir="+input.project+"\n\n")

    print("# DONE!")

