################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from string import Template
from shutil import copyfile
import numpy as np
import filenames
from glob import glob

import sys
if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running createCM-PKScript.py with {}".format(sys.argv[1]))

if input.cat_type is not 'pinocchio':
    print("ERROR: I need to process a set of catalogs, cat_type should be pinocchio")
    sys.exit(0)

if input.pinocchio_last_run is None:
    print("ERROR: I need to process a set of catalogs, pinocchio_last_run should be specified")
    sys.exit(0)
    
# Open a file: params
f = open('Templates/parameters_CM-PK_template.ini',mode='r')
params = Template( f.read() )
f.close()

r1=input.pinocchio_first_run
r2=input.pinocchio_last_run
toprocess=range(r1,r2+1)

# script
scriptfname=filenames.estimator_script(input,'CM-PK',0)
copyfile("Templates/eden_template.sh",scriptfname)
print("Writing {}".format(scriptfname))
script=open(scriptfname,"a")

for z1, z2 in input.finalCatZShell:

    # parameter file
    paramfname=filenames.estimator_params(input,'CM-PK',z1,z2)
    print("writing %s"%paramfname)

    if input.Lbox is None:
        clbox='true'
        lbox=0.
    else:
        clbox='false'
        lbox=input.Lbox
    f = open(paramfname, "w")
    f.write(params.substitute(GRID   = input.ngrid,
                              PKDIR  = '',
                              OUTPUT = 'CM-PK/Measures/'+filenames.exclude_dir(filenames.estimator_measure(input, 'CM-PK', z1, z2, None)),
                              LIST   = 'CM-PK/Params/'+filenames.exclude_dir(filenames.estimator_list(input, 'CM-PK', z1, z2)) ))
    f.close()

    # json list
    json_fname=filenames.estimator_list(input,'CM-PK',z1,z2)
    string='['
    for myrun in range(r1,r2+1):
        fname = glob(filenames.estimator_measure(input,'PK',z1,z2,myrun)+"/EUC_LE3_GCL_PK_*0Z_1D.fits")[0]
        string += '"'+fname[len(input.project):]+'"'
        if myrun is not r2:
            string += ', '
    string += ']'
    print("writing {}".format(json_fname))
    f=open(json_fname,"w")
    f.write(string)
    f.close()


    script.write("\n")
    script.write("mkdir -p "+filenames.estimator_measure(input,'CM-PK',z1,z2)+"\n")
    script.write("E-Run LE3_GC_PowerSpectrumCovariance  LE3_GC_ComputePowerSpectrumCovariance --log-level=DEBUG --parfile=CM-PK/Params/"+filenames.exclude_dir(filenames.estimator_params(input,'CM-PK',z1,z2))+" --workdir="+input.project+"\n\n")

script.write("echo '### {} DONE! ###'\n".format(scriptfname))
script.close()

print("# DONE!")

