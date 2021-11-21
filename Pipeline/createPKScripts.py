################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from string import Template
from shutil import copyfile
import numpy as np
import filenames

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

print("# Running createPKScript.py with {}".format(sys.argv[1]))

# Open a file: params
f = open('Templates/parameters_PK_template.ini',mode='r')
params = Template( f.read() )
f.close()

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

for myrun in toprocess:
    for z1, z2 in input.finalCatZShell:

        paramfname=filenames.estimator_params(input,'PK',z1,z2,myrun)
        print("writing %s"%paramfname)

        if input.Lbox is None:
            clbox='true'
            lbox=0.
        else:
            clbox='false'
            lbox=input.Lbox
        f = open(paramfname, "w")
        f.write(params.substitute(GRID   = input.ngrid,
                                  DATA   = filenames.exclude_dir(filenames.LE3_data(input, z1, z2, myrun)),
                                  RANDOM = filenames.exclude_dir(filenames.LE3_random(input, z1, z2)),
                                  OUTPUT = 'PK/Measures/'+filenames.exclude_dir(filenames.estimator_measure(input, 'PK', z1, z2, myrun)),
                                  CLBOX  = clbox, 
                                  LBOX   = lbox))
        f.close()

        if count==0:

            scriptfname=filenames.estimator_script(input,'PK',nscript)
            copyfile("Templates/eden_template.sh",scriptfname)
            print("Writing {}".format(scriptfname))
            script=open(scriptfname,"a")

        script.write("\n")
        script.write("mkdir -p "+filenames.estimator_measure(input,'PK',z1,z2,myrun)+"\n")
        script.write("E-Run LE3_GC_PowerSpectrum  LE3_GC_ComputePowerSpectrum --log-level=DEBUG --parfile=PK/Params/"+filenames.exclude_dir(filenames.estimator_params(input,'PK',z1,z2,myrun))+" --workdir="+input.project+"\n\n")

        count += 1
        if count==input.max_PKs_in_script:
            count=0
            script.write("echo '### {} DONE! ###\n'".format(scriptfname))
            script.close()
            nscript+=1

script.write("echo '### {} DONE! ###'\n".format(scriptfname))
script.close()

print("# DONE!")

