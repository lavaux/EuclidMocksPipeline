
# This script creates the input files for a set of options

from string import Template
from shutil import copyfile
import numpy as np
import sys

if len(sys.argv)<5:
    print("Usage: python {} tag selection_data_tag selection_random_tag apply_dataselection_to_random")
    sys.exit(0)

print("tag:                           {}".format(sys.argv[1]))
print("selection_data_tag:            {}".format(sys.argv[2]))
print("selection_random_tag:          {}".format(sys.argv[3]))
print("apply_dataselection_to_random: {}".format(sys.argv[4]))

f=open('Templates/input_22pin_template.py',mode='r')
myinput = Template( f.read() )
f.close

if sys.argv[2]=='None':
    mySDT=None
else:
    mySDT="'"+sys.argv[2]+"'"

if sys.argv[3]=='None':
    mySRT=None
else:
    mySRT="'"+sys.argv[3]+"'"

f=open('input_22pin{}.py'.format(sys.argv[1]),mode='w')
f.write(myinput.substitute(SDT = mySDT, SRT = mySRT, ADTR = sys.argv[4]))
f.close()

print("Written file input_22pin{}.py".format(sys.argv[1]))
