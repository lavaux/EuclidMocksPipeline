from __future__ import print_function
import pkg_resources
import argparse as ap
import re
import os

def run():
    input_template = pkg_resources.resource_string(__name__,
                                                   "templates/input.py")

    parser = ap.ArgumentParser(
        description="generate config file for Euclid mock")
    parser.add_argument("--lf_model", type=int, default=1)
    parser.add_argument("--type",
                        choices=["sdhod", "flagship", "pinocchio", "box"],
                        type=str,
                        default="sdhod")
    parser.add_argument("--rsd",choices=["false","true"],default="false")
    parser.add_argument("--outdir",type=str,default=os.getcwd())
    parser.add_argument("-o",type=str)

    args = parser.parse_args()

    replacements = {
        "OUTDIR": args.outdir,
        "FOOTTAG": "None",
        "LFMODEL": repr(args.lf_model),
        "CATTYPE": args.type,
        "SHUFFLE": "False",
        "SELDATA": "None",
        "SELRAND": "None",
        "RSDFLAG": repr(args.rsd == "true")
    }

    print("Output directory is {}".format(args.outdir))

    for k, v in replacements.items():
        input_template = re.sub(k, v, input_template)

    with open(args.o, mode="wt") as f:
        f.write(input_template)
