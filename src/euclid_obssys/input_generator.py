from __future__ import print_function
import pkg_resources
import argparse as ap
import re

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

    args = parser.parse_args()

    replacements = {
        "FOOTTAG": "None",
        "LFMODEL": repr(args.lf_model),
        "CATTYPE": args.type,
        "SHUFFLE": "False",
        "SELDATA": "None",
        "SELRAND": "None",
        "RSDFLAG": repr(args.rsd == "true")
    }

    for k, v in replacements.items():
        input_template = re.sub(k, v, input_template)

    print(input_template)
