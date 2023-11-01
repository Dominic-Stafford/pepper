#!/usr/bin/env python3

import os
import sys
from argparse import ArgumentParser

import uproot

from pepper import HistCollection


parser = ArgumentParser(
    description="Generate a ROOT file containing efficiency histograms needed "
    "for computing b-tagging scale factors. This requires a specific "
    "histogram to have been computed, see the btageff histogram in "
    "example/hist_config.json")
parser.add_argument(
    "histsfile", help="Output Pepper histogram JSON file inside the directory "
    "of the output histograms and usually named 'hists.json'")
parser.add_argument("output", help="Output ROOT file")
parser.add_argument(
    "--cut", default="HasJets", help="Name of the cut before the b-tag "
                                     "requirement. (Default 'HasJets')")
parser.add_argument(
    "--histname", default="btageff", help="Name of the b-tagging efficiency "
                                          "histogram. (Default 'btageff')")
parser.add_argument(
    "--central", action="store_true",
    help="Only output the central efficiency and no additional ones for "
         "systematic variations")
args = parser.parse_args()

if os.path.exists(args.output):
    if input(f"{args.output} exists. Overwrite [y/n] ") != "y":
        sys.exit(1)

with open(args.histsfile) as f:
    hists = HistCollection.from_json(f)

with uproot.recreate(args.output) as f:
    full_hist = hists.load({"cut": args.cut, "hist": args.histname})
    if "sys" in [ax.name for ax in full_hist.axes]:
        hists = {
            sysname: full_hist[{"sys": sysname}]
            for sysname in full_hist.axes["sys"]
        }
    else:
        hists = {"nominal": full_hist}
    for sysname, hist in hists.items():
        hist = hist.project("flav", "pt", "abseta", "btagged")
        eff = hist[{"btagged": "yes"}] / hist[{"btagged": sum}].values()
        if sysname == "nominal":
            f["central"] = eff
            print("Nominal scale factors:")
            print(eff.values())
        elif args.central:
            continue
        elif any(sysname.endswith(x) for x in (
                "XS_down", "XS_up", "lumi_down", "lumi_up")):
            # These aren't shape uncertainties and also do not have much effect
            # on the efficiency
            continue
        else:
            f[sysname] = eff
