#!/usr/bin/env python3

import os
import numpy as np
from argparse import ArgumentParser
import hist as hi
import coffea.util
import hjson


def round_to_4(x):
    return round(x, 3-int(np.floor(np.log10(abs(x)))))


parser = ArgumentParser(
    description="Calculate factors to stitch MC samples from an already "
    "produced histogram (currently only 1D histograms supported)")
parser.add_argument(
    "plot_config", help="Path to a configuration file for plotting, where "
    "the 'MC_bkgd' contain the samples to be stitched (including the "
    "inclusive), and the 'Data' specifies the inclusive sample(s). Can "
    "also include a rebinning to be performed on histogram before "
    "calculating factors")
parser.add_argument(
    "histfile", help="Coffea file with a single histogram produced by "
    "select_events.py with the binning of the files to stitched")
args = parser.parse_args()

with open(args.plot_config) as f:
    config = hjson.load(f)

print("Processing {}".format(args.histfile))
srcdir = os.path.dirname(args.histfile)
if os.path.exists(os.path.join(srcdir, "hists.json")):
    with open(os.path.join(srcdir, "hists.json")) as f:
        histmap = {tuple(k): v for k, v in zip(*hjson.load(f))}
    histmap_inv = dict(zip(histmap.values(), histmap.keys()))
    histkey = histmap_inv[os.path.relpath(args.histfile, srcdir)]
else:
    histmap = None
    histkey = None
if args.plotdir:
    os.makedirs(args.outdir, exist_ok=True)
    namebase, fileext = os.path.splitext(os.path.basename(args.histfile))
hist = coffea.util.load(args.histfile)
dsaxis = hist.axis("dataset")

data_hist = hist.integrate("dataset", config["Data"])
mc_hist = hist.integrate("dataset", config["Labels"])

dense_axes = [
    ax
    for ax in hist.axes
    if isinstance(ax, (hi.axis.Regular, hi.axis.Variable, hi.axis.Integer))
]
if len(dense_axes) != 1:
    raise ValueError("Can only calculate stitching factors for 1d histograms!")
dense = dense_axes[0]
mc_preped = mc_hist.project(dense)
data_preped = data_hist.project(dense)
sfs = (data_preped.values(flow=True) / mc_preped.values(flow=True))
print([round_to_4(sf) for sf in sfs.tolist()[1:-2]])
