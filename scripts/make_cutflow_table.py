#!/usr/bin/env python3
"""
Write the cutflow out as a table, split into W pt bins.

For every dataset and every cut, gives the number of events surviving in
total and in each W pt bin. Two numbers are reported per cell: the sum of
weights, which is what pepper prints during a run, and the effective number
of events, which is what the statistical uncertainty depends on.

The per bin numbers come from the Wplus_pt_bins and Wminus_pt_bins
histograms, which are filled at every cut. Each event contributes twice,
once for its W+ and once for its W-, so the per bin columns count W bosons
rather than events and their total is twice the event count.

Usage
-----
    python scripts/make_cutflow_table.py output/hists/hists.json \\
        -o output/cutflow_table.csv
"""
import os
import csv
import json
from argparse import ArgumentParser

import numpy as np
import hist

from pepper import HistCollection

# Must match the groups used in the plotting scripts
WPT_GROUPS = [(0, 60), (60, 120), (120, 200), (200, None)]
# Histograms holding the W pt at every cut, one per charge
WPT_VARIABLES = ["Wplus_pt_bins", "Wminus_pt_bins"]


def group_label(lo, hi):
    return f"pt>{lo:.0f}" if hi is None else f"pt{lo:.0f}-{hi:.0f}"


def integrate_categories(h, dataset=None):
    """Sum over every category axis, keeping only the nominal systematic.
    Returns None if `dataset` never filled this histogram, which happens for
    the Higgs columns in a sample that has no Higgs."""
    if dataset is not None and "dataset" in h.axes.name:
        if dataset not in list(h.axes["dataset"]):
            return None
        h = h[{"dataset": dataset}]
    if "sys" in h.axes.name:
        h = h[{"sys": "nominal"}]
    for ax in list(h.axes):
        if isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory)):
            h = h.integrate(ax.name, list(ax))
    return h


def wpt_totals(hists, cutname, dataset):
    """Sum of weights and sum of squared weights per W pt group, adding the
    two W charges together. Returns None if the histograms are absent."""
    values = variances = axis = None
    for variable in WPT_VARIABLES:
        key = (cutname, variable)
        if key not in hists.keys():
            continue
        h = integrate_categories(hists.load(key), dataset)
        if h is None:
            continue
        dense = [a for a in h.axes
                 if not isinstance(a, (hist.axis.StrCategory,
                                       hist.axis.IntCategory))]
        if len(dense) != 1:
            continue
        axis = dense[0]
        values = h.values() if values is None else values + h.values()
        variances = (h.variances() if variances is None
                     else variances + h.variances())
    if axis is None:
        return None

    edges = axis.edges
    out = {}
    for lo, hi in WPT_GROUPS:
        top = edges[-1] if hi is None else hi
        keep = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= top + 1e-9)
        sumw = values[keep].sum()
        sumw2 = variances[keep].sum()
        out[group_label(lo, hi)] = (sumw, sumw2)
    return out


def main():
    parser = ArgumentParser(
        description="Write the cutflow as a table split into W pt bins")
    parser.add_argument("histfile", help="hists.json written by runproc")
    parser.add_argument("-c", "--cutflows", help="cutflows.json. Defaults to "
                        "the one next to the hists directory")
    parser.add_argument("-o", "--outfile", help="Where to write the csv. "
                        "Defaults to cutflow_table.csv beside cutflows.json")
    args = parser.parse_args()

    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    histdir = os.path.dirname(os.path.realpath(args.histfile))
    cutflow_path = args.cutflows
    if cutflow_path is None:
        cutflow_path = os.path.join(os.path.dirname(histdir), "cutflows.json")
    cutflows = {}
    if os.path.exists(cutflow_path):
        with open(cutflow_path) as f:
            cutflows = json.load(f)
    else:
        print(f"No cutflows at {cutflow_path}, the total column will be blank")

    outfile = args.outfile
    if outfile is None:
        outfile = os.path.join(os.path.dirname(histdir), "cutflow_table.csv")

    # Cut order, taken from the collection userdata where pepper stores it
    cuts = None
    if isinstance(hists.userdata, dict) and "cuts" in hists.userdata:
        cuts = list(hists.userdata["cuts"])
    if not cuts:
        cuts = sorted({k[0] for k in hists.keys()})
        print("Cut order not found in the hists file, falling back to "
              "alphabetical order")

    datasets = sorted(cutflows) if cutflows else [None]
    datasets = [d for d in datasets if d != "all"] or [None]

    group_names = [group_label(lo, hi) for lo, hi in WPT_GROUPS]
    header = (["dataset", "cut", "total_sumw", "total_eff_events"]
              + [f"{g}_{s}" for g in group_names
                 for s in ("sumw", "eff_events")])

    rows = []
    for dataset in datasets:
        for cut in cuts:
            total = ""
            if dataset is not None and dataset in cutflows:
                entry = cutflows[dataset]
                if isinstance(entry, dict) and "all" in entry:
                    entry = entry["all"]
                total = entry.get(cut, "")
            row = [dataset if dataset is not None else "all", cut,
                   f"{total:.1f}" if isinstance(total, (int, float)) else "",
                   ""]
            per_bin = wpt_totals(hists, cut, dataset)
            if per_bin is None:
                row += [""] * (2 * len(group_names))
            else:
                grand_w = sum(v[0] for v in per_bin.values())
                grand_w2 = sum(v[1] for v in per_bin.values())
                # Each event enters twice, once per W
                if grand_w2 > 0:
                    row[3] = f"{grand_w**2 / grand_w2 / 2:.0f}"
                for g in group_names:
                    sumw, sumw2 = per_bin[g]
                    n_eff = sumw**2 / sumw2 if sumw2 > 0 else 0.
                    row += [f"{sumw:.1f}", f"{n_eff:.0f}"]
            rows.append(row)

    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    # Same thing on screen, so it is readable without opening the file
    widths = [max(len(str(r[i])) for r in [header] + rows)
              for i in range(len(header))]
    for r in [header] + rows:
        print("  " + "  ".join(str(c).ljust(w) for c, w in zip(r, widths)))
    print(f"\nWrote {outfile}")


if __name__ == "__main__":
    main()
