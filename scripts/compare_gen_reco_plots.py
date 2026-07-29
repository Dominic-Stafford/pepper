#!/usr/bin/env python3
"""
Overlay gen-level and reco-level pull angle distributions on the same axes.

For every histogram pair (`X`, `reco_X`) found in the pepper histogram
output, this draws both as step lines, each normalised to unit area, with a
reco/gen ratio panel underneath.

Usage
-----
    python scripts/compare_gen_reco_plots.py configs/plot_config.json \\
        output/hists/hists.json -o output/comparison

The plot config is only used for cosmetics (year, luminosity, com_energy,
cmslabel) and for the list of datasets to include; the stacking/background
machinery of make_plots.py is not used here.
"""
import os
import json
from argparse import ArgumentParser
from warnings import filterwarnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep
import hist
import tqdm

from pepper import Config, HistCollection

filterwarnings("ignore", "invalid value encountered in divide")
filterwarnings("ignore", "divide by zero encountered in divide")
filterwarnings("ignore", "kwarg `label` in function")

matplotlib.use("Agg")
plt.set_loglevel("error")
plt.style.use(mplhep.style.CMS)


# Curves to draw for each gen-level variable `X`. The prefix is prepended to
# the variable name to find the corresponding histogram; a curve is simply
# skipped if its histogram is not in the collection.
# (prefix, legend label, colour, line style, prefix of the curve to divide by
#  in the ratio panel - None means this curve is itself a reference)
CURVES = [
    ("",        "Gen level",          "tab:blue",   "-",  None),
    ("reco_",   "Reco level",         "tab:red",    "-",  ""),
    ("puppi_",  "Reco level (PUPPI)", "tab:green",  "-",  ""),
    ("genwta_", "Gen level (WTA)",    "tab:cyan",   "--", None),
    ("wta_",    "Reco level (WTA)",   "tab:purple", "-",  "genwta_"),
]

# Any variable starting with one of these is a counterpart, not a gen variable
COUNTERPART_PREFIXES = tuple(c[0] for c in CURVES if c[0])

# Curve that must always be present for a plot to be made
BASE_PREFIX = ""


def integrate_categories(h):
    """Sum a histogram over every category axis (dataset, sys, ...) so that
    only the dense axis is left.

    For the `sys` axis only the nominal is kept, if it is present.
    """
    if "sys" in h.axes.name:
        h = h[{"sys": "nominal"}]
    for ax in list(h.axes):
        if isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory)):
            h = h.integrate(ax.name, list(ax))
    return h


def get_dense_axis(h):
    """Return the single dense axis of `h`, or None if there isn't exactly
    one."""
    dense = [
        ax for ax in h.axes
        if not isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory))
    ]
    if len(dense) != 1:
        return None
    return dense[0]


def normalise(values, variances, widths):
    """Normalise to unit area. Returns (values, errors), both nan-safe."""
    area = np.sum(values * widths)
    if area == 0:
        return np.zeros_like(values), np.zeros_like(values)
    return values / area, np.sqrt(variances) / area


def plot_comparison(curves, dense_axis, outfile, exts, config, log=False):
    """Draw an overlay of all `curves` plus a ratio panel against the
    reference curve.

    `curves` is a list of (prefix, label, color, hist) tuples.
    """
    fig, (ax1, ax2) = plt.subplots(
        nrows=2, sharex=True, gridspec_kw={"height_ratios": [2, 1]})

    edges = dense_axis.edges
    widths = np.diff(edges)

    normalised = {}
    for prefix, label, color, style, ref_prefix, h in curves:
        vals, errs = normalise(h.values(), h.variances(), widths)
        normalised[prefix] = (vals, errs)
        # --- upper panel: step lines ---
        mplhep.histplot(
            vals, edges, yerr=errs,
            histtype="step", edges=False, linewidth=1.5, linestyle=style,
            color=color, label=label, ax=ax1)

    # --- lower panel: each curve divided by its own reference ---
    all_ratios = []
    for prefix, label, color, style, ref_prefix, h in curves:
        if ref_prefix is None or ref_prefix not in normalised:
            continue
        ref_vals, ref_errs = normalised[ref_prefix]
        vals, errs = normalised[prefix]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref_vals > 0, vals / ref_vals, np.nan)
            # Error propagation treating the two as independent
            rel_ref = np.where(ref_vals > 0, ref_errs / ref_vals, 0.)
            rel_this = np.where(vals > 0, errs / vals, 0.)
            ratio_errs = np.abs(ratio) * np.sqrt(rel_ref**2 + rel_this**2)
        all_ratios.append(ratio)
        mplhep.histplot(
            np.nan_to_num(ratio, nan=0.), edges,
            yerr=np.nan_to_num(ratio_errs),
            histtype="errorbar", color=color, markersize=8, ax=ax2)

    ratio = np.concatenate(all_ratios) if all_ratios else np.array([1.])

    ax2.hlines(1, edges[0], edges[-1], color="black", linestyle="dashed",
               linewidth=1., zorder=-90)

    # --- axis cosmetics ---
    ax1.set_xlabel("")
    ax1.set_ylabel("a.u.")
    ax1.set_ylim(bottom=0)
    if log:
        ax1.set_yscale("log")
        ax1.autoscale("y")
    ax1.legend()

    ax2.set_ylabel("Reco / gen")
    ax2.set_xlabel(dense_axis.label)
    ax2.set_xlim(edges[0], edges[-1])
    finite = ratio[np.isfinite(ratio)]
    if len(finite) > 0:
        spread = np.max(np.abs(finite - 1.))
        spread = min(max(spread, 0.05), 1.0)
    else:
        spread = 0.5
    ax2.set_ylim(1. - spread * 1.2, 1. + spread * 1.2)

    if "year" in config:
        # With data=False mplhep writes "Simulation" itself, so passing a
        # cmslabel of "Simulation" from the config would print it twice.
        # Only forward the config label if it says something else.
        label_kwargs = {}
        cmslabel = config["cmslabel"] if "cmslabel" in config else None
        if cmslabel is not None and cmslabel.strip().lower() != "simulation":
            label_kwargs["label"] = cmslabel
        mplhep.cms.label(
            ax=ax1,
            data=False,
            year=config["year"],
            lumi=f"{config['luminosity']:.1f}"
                 if "luminosity" in config else None,
            com=config["com_energy"] if "com_energy" in config else 13,
            **label_kwargs)

    fig.align_ylabels()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0)

    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def main():
    parser = ArgumentParser(
        description="Overlay gen-level and reco-level distributions")
    parser.add_argument(
        "plot_config", help="Path to a plotting config file (used for the "
        "CMS label, year and luminosity only)")
    parser.add_argument(
        "histfile", help="JSON file containing histogram info, i.e. the "
        "hists.json written by runproc")
    parser.add_argument(
        "-o", "--outdir", help="Output directory. Defaults to the directory "
        "containing histfile")
    parser.add_argument(
        "--log", action="store_true", help="Make logarithmic plots")
    parser.add_argument(
        "--ext", choices=["pdf", "svg", "png"], default=None, action="append",
        help="Output file format, can be given multiple times")
    parser.add_argument(
        "-c", "--cut", default=None, help="Only plot this cut")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]

    config = Config(args.plot_config)

    if not args.histfile.endswith(".json"):
        raise SystemExit(
            "This script needs the hists.json file, not a single histogram")
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    # Collect the gen-level variables that have at least one counterpart
    available = set(hists.keys())
    todo = []
    for key in sorted(available):
        cutname, variable = key[0], key[1]
        if variable.startswith(COUNTERPART_PREFIXES):
            continue
        keys = {}
        for curve in CURVES:
            cand = (cutname, curve[0] + variable) + tuple(key[2:])
            if cand in available:
                keys[curve[0]] = cand
        # Need the gen curve plus at least one thing to compare it against
        if BASE_PREFIX in keys and len(keys) > 1:
            todo.append((cutname, variable, keys))

    if len(todo) == 0:
        raise SystemExit(
            f"Found no variable with a matching "
            f"{list(COUNTERPART_PREFIXES)} counterpart in {args.histfile}")

    for cutname, variable, keys in tqdm.tqdm(todo):
        if args.cut is not None and cutname != args.cut:
            continue

        curves = []
        dense_axis = None
        for prefix, label, color, style, ref_prefix in CURVES:
            if prefix not in keys:
                continue
            h = integrate_categories(hists.load(keys[prefix]))
            axis = get_dense_axis(h)
            if axis is None:
                print(f"Skipping {prefix}{variable} ({cutname}): not exactly "
                      f"one dense axis")
                continue
            if dense_axis is None:
                dense_axis = axis
            elif axis.edges.shape != dense_axis.edges.shape:
                print(f"Skipping {prefix}{variable} ({cutname}): binning "
                      f"differs from the gen-level histogram")
                continue
            curves.append((prefix, label, color, style, ref_prefix, h))

        if len(curves) < 2:
            continue

        cut_outdir = os.path.join(outdir, cutname)
        os.makedirs(cut_outdir, exist_ok=True)
        outfile = os.path.join(cut_outdir, f"genreco_{variable}_{cutname}")

        plot_comparison(curves, dense_axis, outfile, args.ext,
                        config, log=args.log)

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
