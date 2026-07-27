#!/usr/bin/env python3
"""
Overlay gen-level and reco-level pull angle distributions on the same axes.

For every histogram pair (`X`, `reco_X`) found in the pepper histogram
output, this draws both as step lines, each normalised to unit area, with a
reco/gen ratio panel underneath.
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


# Prefix that marks the reco-level counterpart of a gen-level variable
RECO_PREFIX = "reco_"

GEN_COLOR = "tab:blue"
RECO_COLOR = "tab:red"


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


def plot_comparison(h_gen, h_reco, dense_axis, outfile, exts, config,
                    log=False):
    """Draw the gen/reco overlay plus a reco/gen ratio panel."""
    fig, (ax1, ax2) = plt.subplots(
        nrows=2, sharex=True, gridspec_kw={"height_ratios": [2, 1]})

    edges = dense_axis.edges
    widths = np.diff(edges)

    gen_vals, gen_errs = normalise(
        h_gen.values(), h_gen.variances(), widths)
    reco_vals, reco_errs = normalise(
        h_reco.values(), h_reco.variances(), widths)

    # --- upper panel: both distributions as step lines ---
    mplhep.histplot(
        gen_vals, edges, yerr=gen_errs,
        histtype="step", edges=False, linewidth=1.5,
        color=GEN_COLOR, label="Gen level", ax=ax1)
    mplhep.histplot(
        reco_vals, edges, yerr=reco_errs,
        histtype="step", edges=False, linewidth=1.5,
        color=RECO_COLOR, label="Reco level", ax=ax1)

    # --- lower panel: reco / gen ---
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(gen_vals > 0, reco_vals / gen_vals, np.nan)
        # Error propagation treating the two as independent
        rel_gen = np.where(gen_vals > 0, gen_errs / gen_vals, 0.)
        rel_reco = np.where(reco_vals > 0, reco_errs / reco_vals, 0.)
        ratio_errs = np.abs(ratio) * np.sqrt(rel_gen**2 + rel_reco**2)

    mplhep.histplot(
        np.nan_to_num(ratio, nan=0.), edges, yerr=np.nan_to_num(ratio_errs),
        histtype="errorbar", color="black", markersize=8, ax=ax2)

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

    ax2.set_ylabel("Reco / Gen")
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

    # Build the list of (cut, gen_variable) pairs that have a reco counterpart
    available = set(hists.keys())
    pairs = []
    for key in sorted(available):
        cutname, variable = key[0], key[1]
        if variable.startswith(RECO_PREFIX):
            continue
        reco_key = (cutname, RECO_PREFIX + variable) + tuple(key[2:])
        if reco_key in available:
            pairs.append((key, reco_key))

    if len(pairs) == 0:
        raise SystemExit(
            f"Found no variable with a matching '{RECO_PREFIX}' counterpart "
            f"in {args.histfile}")

    for gen_key, reco_key in tqdm.tqdm(pairs):
        cutname, variable = gen_key[0], gen_key[1]
        if args.cut is not None and cutname != args.cut:
            continue

        h_gen = integrate_categories(hists.load(gen_key))
        h_reco = integrate_categories(hists.load(reco_key))

        dense_axis = get_dense_axis(h_gen)
        if dense_axis is None:
            print(f"Skipping {variable} ({cutname}): not exactly one dense "
                  f"axis")
            continue
        if get_dense_axis(h_reco).edges.shape != dense_axis.edges.shape:
            print(f"Skipping {variable} ({cutname}): gen and reco binning "
                  f"differ")
            continue

        cut_outdir = os.path.join(outdir, cutname)
        os.makedirs(cut_outdir, exist_ok=True)
        outfile = os.path.join(cut_outdir, f"genreco_{variable}_{cutname}")

        plot_comparison(h_gen, h_reco, dense_axis, outfile, args.ext,
                        config, log=args.log)

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()