#!/usr/bin/env python3
"""
Compare the radiation in the corridor between two jets for colour connected
and for unconnected jet pairs.

Makes one plot per cut and per observable:
  corridor_<obs>_<cut>   the three distributions, each normalised to unit area

The separation between connected and unconnected pairs is also printed as an
area-under-ROC number for each observable, so they can be ranked.

Usage
-----
    python scripts/plot_corridor_activity.py configs/plot_config.json \\
        output/hists/hists.json -o output/plots
"""
import os
from argparse import ArgumentParser
from warnings import filterwarnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep
import hist

from pepper import Config, HistCollection

filterwarnings("ignore", "kwarg `label` in function")
filterwarnings("ignore", "invalid value encountered in divide")
matplotlib.use("Agg")
plt.set_loglevel("error")
plt.style.use(mplhep.style.CMS)

# (observable suffix, x axis label)
OBSERVABLES = [
    ("ptdens", "Corridor $\\Sigma p_{T}$ per unit area [GeV]"),
    ("mass", "Corridor invariant mass [GeV]"),
    ("combined", "Corridor $(\\Sigma p_{T} + m)$ per unit area [GeV]"),
]
# (group suffix, legend label, colour). The first is the signal for the ROC,
# every later one is compared against it.
GROUPS = [
    ("connected", "Colour connected (same W)", "tab:blue"),
    ("crossW", "Not connected ($W^{+}$ with $W^{-}$)", "tab:red"),
    ("bb", "$b\\bar{b}$ pair", "tab:green"),
]


def integrate_categories(h):
    """Sum over every category axis, keeping only the nominal systematic."""
    if "sys" in h.axes.name:
        h = h[{"sys": "nominal"}]
    for ax in list(h.axes):
        if isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory)):
            h = h.integrate(ax.name, list(ax))
    return h


def roc_auc(signal, background):
    """How well a cut on this observable separates the two, as a single
    number. 0.5 means no separation, 1 means perfect. Printed only, so that
    the three observables can be ranked without needing another plot."""
    sig = signal / signal.sum()
    bkg = background / background.sum()
    # Fraction of each kept by a cut placed at every bin edge, from the top
    kept_sig = np.concatenate([[0.], np.cumsum(sig[::-1])])
    kept_bkg = np.concatenate([[0.], np.cumsum(bkg[::-1])])
    return np.trapezoid(kept_sig, kept_bkg) if hasattr(np, "trapezoid") \
        else np.trapz(kept_sig, kept_bkg)


def cms_label(ax, config):
    if "year" not in config:
        return
    label_kwargs = {}
    cmslabel = config["cmslabel"] if "cmslabel" in config else None
    if cmslabel is not None and cmslabel.strip().lower() != "simulation":
        label_kwargs["label"] = cmslabel
    mplhep.cms.label(
        ax=ax, data=False, year=config["year"],
        lumi=f"{config['luminosity']:.1f}" if "luminosity" in config else None,
        com=config["com_energy"] if "com_energy" in config else 13,
        **label_kwargs)


def plot_distributions(entries, axis, xlabel, outfile, exts, config, log):
    """Overlay the normalised distributions as step lines."""
    fig, ax = plt.subplots()
    edges = axis.edges
    widths = np.diff(edges)
    for label, color, counts, variances in entries:
        area = np.sum(counts * widths)
        if area == 0:
            continue
        mplhep.histplot(counts / area, edges, yerr=np.sqrt(variances) / area,
                        histtype="step", edges=False, linewidth=1.5,
                        color=color, label=label, ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalised pairs")
    ax.set_xlim(edges[0], edges[-1])
    if log:
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.legend(fontsize=16)
    cms_label(ax, config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def main():
    parser = ArgumentParser(
        description="Compare corridor radiation for connected and "
                    "unconnected jet pairs")
    parser.add_argument("plot_config", help="Plotting config, used for the "
                        "CMS label only")
    parser.add_argument("histfile", help="hists.json written by runproc")
    parser.add_argument("-o", "--outdir", help="Output directory. Defaults to "
                        "the directory containing histfile")
    parser.add_argument("--ext", choices=["pdf", "svg", "png"], default=None,
                        action="append", help="Output file format")
    parser.add_argument("--log", action="store_true",
                        help="Logarithmic y axis on the distribution plots")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]

    config = Config(args.plot_config)
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    wanted = {f"corridor_{obs}_{grp}"
              for obs, _ in OBSERVABLES for grp, _, _ in GROUPS}
    cuts = sorted({k[0] for k in hists.keys() if k[1] in wanted})
    if len(cuts) == 0:
        raise SystemExit(f"No corridor histograms found in {args.histfile}")

    ranking = {}
    for cutname in cuts:
        for obs, xlabel in OBSERVABLES:
            entries, axis = [], None
            for group, label, color in GROUPS:
                key = (cutname, f"corridor_{obs}_{group}")
                if key not in hists.keys():
                    continue
                h = integrate_categories(hists.load(key))
                dense = [a for a in h.axes
                         if not isinstance(a, (hist.axis.StrCategory,
                                               hist.axis.IntCategory))]
                if len(dense) != 1:
                    continue
                if axis is None:
                    axis = dense[0]
                elif dense[0].edges.shape != axis.edges.shape:
                    print(f"Skipping {key}: binning differs")
                    continue
                entries.append((label, color, h.values(), h.variances()))

            if len(entries) < 2:
                continue

            cut_outdir = os.path.join(outdir, cutname)
            os.makedirs(cut_outdir, exist_ok=True)

            plot_distributions(
                entries, axis, xlabel,
                os.path.join(cut_outdir, f"corridor_{obs}_{cutname}"),
                args.ext, config, args.log)

            # Separation of the first group from each of the others, as
            # numbers only - no plot
            sig_label, _, sig_counts, _ = entries[0]
            print(f"[{cutname}] corridor {obs}")
            print(f"  mean, {sig_label:<34}: "
                  f"{np.average(axis.centers, weights=sig_counts):.3f}")
            for label, color, counts, variances in entries[1:]:
                auc = roc_auc(sig_counts, counts)
                print(f"  mean, {label:<34}: "
                      f"{np.average(axis.centers, weights=counts):.3f}")
                print(f"  AUC  vs {label:<34}: {auc:.4f}")
                ranking.setdefault((cutname, label), []).append((auc, obs))

    for (cutname, label), aucs in ranking.items():
        print(f"\nMost effective observable, {cutname}, against {label}:")
        for auc, obs in sorted(aucs, reverse=True):
            print(f"    {obs:<10} AUC = {auc:.4f}")

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
