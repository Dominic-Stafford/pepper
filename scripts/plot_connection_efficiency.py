#!/usr/bin/env python3
"""
Plot how often the smallest pull angle picks out the true colour connected
partner among the four W decay jets.

Reads the `n_correct_connections` histogram, which counts per event how many
of the four jets were assigned their true partner (so 0 to 4).

Usage
-----
    python scripts/plot_connection_efficiency.py configs/plot_config.json \\
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
matplotlib.use("Agg")
plt.set_loglevel("error")
plt.style.use(mplhep.style.CMS)

VARIABLE = "n_correct_connections"
N_JETS = 4


def integrate_categories(h):
    """Sum over every category axis, keeping only the nominal systematic."""
    if "sys" in h.axes.name:
        h = h[{"sys": "nominal"}]
    for ax in list(h.axes):
        if isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory)):
            h = h.integrate(ax.name, list(ax))
    return h


def plot_efficiency(counts, variances, outfile, exts, config):
    """Bar chart of the fraction of events with n correctly assigned jets."""
    n_events = counts.sum()
    fractions = counts / n_events
    n_values = np.arange(len(counts))

    # Per-jet efficiency: total correct assignments over total assignments made
    n_correct_total = np.dot(counts, n_values)
    efficiency = n_correct_total / (N_JETS * n_events)
    # The entries are weighted, so counts.sum() is a sum of weights and not a
    # number of events. Use the effective sample size for the uncertainty.
    n_eff = n_events**2 / variances.sum()
    eff_err = np.sqrt(efficiency * (1 - efficiency) / (N_JETS * n_eff))
    all_right = fractions[-1]

    fig, ax = plt.subplots()

    ax.bar(n_values, fractions, width=0.7, color="tab:blue", alpha=0.75,
           edgecolor="tab:blue", linewidth=1.5)

    # Print the fraction on top of each bar
    for n, frac in zip(n_values, fractions):
        if frac > 0:
            ax.text(n, frac + 0.015, f"{frac:.1%}", ha="center", va="bottom",
                    fontsize=15)

    ax.set_xlabel("Correctly assigned jets per event")
    ax.set_ylabel("Fraction of events")
    ax.set_xticks(n_values)
    ax.set_ylim(0, max(fractions) * 1.30)

    summary = (f"Per-jet efficiency: {efficiency:.1%} $\\pm$ {eff_err:.1%}\n"
               f"All {N_JETS} correct: {all_right:.1%}\n"
               f"Random guess: {1/(N_JETS-1):.1%}")
    ax.text(0.03, 0.97, summary, transform=ax.transAxes, va="top", ha="left",
            fontsize=17)

    if "year" in config:
        label_kwargs = {}
        cmslabel = config["cmslabel"] if "cmslabel" in config else None
        if cmslabel is not None and cmslabel.strip().lower() != "simulation":
            label_kwargs["label"] = cmslabel
        mplhep.cms.label(
            ax=ax, data=False, year=config["year"],
            lumi=f"{config['luminosity']:.1f}" if "luminosity" in config
                 else None,
            com=config["com_energy"] if "com_energy" in config else 13,
            **label_kwargs)

    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)

    return efficiency, eff_err, all_right


def main():
    parser = ArgumentParser(
        description="Plot the colour connection assignment efficiency")
    parser.add_argument("plot_config", help="Plotting config, used for the "
                        "CMS label only")
    parser.add_argument("histfile", help="hists.json written by runproc")
    parser.add_argument("-o", "--outdir", help="Output directory. Defaults to "
                        "the directory containing histfile")
    parser.add_argument("--ext", choices=["pdf", "svg", "png"], default=None,
                        action="append", help="Output file format")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]

    config = Config(args.plot_config)
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    keys = [k for k in hists.keys() if k[1] == VARIABLE]
    if len(keys) == 0:
        raise SystemExit(f"No '{VARIABLE}' histogram in {args.histfile}")

    for key in sorted(keys):
        cutname = key[0]
        h = integrate_categories(hists.load(key))
        counts = h.values()
        variances = h.variances()
        if len(counts) != N_JETS + 1:
            print(f"Skipping {cutname}: expected {N_JETS+1} bins, "
                  f"got {len(counts)}")
            continue

        cut_outdir = os.path.join(outdir, cutname)
        os.makedirs(cut_outdir, exist_ok=True)
        outfile = os.path.join(cut_outdir, f"{VARIABLE}_{cutname}")

        eff, eff_err, all_right = plot_efficiency(
            counts, variances, outfile, args.ext, config)

        print(f"[{cutname}]")
        print(f"  sum of weights         : {counts.sum():.0f}")
        print(f"  effective n events     : "
              f"{counts.sum()**2/variances.sum():.0f}")
        print(f"  counts per n correct   : {counts.astype(int)}")
        print(f"  per-jet efficiency     : {eff:.4f} +- {eff_err:.4f}")
        print(f"  all {N_JETS} correct          : {all_right:.4f}")
        print(f"  random-guess baseline  : {1/(N_JETS-1):.4f}")

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
