#!/usr/bin/env python3
"""
Plot how often the smallest pull angle picks out the true colour connected
partner among the four W decay jets, at gen and at reco level.

Reads the `n_correct_connections` and `reco_n_correct_connections` histograms,
which count per event how many of the four jets were assigned their true
partner (so 0 to 4).

Usage
-----
    python scripts/plot_connection_efficiency.py configs/plot_config.json \\
        output/hists/hists.json -o output/plots
"""
import os
from math import comb
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

# (variable name, legend label, colour). The first one is the method the
# binomial expectation is drawn for.
LEVELS = [
    ("n_correct_connections", "Pull angle", "tab:blue"),
    ("n_correct_connections_mass_sum", "W mass, sum", "tab:orange"),
    ("n_correct_connections_mass_max", "W mass, max", "tab:green"),
]
N_JETS = 4


def integrate_categories(h):
    """Sum over every category axis, keeping only the nominal systematic."""
    if "sys" in h.axes.name:
        h = h[{"sys": "nominal"}]
    for ax in list(h.axes):
        if isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory)):
            h = h.integrate(ax.name, list(ax))
    return h


def get_efficiency(counts, variances):
    """Per-jet assignment efficiency and its uncertainty.

    Each event contributes a number n of correctly assigned jets out of
    N_JETS, so the efficiency is the mean of n divided by N_JETS. The
    uncertainty is the standard error on that mean, using the spread of n
    across events. Doing it that way rather than assuming N_JETS independent
    binomial trials per event matters, because the assignments within an
    event are correlated - the jets pick from the same pool.

    The histogram entries are weighted, so counts.sum() is a sum of weights
    and not a number of events. The effective sample size (sum w)^2 / sum w^2
    is used instead.
    """
    n_values = np.arange(len(counts))
    sum_w = counts.sum()
    fractions = counts / sum_w

    mean_n = np.dot(fractions, n_values)
    var_n = np.dot(fractions, n_values**2) - mean_n**2
    n_eff = sum_w**2 / variances.sum()

    efficiency = mean_n / N_JETS
    eff_err = np.sqrt(var_n / n_eff) / N_JETS
    return efficiency, eff_err, fractions, n_eff


def binomial_pmf(p):
    """Probability of getting exactly k out of N_JETS right, if each jet were
    an independent trial with success probability p."""
    return np.array([comb(N_JETS, k) * p**k * (1 - p)**(N_JETS - k)
                     for k in range(N_JETS + 1)])


def plot_efficiency(entries, curves, outfile, exts, config):
    """Grouped bar chart of the fraction of events with n correct jets, with
    binomial expectations overlaid as lines.

    `curves` is a list of (p, label, colour, linestyle).
    """
    fig, ax = plt.subplots()

    n_bars = len(entries)
    width = 0.8 / n_bars
    n_values = np.arange(N_JETS + 1)
    highest = 0.

    for i, (label, color, fractions, efficiency, eff_err) in \
            enumerate(entries):
        offset = (i - (n_bars - 1) / 2) * width
        ax.bar(n_values + offset, fractions, width=width * 0.92,
               color=color, alpha=0.75, edgecolor=color, linewidth=1.5,
               label=f"{label}: {efficiency:.1%} $\\pm$ {eff_err:.1%}")
        highest = max(highest, fractions.max())

    # Binomial expectations, i.e. what four independent jets would give
    for p, label, color, linestyle in curves:
        pmf = binomial_pmf(p)
        ax.plot(n_values, pmf, linestyle=linestyle, color=color,
                linewidth=2, marker="o", markersize=9, label=label)
        highest = max(highest, pmf.max())

    ax.set_xlabel("Correctly assigned jets per event")
    ax.set_ylabel("Fraction of events")
    ax.set_xticks(n_values)
    ax.set_ylim(0, highest * 1.45)
    ax.legend(loc="upper left", fontsize=16)

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
    parser.add_argument("--binomial", type=float, default=None, action="append",
                        metavar="P", help="Overlay a binomial expectation for "
                        "this per-jet probability. Can be given several times. "
                        "The random-guess baseline and the measured efficiency "
                        "are always drawn.")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]

    config = Config(args.plot_config)
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    variables = [v for v, _, _ in LEVELS]
    cuts = sorted({k[0] for k in hists.keys() if k[1] in variables})
    if len(cuts) == 0:
        raise SystemExit(
            f"None of {variables} found in {args.histfile}")

    for cutname in cuts:
        entries = []
        for variable, label, color in LEVELS:
            key = (cutname, variable)
            if key not in hists.keys():
                continue
            counts = integrate_categories(hists.load(key))
            variances = counts.variances()
            counts = counts.values()
            if len(counts) != N_JETS + 1:
                print(f"Skipping {variable} ({cutname}): expected "
                      f"{N_JETS+1} bins, got {len(counts)}")
                continue
            eff, eff_err, fractions, n_eff = get_efficiency(counts, variances)
            entries.append((label, color, fractions, eff, eff_err))

            # Spread of the observed distribution against the binomial one at
            # the same per-jet efficiency. A ratio above 1 means the four
            # assignments in an event are positively correlated.
            n_values = np.arange(N_JETS + 1)
            var_obs = (np.dot(fractions, n_values**2)
                       - np.dot(fractions, n_values)**2)
            var_binom = N_JETS * eff * (1 - eff)

            print(f"[{cutname}] {variable}")
            print(f"  sum of weights        : {counts.sum():.0f}")
            print(f"  effective n events    : {n_eff:.0f}")
            print(f"  fraction per n correct: "
                  f"{np.array2string(fractions, precision=4)}")
            print(f"  binomial at this eff. : "
                  f"{np.array2string(binomial_pmf(eff), precision=4)}")
            print(f"  per-jet efficiency    : {eff:.4f} +- {eff_err:.4f}")
            print(f"  all {N_JETS} correct         : {fractions[-1]:.4f}")
            print(f"  random-guess baseline : {1/(N_JETS-1):.4f}")
            print(f"  var(n) observed       : {var_obs:.4f}")
            print(f"  var(n) binomial       : {var_binom:.4f}")
            print(f"  ratio (>1 = correlated): {var_obs/var_binom:.4f}")

        if len(entries) == 0:
            continue

        # Always show the random-guess baseline and the measured efficiency
        curves = [(1 / (N_JETS - 1),
                   f"Binomial, p = {1/(N_JETS-1):.1%} (random)",
                   "tab:gray", "--")]
        # Only for the first method - for the mass pairing the result is 0 or
        # 4 by construction, so a binomial comparison says nothing
        label, color, fractions, eff, eff_err = entries[0]
        curves.append((eff, f"Binomial, p = {eff:.1%} ({label})", "black", ":"))
        for p in (args.binomial or []):
            curves.append((p, f"Binomial, p = {p:.1%}", "tab:green", "-."))

        cut_outdir = os.path.join(outdir, cutname)
        os.makedirs(cut_outdir, exist_ok=True)
        outfile = os.path.join(cut_outdir, f"connection_efficiency_{cutname}")
        plot_efficiency(entries, curves, outfile, args.ext, config)

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
