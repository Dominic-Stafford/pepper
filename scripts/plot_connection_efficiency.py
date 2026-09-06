#!/usr/bin/env python3
"""
How often does a given method pick out the true colour connected partner
among the four W decay jets?

Makes two plots per cut:
  pull_angle_connection_efficiency_<cut>
      distribution of how many of the four jets the pull angle assigned
      correctly, from 0 to 4, with binomial expectations overlaid
  mass_connection_efficiency_<cut>
      success and failure rate of the two W mass pairing scores, which are
      all-or-nothing and so have no such distribution

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

# (variable, legend label, colour, has a spread over 0 to 4)
# The mass pairing picks one of three ways to split the four jets into two
# pairs, so it is either all right or all wrong - it has no spread, and it
# is left out of the distribution plot.
LEVELS = [
    ("n_correct_connections", "Pull angle", "tab:blue", True),
    ("n_correct_connections_mass_sum", "W mass, sum", "tab:orange", False),
    ("n_correct_connections_mass_max", "W mass, max", "tab:green", False),
]
N_JETS = 4
RANDOM = 1 / (N_JETS - 1)

# The eight jet problem: three connected pairs picked out of eight jets, with
# two left over, and one of the three assigned to the Higgs.
# C(8,2) * 15 * 3 = 1260 options, so guessing gets it right 1/1260 of the time.
EIGHT_OPTIONS = 1260
EIGHT_METHODS = [
    ("eight_pairing_correct_jetmass", "Jet mass", "tab:blue"),
    ("eight_pairing_correct_capsule_r02", "Capsule mass, R = 0.2", "tab:orange"),
    ("eight_pairing_correct_capsule_r03", "Capsule mass, R = 0.3", "tab:green"),
]


def load_binary(hists, cutname, variable, dataset=None):
    """Fraction of entries in the upper bin of a two bin histogram, with its
    uncertainty from the effective number of entries."""
    key = (cutname, variable)
    if key not in hists.keys():
        return None
    h = integrate_categories(hists.load(key), dataset)
    if h is None:
        return None
    counts, variances = h.values(), h.variances()
    total = counts.sum()
    if total <= 0 or len(counts) != 2:
        return None
    frac = counts[1] / total
    n_eff = total**2 / variances.sum()
    return frac, np.sqrt(max(frac * (1 - frac), 0) / n_eff)


def integrate_categories(h, dataset=None):
    """Sum over every category axis, keeping only the nominal systematic.
    If `dataset` is given, keep only that dataset instead of summing them.
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


def binomial_pmf(p):
    """Probability of exactly k of N_JETS right, if each jet were an
    independent trial with success probability p."""
    return np.array([comb(N_JETS, k) * p**k * (1 - p)**(N_JETS - k)
                     for k in range(N_JETS + 1)])


def get_efficiency(counts, variances):
    """Per-jet efficiency and its uncertainty.

    Each event gives a number n of correctly assigned jets out of N_JETS, so
    the efficiency is mean(n) / N_JETS and the uncertainty is the standard
    error on that mean. Using the observed spread of n rather than a binomial
    formula matters, because the assignments within an event are correlated.

    Entries are weighted, so counts.sum() is a sum of weights and not a number
    of events - the effective sample size (sum w)^2 / sum w^2 is used instead.
    """
    n_values = np.arange(len(counts))
    sum_w = counts.sum()
    fractions = counts / sum_w
    mean_n = np.dot(fractions, n_values)
    var_n = np.dot(fractions, n_values**2) - mean_n**2
    n_eff = sum_w**2 / variances.sum()
    return mean_n / N_JETS, np.sqrt(var_n / n_eff) / N_JETS, \
        fractions, var_n, n_eff


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


def plot_ndist(entries, outfile, exts, config, extra_p):
    """Distribution of the number of correctly assigned jets per event, with
    binomial expectations overlaid."""
    fig, ax = plt.subplots()
    n_values = np.arange(N_JETS + 1)
    highest = 0.

    for label, color, fractions, eff, eff_err in entries:
        ax.bar(n_values, fractions, width=0.65, color=color, alpha=0.75,
               edgecolor=color, linewidth=1.5,
               label=f"{label}: {eff:.1%} $\\pm$ {eff_err:.1%}")
        highest = max(highest, fractions.max())

    curves = [(RANDOM, f"Binomial, p = {RANDOM:.1%} (random)", "tab:gray", "--")]
    for label, color, fractions, eff, eff_err in entries:
        curves.append((eff, f"Binomial, p = {eff:.1%} (uncorrelated)",
                       "black", ":"))
    for p in extra_p:
        curves.append((p, f"Binomial, p = {p:.1%}", "tab:red", "-."))

    for p, label, color, linestyle in curves:
        pmf = binomial_pmf(p)
        ax.plot(n_values, pmf, linestyle=linestyle, color=color, linewidth=2,
                marker="o", markersize=8, label=label)
        highest = max(highest, pmf.max())

    ax.set_xlabel("Correctly assigned jets per event")
    ax.set_ylabel("Fraction of events")
    ax.set_xticks(n_values)
    ax.set_ylim(0, highest * 1.45)
    ax.legend(loc="upper left", fontsize=15)
    cms_label(ax, config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def plot_pairing(results, outfile, exts, config):
    """Success and failure rate of each all-or-nothing pairing method.

    `results` is a list of (label, fraction correct, uncertainty).
    """
    fig, ax = plt.subplots()

    labels = [r[0] for r in results]
    correct = np.array([r[1] for r in results])
    errors = np.array([r[2] for r in results])
    x = np.arange(len(results))
    width = 0.34

    for offset, values, color, name in [
            (-width / 2, correct, "tab:green", "Correct pairing"),
            (+width / 2, 1 - correct, "tab:red", "Wrong pairing")]:
        ax.bar(x + offset, values, width=width * 0.92, yerr=errors,
               color=color, alpha=0.8, edgecolor=color, linewidth=1.5,
               label=name, error_kw={"ecolor": "black", "capsize": 5})
        for xi, v in zip(x + offset, values):
            ax.text(xi, v + 0.03, f"{v:.1%}", ha="center", va="bottom",
                    fontsize=17)

    #ax.axhline(RANDOM, color="black", linestyle="--", linewidth=1.5)
    ax.text(len(results) - 0.45, RANDOM + 0.015, f"random = {RANDOM:.1%}",
            fontsize=15, va="bottom", ha="right")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of events")
    ax.set_ylim(0, 1.22)
    ax.legend(loc="upper center", ncol=2, fontsize=17)
    cms_label(ax, config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def plot_eight(results, outfile, exts, config, level):
    """Success rate of the eight jet assignment, one bar per method.

    This is a much harder problem than the four jet one: 1260 combinations
    rather than 3, so the baseline sits at 0.08 percent rather than 33.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    labels = [r[0] for r in results]
    colors = [r[1] for r in results]
    fracs = np.array([r[2] for r in results])
    errs = np.array([r[3] for r in results])
    x = np.arange(len(results))

    ax.bar(x, fracs, yerr=errs, width=0.6, color=colors, alpha=0.8,
           edgecolor=colors, linewidth=1.5,
           error_kw={"ecolor": "black", "capsize": 5})
    for xi, f, e in zip(x, fracs, errs):
        ax.text(xi, f + e + max(fracs) * 0.03, f"{f:.1%}", ha="center",
                va="bottom", fontsize=17)
    ax.axhline(1 / EIGHT_OPTIONS, color="black", linestyle="--", linewidth=1.5)
    ax.text(len(results) - 0.5, 1 / EIGHT_OPTIONS, 
            f"  random = {1/EIGHT_OPTIONS:.2%}", fontsize=14, va="bottom",
            ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=15)
    ax.set_ylabel("Events with all three pairs correct", fontsize=16)
    ax.set_ylim(0, max(fracs) * 1.35 if max(fracs) > 0 else 1.)
    cms_label(ax, config)
    fig.suptitle(f"Three connected pairs from eight jets, {level} level",
                 fontsize=16, y=0.995)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
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
                        metavar="P", help="Overlay an extra binomial "
                        "expectation for this per-jet probability")
    parser.add_argument("--split-datasets", action="store_true",
                        help="Make a separate set of plots for every dataset "
                        "instead of summing them, in a subdirectory per "
                        "dataset")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]

    config = Config(args.plot_config)
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    variables = [v for v, _, _, _ in LEVELS]
    cuts = sorted({k[0] for k in hists.keys() if k[1] in variables})
    if len(cuts) == 0:
        raise SystemExit(f"None of {variables} found in {args.histfile}")

    # None means sum over all datasets
    datasets = [None]
    if args.split_datasets:
        # Take the union over the histograms this script uses: some of them
        # are filled for only one sample, so one histogram alone would not
        # list every dataset.
        wanted = set(variables) | {p + v for p in ("", "reco_")
                                   for v, _, _ in EIGHT_METHODS}
        found = []
        for key in hists.keys():
            if key[1] not in wanted:
                continue
            axes = hists.load(key).axes
            if "dataset" not in axes.name:
                continue
            found += [d for d in axes["dataset"] if d not in found]
        if found:
            datasets = found
            print(f"Making a set of plots for each of: {datasets}")

    for cutname, dataset in ((c, d) for c in cuts for d in datasets):
        dist_entries = []      # pull angle, for the 0 to 4 distribution plot
        pairing_results = []   # W mass, for the success/failure plot

        for variable, label, color, has_spread in LEVELS:
            key = (cutname, variable)
            if key not in hists.keys():
                continue
            h = integrate_categories(hists.load(key), dataset)
            if h is None:
                continue
            counts, variances = h.values(), h.variances()
            if len(counts) != N_JETS + 1:
                print(f"Skipping {variable} ({cutname}): expected "
                      f"{N_JETS+1} bins, got {len(counts)}")
                continue
            eff, eff_err, fractions, var_n, n_eff = \
                get_efficiency(counts, variances)

            tag = cutname if dataset is None else f"{cutname}, {dataset}"
            print(f"[{tag}] {label}")
            print(f"  effective n events    : {n_eff:.0f}")
            print(f"  per-jet efficiency    : {eff:.4f} +- {eff_err:.4f}"
                  f"   (random = {RANDOM:.4f})")
            if has_spread:
                dist_entries.append((label, color, fractions, eff, eff_err))
                var_binom = N_JETS * eff * (1 - eff)
                print(f"  fraction per n correct: "
                      f"{np.array2string(fractions, precision=4)}")
                print(f"  if uncorrelated       : "
                      f"{np.array2string(binomial_pmf(eff), precision=4)}")
                print(f"  var(n) obs / binomial : "
                      f"{var_n:.4f} / {var_binom:.4f} = "
                      f"{var_n/var_binom:.3f}")
            else:
                pairing_results.append((label, fractions[-1], eff_err))
                print(f"  all-or-nothing method, correct pairing in "
                      f"{fractions[-1]:.1%} of events")

        # --- the eight jet assignment, jet mass against capsule mass ---
        for level, pre in [("gen", ""), ("reco", "reco_")]:
            eight = []
            for variable, label, color in EIGHT_METHODS:
                got = load_binary(hists, cutname, pre + variable, dataset)
                if got is None:
                    continue
                frac, err = got
                eight.append((label, color, frac, err))
                tag = cutname if dataset is None else f"{cutname}, {dataset}"
                print(f"[{tag}] {level} eight jet, {label:<22}: "
                      f"{frac:.2%} +- {err:.2%}   "
                      f"(random {1/EIGHT_OPTIONS:.2%})")
            if len(eight) > 0:
                out = os.path.join(outdir, cutname)
                if dataset is not None:
                    out = os.path.join(out, dataset)
                os.makedirs(out, exist_ok=True)
                plot_eight(eight,
                           os.path.join(out, f"{pre}eight_pairing_{cutname}"),
                           args.ext, config, level)

        cut_outdir = os.path.join(outdir, cutname)
        if dataset is not None:
            cut_outdir = os.path.join(cut_outdir, dataset)
        if len(dist_entries) > 0 or len(pairing_results) > 0:
            os.makedirs(cut_outdir, exist_ok=True)

        if len(dist_entries) > 0:
            plot_ndist(dist_entries,
                       os.path.join(
                           cut_outdir,
                           f"pull_angle_connection_efficiency_{cutname}"),
                       args.ext, config, args.binomial or [])
        if len(pairing_results) > 0:
            plot_pairing(pairing_results,
                         os.path.join(
                             cut_outdir,
                             f"mass_connection_efficiency_{cutname}"),
                         args.ext, config)

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
