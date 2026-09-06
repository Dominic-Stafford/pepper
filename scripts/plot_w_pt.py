#!/usr/bin/env python3
"""
Distribution of the W transverse momentum, one curve for the W+ and one for
the W-.

The W four-momentum is taken at gen level as the sum of its two decay quarks,
so this is the true W pt and not a reconstructed one.

Makes one plot per cut:
  w_pt_<cut>   the two distributions, each normalised to unit area

With --split-datasets a separate plot is made for every sample, in a
subdirectory per dataset, so tt and ttH can be compared.

Usage
-----
    python scripts/plot_w_pt.py configs/plot_config.json \\
        output/hists/hists.json -o output/plots --split-datasets
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

# (variable, legend label, colour, include in the combined curve)
CURVES = [
    ("Wplus_pt", "$W^{+}$", "tab:blue", True),
    ("Wminus_pt", "$W^{-}$", "tab:red", True),
    ("H_pt", "$H$", "tab:green", False),
]
# Two dimensional pt against pull angle, summed over the two W charges. The
# same gen level W pt is used for every variant so they stay comparable.
# variant name -> (histogram prefix, axis label)
PULL_VARIANTS = {
    "gen":  ("",      "Gen pull angle"),
    "reco": ("reco_", "Reco pull angle"),
}

# W pt groups to draw one curve for. The histogram is filled with a finer
# binning than this, so these can be changed freely without rerunning the
# processor, as long as every boundary falls on one of the fine bin edges.
# None as the upper edge means everything above.
WPT_GROUPS = [(0, 60), (60, 120), (120, 200), (200, None)]


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


def load(hists, cutname, variable, dataset=None):
    """Return (dense axis, counts, variances), or None if not present."""
    key = (cutname, variable)
    if key not in hists.keys():
        return None
    h = integrate_categories(hists.load(key), dataset)
    if h is None:
        return None
    dense = [a for a in h.axes
             if not isinstance(a, (hist.axis.StrCategory,
                                   hist.axis.IntCategory))]
    if len(dense) != 1:
        return None
    return dense[0], h.values(), h.variances()


def load_2d(hists, cutname, variable, dataset=None):
    """Return (wpt axis, phi axis, values, variances) for a two dimensional
    histogram, or None if it is not there."""
    key = (cutname, variable)
    if key not in hists.keys():
        return None
    h = integrate_categories(hists.load(key), dataset)
    if h is None:
        return None
    if set(h.axes.name) != {"wpt", "phi"}:
        return None
    h = h.project("wpt", "phi")
    return h.axes["wpt"], h.axes["phi"], h.values(), h.variances()


def group_wpt(wpt_axis, values, variances):
    """Merge the fine W pt bins into the groups given by WPT_GROUPS.

    Returns a list of (label, counts, variances). A warning is printed if a
    group boundary does not line up with a fine bin edge, since the group
    would then not cover what it claims to.
    """
    edges = wpt_axis.edges
    grouped = []
    for lo, hi in WPT_GROUPS:
        top = edges[-1] if hi is None else hi
        for value in (lo, top):
            if not np.any(np.isclose(edges, value)):
                print(f"  WARNING: {value:g} GeV is not a bin edge, the "
                      f"group will be rounded to the nearest ones")
        keep = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= top + 1e-9)
        if not keep.any():
            continue
        label = (f"$>$ {lo:.0f} GeV" if hi is None
                 else f"{lo:.0f} to {hi:.0f} GeV")
        grouped.append((label, values[keep].sum(axis=0),
                        variances[keep].sum(axis=0)))
    return grouped


def plot_pull_vs_wpt(grouped, phi_axis, outfile, exts, config, xlabel):
    """One panel per W pt group, laid out as a grid so the four are easy to
    compare without eight lines on top of each other."""
    n = len(grouped)
    n_cols = 2 if n > 1 else 1
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, squeeze=False, sharex=True,
                             sharey=True, figsize=(6.2 * n_cols, 5.2 * n_rows))
    edges = phi_axis.edges
    widths = np.diff(edges)
    # Hex strings, not RGBA tuples: mplhep treats a sequence valued kwarg as
    # one entry per histogram and would index off the end of a single one
    colors = [matplotlib.colors.to_hex(c) for c in
              plt.cm.viridis(np.linspace(0., 0.85, n))]
    highest = 0.
    for i in range(n_rows * n_cols):
        ax = axes[i // n_cols][i % n_cols]
        if i >= n:
            ax.axis("off")
            continue
        label, counts, var = grouped[i]
        area = np.sum(counts * widths)
        if area == 0:
            continue
        mplhep.histplot(counts / area, edges, yerr=np.sqrt(var) / area,
                        histtype="step", edges=False, linewidth=1.5,
                        color=colors[i], ax=ax, label=label)
        highest = max(highest, (counts / area).max())
        ax.set_xlim(edges[0], edges[-1])
        ax.legend(title="$p_{T}(W)$", fontsize=13, title_fontsize=13)
        ax.tick_params(labelsize=13)
        if i // n_cols == n_rows - 1:
            ax.set_xlabel(xlabel, fontsize=15)
        if i % n_cols == 0:
            ax.set_ylabel("Normalised events", fontsize=15)
    axes[0][0].set_ylim(0, highest * 1.35)
    cms_label(axes[0][0], config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def summarise(centers, counts):
    """Mean and median of a distribution given as a histogram."""
    total = counts.sum()
    cdf = np.cumsum(counts) / total
    return {"mean": np.average(centers, weights=counts),
            "median": np.interp(0.5, cdf, centers),
            "overflow_free": total}


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


def plot(entries, axis, outfile, exts, config, log):
    """Overlay the two distributions as step lines, plus their sum, with the
    means in the legend."""
    fig, ax = plt.subplots()
    edges = axis.edges
    widths = np.diff(edges)
    combine = [e for e in entries if e[0] in ("$W^{+}$", "$W^{-}$")]
    if len(combine) > 1:
        entries = entries + [("$W^{+}$ and $W^{-}$", "black",
                              sum(e[2] for e in combine),
                              sum(e[3] for e in combine))]
    for label, color, counts, variances in entries:
        area = np.sum(counts * widths)
        if area == 0:
            continue
        mean = np.average(axis.centers, weights=counts)
        mplhep.histplot(counts / area, edges, yerr=np.sqrt(variances) / area,
                        histtype="step", edges=False, linewidth=1.5,
                        color=color, ax=ax,
                        label=f"{label}, mean {mean:.1f} GeV")
    ax.set_xlabel("$p_{T}(W)$ [GeV]")
    ax.set_ylabel("Normalised events")
    ax.set_xlim(edges[0], edges[-1])
    if log:
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.legend(fontsize=17)
    cms_label(ax, config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def main():
    parser = ArgumentParser(description="Plot the W transverse momentum")
    parser.add_argument("plot_config", help="Plotting config, used for the "
                        "CMS label only")
    parser.add_argument("histfile", help="hists.json written by runproc")
    parser.add_argument("-o", "--outdir", help="Output directory. Defaults to "
                        "the directory containing histfile")
    parser.add_argument("--ext", choices=["pdf", "svg", "png"], default=None,
                        action="append", help="Output file format")
    parser.add_argument("--log", action="store_true",
                        help="Logarithmic y axis, useful for the high pt tail")
    parser.add_argument("--split-datasets", action="store_true",
                        help="Make a separate plot for every dataset instead "
                        "of summing them, in a subdirectory per dataset")
    parser.add_argument("--variant", choices=sorted(PULL_VARIANTS),
                        default=None, action="append",
                        help="Which pull angle to bin in W pt. Can be given "
                        "several times. Defaults to all that are present.")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]

    config = Config(args.plot_config)
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    wanted = {v for v, _, _, _ in CURVES}
    cuts = sorted({k[0] for k in hists.keys() if k[1] in wanted})
    if len(cuts) == 0:
        raise SystemExit(f"None of {sorted(wanted)} found in {args.histfile}")

    # None means sum over all datasets
    datasets = [None]
    if args.split_datasets:
        # Union over the histograms used here: H_pt is filled for ttH only,
        # so one histogram alone would not list every dataset.
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
            print(f"Making a plot for each of: {datasets}")

    for cutname, dataset in ((c, d) for c in cuts for d in datasets):
        entries, axis = [], None
        for variable, label, color, _ in CURVES:
            got = load(hists, cutname, variable, dataset)
            if got is None:
                print(f"  MISSING {variable}, skipping")
                continue
            axis, counts, variances = got
            entries.append((label, color, counts, variances))

        if len(entries) == 0:
            continue

        print(f"[{cutname}]" if dataset is None
              else f"[{cutname}, {dataset}]")
        for label, _, counts, _ in entries:
            s = summarise(axis.centers, counts)
            print(f"  {label:<8} mean {s['mean']:7.1f} GeV   "
                  f"median {s['median']:7.1f} GeV")

        cut_outdir = os.path.join(outdir, cutname)
        if dataset is not None:
            cut_outdir = os.path.join(cut_outdir, dataset)
        os.makedirs(cut_outdir, exist_ok=True)
        plot(entries, axis, os.path.join(cut_outdir, f"w_pt_{cutname}"),
             args.ext, config, args.log)

        # Pull angle split into W pt bins, the two W charges added together
        for variant in (args.variant or sorted(PULL_VARIANTS)):
            prefix, xlabel = PULL_VARIANTS[variant]
            names = [f"{prefix}pull_angle_W{c}_vs_wpt"
                     for c in ("plus", "minus")]
            got = [load_2d(hists, cutname, v, dataset) for v in names]
            got = [g for g in got if g is not None]
            if len(got) == 0:
                continue
            wpt_axis, phi_axis = got[0][0], got[0][1]
            values = sum(g[2] for g in got)
            variances = sum(g[3] for g in got)
            grouped = group_wpt(wpt_axis, values, variances)
            plot_pull_vs_wpt(
                grouped, phi_axis,
                os.path.join(cut_outdir,
                             f"{prefix}pull_angle_vs_wpt_{cutname}"),
                args.ext, config, xlabel)
            print(f"  {xlabel}:")
            for label, counts, _ in grouped:
                if counts.sum() == 0:
                    continue
                central = (counts[np.abs(phi_axis.centers) < 1.].sum()
                           / counts.sum())
                print(f"    pt {label:<18}: "
                      f"{counts.sum()/values.sum():6.1%} of events, "
                      f"|pull angle| < 1 for {central:.1%}")

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
