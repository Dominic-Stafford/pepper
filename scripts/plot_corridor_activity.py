#!/usr/bin/env python3
"""
Everything about the corridor between two jets.

Makes two plots per cut:
  corridor_mass_<cut>    invariant mass of the particles in the corridor, for
                         colour connected and for unconnected jet pairs
  capsule_mass_<cut>     m(two jets) against m(two jets + corridor) for
                         connected pairs, to see whether sweeping in the soft
                         particles sharpens the W peak

Also prints, for each cut:
  - how large a fraction of pairs was vetoed for having a third jet in the
    corridor
  - the separation between connected and unconnected, as an area under ROC
  - the peak, width and W window fraction of the two mass distributions

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

W_MASS = 80.4
WINDOW = 10.0   # half width of the window around the W mass, in GeV

# Corridor pt density, connected against unconnected
CORRIDOR_CURVES = [
    ("corridor_ptdens_connected", "Colour connected", "tab:blue"),
    ("corridor_ptdens_unconnected", "Not colour connected", "tab:red"),
]
# The two jets alone against the two jets plus the corridor
CAPSULE_CURVES = [
    ("corridor_dijet_connected", "Two jets", "tab:blue"),
    ("corridor_capsule_connected", "Two jets + corridor", "tab:orange"),
]
VETO_VARIABLES = [
    ("corridor_vetoed_connected", "connected"),
    ("corridor_vetoed_unconnected", "not connected"),
]

# W pt groups for the binned plots. Only the colour connected pairs are split
# this way, since an unconnected pair has no W to bin it by. Boundaries must
# fall on edges of the fine binning used when filling.
WPT_GROUPS = [(0, 60), (60, 120), (120, 200), (200, None)]
# Each binned plot: (title, list of (variable, line style, legend suffix),
#                    optional inclusive extra curve, x axis label)
BINNED_PLOTS = [
    ("corridor_ptdens_vs_wpt",
     [("corridor_ptdens_connected_vs_wpt", "-", "")],
     ("corridor_ptdens_unconnected", "Not connected, all $p_{T}$"),
     "Corridor $\\Sigma p_{T}$ per unit area [GeV]"),
    ("capsule_mass_vs_wpt",
     [("corridor_capsule_connected_vs_wpt", "-", " (jets + corridor)"),
      ("corridor_dijet_connected_vs_wpt", "--", " (jets only)")],
     None,
     "Invariant mass [GeV]"),
]


def integrate_categories(h, dataset=None):
    """Sum over every category axis, keeping only the nominal systematic.
    If `dataset` is given, keep only that dataset instead of summing them."""
    if dataset is not None and "dataset" in h.axes.name:
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
    dense = [a for a in h.axes
             if not isinstance(a, (hist.axis.StrCategory, hist.axis.IntCategory))]
    if len(dense) != 1:
        return None
    return dense[0], h.values(), h.variances()


def load_2d(hists, cutname, variable, dataset=None):
    """Return (wpt axis, other axis, values, variances), or None."""
    key = (cutname, variable)
    if key not in hists.keys():
        return None
    h = integrate_categories(hists.load(key), dataset)
    names = list(h.axes.name)
    if "wpt" not in names or len(names) != 2:
        return None
    other = [n for n in names if n != "wpt"][0]
    h = h.project("wpt", other)
    return h.axes["wpt"], h.axes[other], h.values(), h.variances()


def group_wpt(wpt_axis, values, variances):
    """Merge the fine W pt bins into the groups given by WPT_GROUPS."""
    edges = wpt_axis.edges
    grouped = []
    for lo, hi in WPT_GROUPS:
        top = edges[-1] if hi is None else hi
        for value in (lo, top):
            if not np.any(np.isclose(edges, value)):
                print(f"  WARNING: {value:g} GeV is not a bin edge")
        keep = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= top + 1e-9)
        if not keep.any():
            continue
        label = (f"$>$ {lo:.0f} GeV" if hi is None
                 else f"{lo:.0f} to {hi:.0f} GeV")
        grouped.append((label, values[keep].sum(axis=0),
                        variances[keep].sum(axis=0)))
    return grouped


def plot_binned(curves, extra, axis, xlabel, outfile, exts, config, wline):
    """Colour for the W pt bin, line style for the kind of curve.

    `curves` is a list of (style, suffix, grouped) and `extra` an optional
    (label, counts, variances) drawn in black over all pt.
    """
    fig, ax = plt.subplots()
    edges = axis.edges
    widths = np.diff(edges)
    n_bins = max(len(g) for _, _, g in curves)
    # Hex strings, not RGBA: mplhep reads a sequence kwarg as one per histogram
    colors = [matplotlib.colors.to_hex(c) for c in
              plt.cm.viridis(np.linspace(0., 0.85, n_bins))]
    for style, suffix, grouped in curves:
        for i, (label, counts, var) in enumerate(grouped):
            area = np.sum(counts * widths)
            if area == 0:
                continue
            mplhep.histplot(counts / area, edges, yerr=np.sqrt(var) / area,
                            histtype="step", edges=False, linewidth=1.5,
                            linestyle=style, color=colors[i], ax=ax,
                            label=f"{label}{suffix}")
    if extra is not None:
        label, counts, var = extra
        area = np.sum(counts * widths)
        if area > 0:
            mplhep.histplot(counts / area, edges, yerr=np.sqrt(var) / area,
                            histtype="step", edges=False, linewidth=2,
                            linestyle=":", color="black", ax=ax, label=label)
    if wline:
        ax.axvline(W_MASS, color="black", linestyle="--", linewidth=1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalised pairs")
    ax.set_xlim(edges[0], edges[-1])
    # Enough headroom for a legend that can run to eight entries
    highest = max((c / np.sum(c * widths)).max()
                  for _, _, g in curves for _, c, _ in g
                  if np.sum(c * widths) > 0)
    ax.set_ylim(0, highest * (1.75 if len(curves) > 1 else 1.45))
    ax.legend(title="$p_{T}(W)$", fontsize=12, title_fontsize=12,
              ncol=2 if len(curves) > 1 else 1, loc="upper left")
    cms_label(ax, config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def roc_auc(signal, background):
    """How well a cut on this observable separates the two, as one number.
    0.5 means no separation, 1 means perfect."""
    sig = signal / signal.sum()
    bkg = background / background.sum()
    kept_sig = np.concatenate([[0.], np.cumsum(sig[::-1])])
    kept_bkg = np.concatenate([[0.], np.cumsum(bkg[::-1])])
    return np.trapezoid(kept_sig, kept_bkg) if hasattr(np, "trapezoid") \
        else np.trapz(kept_sig, kept_bkg)


def summarise(centers, counts):
    """Peak, median, interquartile range and fraction near the W mass.

    Median and IQR rather than mean and RMS, because these distributions have
    a long high mass tail that would dominate an RMS.
    """
    cdf = np.cumsum(counts) / counts.sum()

    def quantile(q):
        return np.interp(q, cdf, centers)

    return {"peak": centers[np.argmax(counts)], "median": quantile(0.5),
            "iqr": quantile(0.75) - quantile(0.25),
            "in_window": counts[np.abs(centers - W_MASS) < WINDOW].sum()
            / counts.sum()}


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


def overlay(entries, axis, xlabel, outfile, exts, config, wline=False,
            show_mean=False):
    """Overlay normalised distributions as step lines. With `show_mean` the
    mean of each distribution is added to its legend entry."""
    fig, ax = plt.subplots()
    edges = axis.edges
    widths = np.diff(edges)
    for label, color, counts, variances in entries:
        area = np.sum(counts * widths)
        if area == 0:
            continue
        if show_mean:
            mean = np.average(axis.centers, weights=counts)
            label = f"{label}, mean {mean:.1f} GeV"
        mplhep.histplot(counts / area, edges, yerr=np.sqrt(variances) / area,
                        histtype="step", edges=False, linewidth=1.5,
                        color=color, label=label, ax=ax)
    if wline:
        ax.axvline(W_MASS, color="black", linestyle="--", linewidth=1.5)
        ax.text(W_MASS + 2, 0.97, f"$m_W$ = {W_MASS} GeV", fontsize=15,
                va="top", ha="left", transform=ax.get_xaxis_transform())
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalised pairs")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=17)
    cms_label(ax, config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def main():
    parser = ArgumentParser(description="Corridor mass and capsule mass plots")
    parser.add_argument("plot_config", help="Plotting config, used for the "
                        "CMS label only")
    parser.add_argument("histfile", help="hists.json written by runproc")
    parser.add_argument("-o", "--outdir", help="Output directory. Defaults to "
                        "the directory containing histfile")
    parser.add_argument("--ext", choices=["pdf", "svg", "png"], default=None,
                        action="append", help="Output file format")
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

    wanted = ({v for v, _, _ in CORRIDOR_CURVES}
              | {v for v, _, _ in CAPSULE_CURVES}
              | {v for v, _ in VETO_VARIABLES})
    cuts = sorted({k[0] for k in hists.keys() if k[1] in wanted})
    if len(cuts) == 0:
        raise SystemExit(f"No corridor histograms found in {args.histfile}")

    # None means sum over all datasets
    datasets = [None]
    if args.split_datasets:
        first = next(k for k in hists.keys() if k[1] in wanted)
        axes = hists.load(first).axes
        if "dataset" in axes.name:
            datasets = list(axes["dataset"])
            print(f"Making a set of plots for each of: {datasets}")

    for cutname, dataset in ((c, d) for c in cuts for d in datasets):
        print(f"[{cutname}]" if dataset is None else f"[{cutname}, {dataset}]")
        cut_outdir = os.path.join(outdir, cutname)
        if dataset is not None:
            cut_outdir = os.path.join(cut_outdir, dataset)
        os.makedirs(cut_outdir, exist_ok=True)

        # --- how many pairs lost a third jet in their corridor ---
        for variable, name in VETO_VARIABLES:
            got = load(hists, cutname, variable, dataset)
            if got is None:
                continue
            _, counts, variances = got
            total = counts.sum()
            if total > 0:
                # Entries are weighted, so use the effective number of pairs
                n_eff = total**2 / variances.sum()
                print(f"  vetoed, {name:<20}: {counts[1]/total:6.1%} of "
                      f"{n_eff:.0f} pairs that had a corridor")

        # --- corridor mass, connected against unconnected ---
        entries, axis = [], None
        for variable, label, color in CORRIDOR_CURVES:
            got = load(hists, cutname, variable, dataset)
            if got is None:
                print(f"  MISSING {variable}, skipping the corridor mass plot")
                continue
            axis, counts, variances = got
            entries.append((label, color, counts, variances))
        if len(entries) == 2:
            overlay(entries, axis,
                    "Corridor $\\Sigma p_{T}$ per unit area [GeV]",
                    os.path.join(cut_outdir, f"corridor_ptdens_{cutname}"),
                    args.ext, config)
            auc = roc_auc(entries[0][2], entries[1][2])
            for label, _, counts, _ in entries:
                print(f"  mean corridor pt density, {label:<22}: "
                      f"{np.average(axis.centers, weights=counts):7.2f} GeV")
            print(f"  separation, area under ROC     : {auc:.4f}"
                  f"   (0.5 = none)")

        # --- two jets against two jets plus corridor ---
        entries, axis = [], None
        for variable, label, color in CAPSULE_CURVES:
            got = load(hists, cutname, variable, dataset)
            if got is None:
                print(f"  MISSING {variable}, skipping the capsule mass plot")
                continue
            axis, counts, variances = got
            entries.append((label, color, counts, variances))
        if len(entries) == 2:
            overlay(entries, axis, "Invariant mass [GeV]",
                    os.path.join(cut_outdir, f"capsule_mass_{cutname}"),
                    args.ext, config, wline=True, show_mean=True)
            stats = {}
            for label, _, counts, _ in entries:
                s = summarise(axis.centers, counts)
                stats[label] = s
                print(f"  {label:<22} peak {s['peak']:6.1f}   "
                      f"median {s['median']:6.1f}   IQR {s['iqr']:6.1f}   "
                      f"within {WINDOW:.0f} GeV of W: {s['in_window']:.1%}")
            first, second = entries[0][0], entries[1][0]
            d_median = stats[second]["median"] - stats[first]["median"]
            d_iqr = stats[second]["iqr"] - stats[first]["iqr"]
            print(f"  -> median moves {d_median:+.1f} GeV, IQR changes "
                  f"{d_iqr:+.1f} GeV "
                  f"({'sharper' if d_iqr < 0 else 'broader'})")

        # --- the same quantities split into W pt bins ---
        for name, specs, extra_spec, xlabel in BINNED_PLOTS:
            curves, axis = [], None
            for variable, style, suffix in specs:
                got = load_2d(hists, cutname, variable, dataset)
                if got is None:
                    print(f"  MISSING {variable}, skipping {name}")
                    continue
                wpt_axis, axis, values, variances = got
                curves.append((style, suffix,
                               group_wpt(wpt_axis, values, variances)))
            if len(curves) == 0:
                continue
            extra = None
            if extra_spec is not None:
                got = load(hists, cutname, extra_spec[0], dataset)
                if got is not None:
                    extra = (extra_spec[1], got[1], got[2])
            plot_binned(curves, extra, axis, xlabel,
                        os.path.join(cut_outdir, f"{name}_{cutname}"),
                        args.ext, config, wline="mass" in name)

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
