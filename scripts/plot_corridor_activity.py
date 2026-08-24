#!/usr/bin/env python3
"""
Everything about the corridor between two jets, at several corridor radii.

The corridor is the rectangle in the (rapidity, phi) plane spanning the two
jet centres, of half width R, with the parts inside either jet cone removed.

Plots per cut, and per dataset with --split-datasets:
  corridor_ptdens_<cut>          one panel per corridor radius, colour
                                 connected pairs split by W pt, with the
                                 unconnected pairs over all pt for comparison
  capsule_mass_<cut>             a grid of panels, one row per corridor radius
                                 and one column per W pt bin, comparing
                                 m(two jets) with m(two jets + corridor)

Also prints, for every radius: the fraction of pairs vetoed for containing a
third jet, the separation between connected and unconnected as an area under
ROC, and the peak, width and W window fraction of the two masses.

Usage
-----
    python scripts/plot_corridor_activity.py configs/plot_config.json \\
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

W_MASS = 80.4
WINDOW = 10.0   # half width of the window around the W mass, in GeV

# Corridor radii to look for, matching corridor_radii in the processor
RADII = [0.2, 0.3, 0.4]

# W pt groups, one column of the capsule grid and one curve of the pt density
# plot each. Boundaries must fall on edges of the fine binning used to fill.
WPT_GROUPS = [(0, 60), (60, 120), (120, 200), (200, None)]


def radius_tag(radius):
    return "r%02d" % round(radius * 10)


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
             if not isinstance(a, (hist.axis.StrCategory,
                                   hist.axis.IntCategory))]
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


def bin_colors(n):
    """One colour per W pt bin. Hex strings, not RGBA tuples, because mplhep
    reads a sequence valued kwarg as one entry per histogram."""
    return [matplotlib.colors.to_hex(c)
            for c in plt.cm.viridis(np.linspace(0., 0.85, n))]


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


def cms_label(ax, config, fontsize=None):
    if "year" not in config:
        return
    label_kwargs = {}
    cmslabel = config["cmslabel"] if "cmslabel" in config else None
    if cmslabel is not None and cmslabel.strip().lower() != "simulation":
        label_kwargs["label"] = cmslabel
    if fontsize is not None:
        label_kwargs["fontsize"] = fontsize
    mplhep.cms.label(
        ax=ax, data=False, year=config["year"],
        lumi=f"{config['luminosity']:.1f}" if "luminosity" in config else None,
        com=config["com_energy"] if "com_energy" in config else 13,
        **label_kwargs)


def step(ax, counts, variances, edges, **kwargs):
    """Draw one normalised step line, doing nothing if it is empty."""
    area = np.sum(counts * np.diff(edges))
    if area == 0:
        return 0.
    mplhep.histplot(counts / area, edges, yerr=np.sqrt(variances) / area,
                    histtype="step", edges=False, ax=ax, **kwargs)
    return (counts / area).max()


def plot_ptdens_panels(per_radius, outfile, exts, config):
    """One panel per corridor radius. In each, the connected pairs split by
    W pt plus the unconnected pairs over all pt.

    `per_radius` is a list of (radius, axis, grouped, unconnected or None).
    """
    n = len(per_radius)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7), squeeze=False)
    for ax, (radius, axis, grouped, unconn) in zip(axes[0], per_radius):
        edges = axis.edges
        colors = bin_colors(len(grouped))
        highest = 0.
        for i, (label, counts, var) in enumerate(grouped):
            highest = max(highest, step(ax, counts, var, edges, linewidth=1.5,
                                        color=colors[i], label=label))
        if unconn is not None:
            highest = max(highest, step(
                ax, unconn[0], unconn[1], edges, linewidth=2, linestyle=":",
                color="black", label="Not connected, all $p_{T}$"))
        ax.set_title(f"Corridor R = {radius}", fontsize=17)
        ax.set_xlabel("Corridor $\\Sigma p_{T}$ per unit area [GeV]",
                      fontsize=17)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0, highest * 1.35)
        ax.tick_params(labelsize=14)
    axes[0][0].set_ylabel("Normalised pairs", fontsize=17)
    axes[0][0].legend(title="$p_{T}(W)$", fontsize=12, title_fontsize=12)
    cms_label(axes[0][0], config, fontsize=16)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def plot_capsule_grid(grid, axis, outfile, exts, config):
    """A panel per (corridor radius, W pt bin), comparing the two masses.

    `grid` maps radius to a list over W pt bins of
    (label, dijet counts, dijet variances, capsule counts, capsule variances).
    """
    radii = list(grid)
    n_rows, n_cols = len(radii), max(len(v) for v in grid.values())
    fig, axes = plt.subplots(n_rows, n_cols, squeeze=False,
                             figsize=(4.6 * n_cols, 4.4 * n_rows),
                             sharex=True)
    edges = axis.edges
    colors = bin_colors(n_cols)
    for r, radius in enumerate(radii):
        for c in range(n_cols):
            ax = axes[r][c]
            if c >= len(grid[radius]):
                ax.axis("off")
                continue
            label, d_counts, d_var, c_counts, c_var = grid[radius][c]
            highest = max(
                step(ax, d_counts, d_var, edges, linewidth=1.5,
                     linestyle="--", color=colors[c], label="Jets only"),
                step(ax, c_counts, c_var, edges, linewidth=1.5,
                     linestyle="-", color=colors[c], label="Jets + corridor"))
            ax.axvline(W_MASS, color="black", linestyle=":", linewidth=1.2)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0, highest * 1.45)
            ax.tick_params(labelsize=12)
            if r == 0:
                ax.set_title(label, fontsize=15)
            if c == 0:
                ax.set_ylabel(f"R = {radius}", fontsize=15)
            if r == n_rows - 1:
                ax.set_xlabel("Invariant mass [GeV]", fontsize=14)
    axes[0][0].legend(fontsize=11)
    fig.suptitle("Columns: $p_{T}(W)$ bin.   Rows: corridor radius",
                 fontsize=15, y=0.998)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def main():
    parser = ArgumentParser(description="Corridor plots at several radii")
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

    cuts = sorted({k[0] for k in hists.keys() if k[1].startswith("corridor_")})
    if len(cuts) == 0:
        raise SystemExit(f"No corridor histograms found in {args.histfile}")

    datasets = [None]
    if args.split_datasets:
        first = next(k for k in hists.keys() if k[1].startswith("corridor_"))
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

        ptdens_panels, capsule_grid, mass_axis = [], {}, None
        for radius in RADII:
            tag = radius_tag(radius)
            print(f"  corridor R = {radius}")

            for group in ("connected", "unconnected"):
                got = load(hists, cutname, f"corridor_vetoed_{group}_{tag}",
                           dataset)
                if got is None:
                    continue
                _, counts, variances = got
                total = counts.sum()
                if total > 0:
                    n_eff = total**2 / variances.sum()
                    print(f"    vetoed, {group:<12}: {counts[1]/total:6.1%} "
                          f"of {n_eff:.0f} pairs with a corridor")

            conn = load(hists, cutname, f"corridor_ptdens_connected_{tag}",
                        dataset)
            unconn = load(hists, cutname,
                          f"corridor_ptdens_unconnected_{tag}", dataset)
            if conn is not None and unconn is not None:
                print(f"    mean pt density, connected  : "
                      f"{np.average(conn[0].centers, weights=conn[1]):7.2f} GeV")
                print(f"    mean pt density, unconnected: "
                      f"{np.average(unconn[0].centers, weights=unconn[1]):7.2f} GeV")
                print(f"    separation, area under ROC  : "
                      f"{roc_auc(conn[1], unconn[1]):.4f}   (0.5 = none)")

            binned = load_2d(hists, cutname,
                             f"corridor_ptdens_connected_{tag}_vs_wpt", dataset)
            if binned is not None:
                wpt_axis, axis, values, variances = binned
                ptdens_panels.append((
                    radius, axis, group_wpt(wpt_axis, values, variances),
                    (unconn[1], unconn[2]) if unconn is not None else None))
            else:
                print(f"    MISSING corridor_ptdens_connected_{tag}_vs_wpt")

            dj = load_2d(hists, cutname,
                         f"corridor_dijet_connected_{tag}_vs_wpt", dataset)
            cp = load_2d(hists, cutname,
                         f"corridor_capsule_connected_{tag}_vs_wpt", dataset)
            if dj is None or cp is None:
                print(f"    MISSING the pt binned masses for {tag}")
                continue
            mass_axis = dj[1]
            d_grouped = group_wpt(dj[0], dj[2], dj[3])
            c_grouped = group_wpt(cp[0], cp[2], cp[3])
            capsule_grid[radius] = [
                (d[0], d[1], d[2], c[1], c[2])
                for d, c in zip(d_grouped, c_grouped)]
            for d, c in zip(d_grouped, c_grouped):
                if d[1].sum() == 0 or c[1].sum() == 0:
                    continue
                sd = summarise(mass_axis.centers, d[1])
                sc = summarise(mass_axis.centers, c[1])
                print(f"    {d[0]:<16} median {sd['median']:6.1f} -> "
                      f"{sc['median']:6.1f} GeV, IQR {sd['iqr']:5.1f} -> "
                      f"{sc['iqr']:5.1f} "
                      f"({'sharper' if sc['iqr'] < sd['iqr'] else 'broader'})")

        if len(ptdens_panels) > 0:
            plot_ptdens_panels(
                ptdens_panels,
                os.path.join(cut_outdir, f"corridor_ptdens_{cutname}"),
                args.ext, config)
        if len(capsule_grid) > 0 and mass_axis is not None:
            plot_capsule_grid(
                capsule_grid, mass_axis,
                os.path.join(cut_outdir, f"capsule_mass_{cutname}"),
                args.ext, config)

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
