#!/usr/bin/env python3
"""
Everything about the corridor between two jets.

The corridor comes in two shapes, both of half width R:
  ctr   the rectangle spanning the two jet centres, minus the parts inside
        either jet cone
  edge  the rectangle running from one cone edge to the other

Plots per cut, and per dataset with --split-datasets:
  corridor_ptdens_<cut>     one panel per radius, colour connected W pairs
                            split by W pt, with the unconnected pairs shown
                            over all pt for comparison
  corridor_shapes_<cut>     one panel per radius, the two corridor shapes
                            overlaid, to see how much the definition matters
  capsule_mass_W_<cut>      rows are corridor radius, columns are W pt bin,
                            m(two jets) against m(two jets + corridor) for the
                            W pairs, with a line at the W mass
  capsule_mass_H_<cut>      the same for the Higgs pair, binned in H pt and
                            with a line at the Higgs mass

Use --level reco for the versions built from reco jets and PF candidates.

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
H_MASS = 125.0
WINDOW = 10.0          # half width of the mass window used in the printout
PTDENS_YMAX = 0.25     # fixed so the panels can be compared by eye

RADII = [0.2, 0.3, 0.4]
DEFINITIONS = ["ctr", "edge"]
LEVEL_PREFIX = {"gen": "", "reco": "reco_"}

# W pt groups, or H pt groups for the Higgs pair. Boundaries must fall on
# edges of the fine binning used when filling.
WPT_GROUPS = [(0, 60), (60, 120), (120, 200), (200, None)]

# group suffix -> (legend name, reference mass, pt symbol)
PARENTS = {"connectedW": ("W", W_MASS, "$p_{T}(W)$"),
           "connectedH": ("H", H_MASS, "$p_{T}(H)$")}


def radius_tag(radius):
    return "r%02d" % round(radius * 10)


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
    """Return (pt axis, other axis, values, variances), or None."""
    key = (cutname, variable)
    if key not in hists.keys():
        return None
    h = integrate_categories(hists.load(key), dataset)
    if h is None:
        return None
    names = list(h.axes.name)
    if "wpt" not in names or len(names) != 2:
        return None
    other = [n for n in names if n != "wpt"][0]
    h = h.project("wpt", other)
    return h.axes["wpt"], h.axes[other], h.values(), h.variances()


def group_pt(pt_axis, values, variances):
    """Merge the fine pt bins into the groups given by WPT_GROUPS."""
    edges = pt_axis.edges
    grouped = []
    for lo, hi in WPT_GROUPS:
        top = edges[-1] if hi is None else hi
        for value in (lo, top):
            if not np.any(np.isclose(edges, value)):
                print(f"    WARNING: {value:g} GeV is not a bin edge")
        keep = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= top + 1e-9)
        if not keep.any():
            continue
        label = (f"$>$ {lo:.0f} GeV" if hi is None
                 else f"{lo:.0f} to {hi:.0f} GeV")
        grouped.append((label, values[keep].sum(axis=0),
                        variances[keep].sum(axis=0)))
    return grouped


def bin_colors(n):
    """Hex strings, not RGBA tuples: mplhep reads a sequence valued kwarg as
    one entry per histogram and would index off the end of a single one."""
    return [matplotlib.colors.to_hex(c)
            for c in plt.cm.viridis(np.linspace(0., 0.85, n))]


def roc_auc(signal, background):
    """Separation as one number. 0.5 means none, 1 means perfect."""
    sig = signal / signal.sum()
    bkg = background / background.sum()
    kept_sig = np.concatenate([[0.], np.cumsum(sig[::-1])])
    kept_bkg = np.concatenate([[0.], np.cumsum(bkg[::-1])])
    return np.trapezoid(kept_sig, kept_bkg) if hasattr(np, "trapezoid") \
        else np.trapz(kept_sig, kept_bkg)


def summarise(centers, counts, reference):
    """Peak, median, interquartile range and fraction near a reference mass.

    Median and IQR rather than mean and RMS, because these distributions have
    a long high mass tail that would dominate an RMS.
    """
    cdf = np.cumsum(counts) / counts.sum()

    def quantile(q):
        return np.interp(q, cdf, centers)

    return {"median": quantile(0.5),
            "iqr": quantile(0.75) - quantile(0.25),
            "in_window": counts[np.abs(centers - reference) < WINDOW].sum()
            / counts.sum()}


def cms_label(ax, config, fontsize=16):
    if "year" not in config:
        return
    label_kwargs = {"fontsize": fontsize}
    cmslabel = config["cmslabel"] if "cmslabel" in config else None
    if cmslabel is not None and cmslabel.strip().lower() != "simulation":
        label_kwargs["label"] = cmslabel
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


def panel_title(ax, text):
    """Above the axes, high enough to clear the CMS header."""
    ax.set_title(text, fontsize=17, pad=28)


def plot_ptdens_panels(per_radius, outfile, exts, config):
    """One panel per radius, connected pairs split by pt plus the unconnected
    ones over all pt. Shared y scale so the panels can be compared."""
    n = len(per_radius)
    fig, axes = plt.subplots(1, n, figsize=(6.5 * n, 6.5), squeeze=False,
                             sharey=True)
    for ax, (radius, axis, grouped, unconn) in zip(axes[0], per_radius):
        edges = axis.edges
        colors = bin_colors(len(grouped))
        for i, (label, counts, var) in enumerate(grouped):
            step(ax, counts, var, edges, linewidth=1.5, color=colors[i],
                 label=label)
        if unconn is not None:
            step(ax, unconn[0], unconn[1], edges, linewidth=2, linestyle=":",
                 color="black", label="Not connected, all $p_{T}$")
        panel_title(ax, f"Corridor R = {radius}")
        ax.set_xlabel("Corridor $\\Sigma p_{T}$ per unit area [GeV]",
                      fontsize=16)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0, PTDENS_YMAX)
        ax.tick_params(labelsize=14)
    axes[0][0].set_ylabel("Normalised pairs", fontsize=16)
    axes[0][0].legend(title="$p_{T}$ of the parent", fontsize=12,
                      title_fontsize=12)
    cms_label(axes[0][0], config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def plot_shape_panels(per_radius, outfile, exts, config):
    """One panel per radius, the two corridor shapes overlaid."""
    n = len(per_radius)
    fig, axes = plt.subplots(1, n, figsize=(6.5 * n, 6.5), squeeze=False,
                             sharey=True)
    styles = {"ctr": ("-", "tab:blue", "Centre to centre, cones removed"),
              "edge": ("--", "tab:orange", "Edge to edge")}
    for ax, (radius, axis, curves) in zip(axes[0], per_radius):
        edges = axis.edges
        for definition, counts, var in curves:
            ls, color, name = styles[definition]
            mean = np.average(axis.centers, weights=counts)
            step(ax, counts, var, edges, linewidth=1.5, linestyle=ls,
                 color=color, label=f"{name}, mean {mean:.2f}")
        panel_title(ax, f"Corridor R = {radius}")
        ax.set_xlabel("Corridor $\\Sigma p_{T}$ per unit area [GeV]",
                      fontsize=16)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0, PTDENS_YMAX)
        ax.tick_params(labelsize=14)
    axes[0][0].set_ylabel("Normalised pairs", fontsize=16)
    axes[0][0].legend(fontsize=12)
    cms_label(axes[0][0], config)
    plt.tight_layout()
    for ext in exts:
        fig.savefig(outfile + "." + ext)
    plt.close(fig)


def plot_capsule_grid(grid, axis, parent, reference, outfile, exts, config):
    """A panel per (corridor radius, pt bin), comparing the two masses."""
    radii = list(grid)
    n_rows, n_cols = len(radii), max(len(v) for v in grid.values())
    fig, axes = plt.subplots(n_rows, n_cols, squeeze=False, sharex=True,
                             figsize=(4.6 * n_cols, 4.4 * n_rows))
    edges = axis.edges
    colors = bin_colors(n_cols)
    for r, radius in enumerate(radii):
        for c in range(n_cols):
            ax = axes[r][c]
            if c >= len(grid[radius]):
                ax.axis("off")
                continue
            label, d_counts, d_var, c_counts, c_var = grid[radius][c]
            step(ax, d_counts, d_var, edges, linewidth=1.5, linestyle="--",
                 color=colors[c], label="Jets only")
            step(ax, c_counts, c_var, edges, linewidth=1.5, linestyle="-",
                 color=colors[c], label="Jets + corridor")
            ax.axvline(reference, color="black", linestyle=":", linewidth=1.2)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(bottom=0)
            ax.tick_params(labelsize=12)
            if r == 0:
                ax.set_title(label, fontsize=15)
            if c == 0:
                ax.set_ylabel(f"R = {radius}", fontsize=15)
            if r == n_rows - 1:
                ax.set_xlabel("Invariant mass [GeV]", fontsize=14)
    axes[0][0].legend(fontsize=11)
    fig.suptitle(f"Jet Mass and Capsule Mass for {parent} jets",
                 fontsize=17, y=0.999)
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
    parser.add_argument("--level", choices=["gen", "reco"], default=None,
                        action="append",
                        help="Which jets to use. Defaults to both.")
    parser.add_argument("--split-datasets", action="store_true",
                        help="Make a separate set of plots for every dataset "
                        "instead of summing them, in a subdirectory per "
                        "dataset")
    args = parser.parse_args()

    if args.ext is None:
        args.ext = ["pdf"]
    levels = args.level or ["gen", "reco"]

    config = Config(args.plot_config)
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.realpath(args.histfile))

    cuts = sorted({k[0] for k in hists.keys() if "corridor_" in k[1]})
    if len(cuts) == 0:
        raise SystemExit(f"No corridor histograms found in {args.histfile}")

    datasets = [None]
    if args.split_datasets:
        # Union over the corridor histograms: the connectedH ones are filled
        # for ttH only, so one histogram alone would not list every dataset.
        found = []
        for key in hists.keys():
            if "corridor_" not in key[1]:
                continue
            axes = hists.load(key).axes
            if "dataset" not in axes.name:
                continue
            found += [d for d in axes["dataset"] if d not in found]
        if found:
            datasets = found
            print(f"Making a set of plots for each of: {datasets}")

    for cutname, dataset in ((c, d) for c in cuts for d in datasets):
        head = cutname if dataset is None else f"{cutname}, {dataset}"
        cut_outdir = os.path.join(outdir, cutname)
        if dataset is not None:
            cut_outdir = os.path.join(cut_outdir, dataset)
        os.makedirs(cut_outdir, exist_ok=True)

        for level in levels:
            pre = LEVEL_PREFIX[level]
            print(f"[{head}] {level} level")
            made_any = False

            # --- pt density, and the two shapes side by side ---
            ptdens_panels, shape_panels = [], []
            for radius in RADII:
                rt = radius_tag(radius)
                unconn = load(hists, cutname,
                              f"{pre}corridor_ptdens_unconnected_ctr_{rt}",
                              dataset)
                conn = load(hists, cutname,
                            f"{pre}corridor_ptdens_connectedW_ctr_{rt}",
                            dataset)
                if conn is not None and unconn is not None:
                    print(f"  R = {radius}  mean density  connected "
                          f"{np.average(conn[0].centers, weights=conn[1]):6.2f}"
                          f"   unconnected "
                          f"{np.average(unconn[0].centers, weights=unconn[1]):6.2f}"
                          f"   AUC {roc_auc(conn[1], unconn[1]):.4f}")
                binned = load_2d(
                    hists, cutname,
                    f"{pre}corridor_ptdens_connectedW_ctr_{rt}_vs_pt", dataset)
                if binned is not None:
                    pt_axis, axis, values, variances = binned
                    ptdens_panels.append((
                        radius, axis, group_pt(pt_axis, values, variances),
                        (unconn[1], unconn[2]) if unconn is not None else None))
                shapes = []
                for definition in DEFINITIONS:
                    got = load(
                        hists, cutname,
                        f"{pre}corridor_ptdens_connectedW_{definition}_{rt}",
                        dataset)
                    if got is not None:
                        shapes.append((definition, got[1], got[2]))
                        axis_s = got[0]
                if len(shapes) == 2:
                    shape_panels.append((radius, axis_s, shapes))

            if ptdens_panels:
                plot_ptdens_panels(
                    ptdens_panels,
                    os.path.join(cut_outdir, f"{pre}corridor_ptdens_{cutname}"),
                    args.ext, config)
                made_any = True
            if shape_panels:
                plot_shape_panels(
                    shape_panels,
                    os.path.join(cut_outdir, f"{pre}corridor_shapes_{cutname}"),
                    args.ext, config)
                made_any = True

            # --- capsule mass, once for the W pairs and once for the H pair ---
            for group, (parent, reference, ptsym) in PARENTS.items():
                grid, mass_axis = {}, None
                for radius in RADII:
                    rt = radius_tag(radius)
                    dj = load_2d(hists, cutname,
                                 f"{pre}corridor_dijet_{group}_ctr_{rt}_vs_pt",
                                 dataset)
                    cp = load_2d(hists, cutname,
                                 f"{pre}corridor_capsule_{group}_ctr_{rt}_vs_pt",
                                 dataset)
                    if dj is None or cp is None:
                        continue
                    mass_axis = dj[1]
                    d_grouped = group_pt(dj[0], dj[2], dj[3])
                    c_grouped = group_pt(cp[0], cp[2], cp[3])
                    grid[radius] = [(d[0], d[1], d[2], c[1], c[2])
                                    for d, c in zip(d_grouped, c_grouped)]
                    for d, c in zip(d_grouped, c_grouped):
                        if d[1].sum() == 0 or c[1].sum() == 0:
                            continue
                        sd = summarise(mass_axis.centers, d[1], reference)
                        sc = summarise(mass_axis.centers, c[1], reference)
                        print(f"  {parent} R={radius} {d[0]:<16} "
                              f"median {sd['median']:6.1f} -> "
                              f"{sc['median']:6.1f}   IQR {sd['iqr']:5.1f} -> "
                              f"{sc['iqr']:5.1f} "
                              f"({'sharper' if sc['iqr'] < sd['iqr'] else 'broader'})")
                if grid and mass_axis is not None:
                    os.makedirs(cut_outdir, exist_ok=True)
                    plot_capsule_grid(
                        grid, mass_axis, parent, reference,
                        os.path.join(cut_outdir,
                                     f"{pre}capsule_mass_{parent}_{cutname}"),
                        args.ext, config)
                    made_any = True

            if not made_any:
                print(f"  nothing found for the {level} level")

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
