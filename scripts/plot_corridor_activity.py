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

# Corridor mass, connected against unconnected
CORRIDOR_CURVES = [
    ("corridor_mass_connected", "Colour connected", "tab:blue"),
    ("corridor_mass_unconnected", "Not colour connected", "tab:red"),
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


def integrate_categories(h):
    """Sum over every category axis, keeping only the nominal systematic."""
    if "sys" in h.axes.name:
        h = h[{"sys": "nominal"}]
    for ax in list(h.axes):
        if isinstance(ax, (hist.axis.StrCategory, hist.axis.IntCategory)):
            h = h.integrate(ax.name, list(ax))
    return h


def load(hists, cutname, variable):
    """Return (dense axis, counts, variances), or None if not present."""
    key = (cutname, variable)
    if key not in hists.keys():
        return None
    h = integrate_categories(hists.load(key))
    dense = [a for a in h.axes
             if not isinstance(a, (hist.axis.StrCategory, hist.axis.IntCategory))]
    if len(dense) != 1:
        return None
    return dense[0], h.values(), h.variances()


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


def overlay(entries, axis, xlabel, outfile, exts, config, wline=False):
    """Overlay normalised distributions as step lines."""
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

    for cutname in cuts:
        print(f"[{cutname}]")
        cut_outdir = os.path.join(outdir, cutname)
        os.makedirs(cut_outdir, exist_ok=True)

        # --- how many pairs lost a third jet in their corridor ---
        for variable, name in VETO_VARIABLES:
            got = load(hists, cutname, variable)
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
            got = load(hists, cutname, variable)
            if got is None:
                continue
            axis, counts, variances = got
            entries.append((label, color, counts, variances))
        if len(entries) == 2:
            overlay(entries, axis, "$m$(corridor particles) [GeV]",
                    os.path.join(cut_outdir, f"corridor_mass_{cutname}"),
                    args.ext, config)
            auc = roc_auc(entries[0][2], entries[1][2])
            for label, _, counts, _ in entries:
                print(f"  mean corridor mass, {label:<22}: "
                      f"{np.average(axis.centers, weights=counts):7.2f} GeV")
            print(f"  separation, area under ROC     : {auc:.4f}"
                  f"   (0.5 = none)")

        # --- two jets against two jets plus corridor ---
        entries, axis = [], None
        for variable, label, color in CAPSULE_CURVES:
            got = load(hists, cutname, variable)
            if got is None:
                continue
            axis, counts, variances = got
            entries.append((label, color, counts, variances))
        if len(entries) == 2:
            overlay(entries, axis, "Invariant mass [GeV]",
                    os.path.join(cut_outdir, f"capsule_mass_{cutname}"),
                    args.ext, config, wline=True)
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

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
