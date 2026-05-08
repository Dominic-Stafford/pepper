"""
Plotting script for pepper cutflow output.
It uses as an input the cutflows.json file produced by pepper.

While it is quite flexible, it is not meant as a one-fits-all
solution, but instead should be extended and adapted to the needs
of your analysis.

Structure of the code:
- A helper functions for cutflow detection and extraction is defined at the top.
- Plotting functions are defined in the second half, including:
    - `plot_stacked_data_mc_ratio`: A comparison plot between Data and stacked MC with a ratio panel.
        You need to manually specify the keys that correspond to the Data samples.
    - `plot_individual_efficiencies`: Individual efficiency plots for each sample.
- File I/O and categorization is done in the `main` function, which reads the JSON, extracts cutflows,
    and calls the plotting functions.
"""

import json
import logging
from argparse import ArgumentParser
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

# --------------------------------------------------


def extract_cutflow(sample_dict, level="all"):

    if level not in sample_dict:
        raise ValueError(f"Level '{level}' not found in sample. Available levels: {list(sample_dict.keys())}")

    level_dict = sample_dict[level]
    cutflow = OrderedDict()

    for cut, hist_dict in level_dict.items():
        try:
            if isinstance(hist_dict, dict):
                total = 0.0
                stack = [hist_dict]

                while stack:
                    item = stack.pop()
                    if isinstance(item, dict):
                        stack.extend(item.values())
                    elif isinstance(item, (int, float)):
                        total += item

                cutflow[cut] = total
            else:
                cutflow[cut] = float(hist_dict)

        except Exception:
            logger.exception(f"Failed to extract cut '{cut}'")
            continue
    return cutflow if cutflow else None

# --------------------------------------------------


def plot_stacked_data_mc_ratio(cutflows_by_sample, outfile: Path, log=True, data_keys=None):

    if data_keys is None:
        raise ValueError("You must specify the keys that correspond to Data samples (e.g., ['EGamma', 'Muon'])")

    first_sample = next(iter(cutflows_by_sample.values()))
    cuts = list(first_sample.keys())
    ncuts = len(cuts)
    x = np.arange(ncuts)

    data_vals = np.zeros(ncuts)
    total_mc_vals = np.zeros(ncuts)
    mc_samples = {}

    for name, cf in cutflows_by_sample.items():
        counts = np.array([cf.get(c, 0.0) for c in cuts])
        if name in data_keys:
            data_vals += counts
        else:
            mc_samples[name] = counts
            total_mc_vals += counts

    if np.sum(data_vals) == 0:
        logger.warning(f"Warning: None of the data keys {data_keys} were found in the input.")

    # Setup Figure (Top panel for bar plot, Bottom panel for ratio)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0.07)

    # --- Top Panel ---
    n_mc = len(mc_samples)
    cmap = plt.get_cmap('viridis')  # Choose your colormap here
    colors = cmap(np.linspace(0, 1, n_mc))

    bottom = np.zeros(ncuts)

    # Stack MC backgrounds
    for i, (name, counts) in enumerate(mc_samples.items()):
        ax1.bar(
            x,
            counts,
            bottom=bottom,
            label=name,
            width=0.8,
            color=colors[i],
            edgecolor='white',
            linewidth=0.5
        )
        bottom += counts

    # Overlay Data
    ax1.errorbar(x, data_vals, yerr=np.sqrt(data_vals), fmt='ko', label='Data', capsize=3, zorder=5)

    ax1.set_ylabel("Events")
    ax1.set_title("Data vs. MC Comparison Cutflow Comparison")

    if log:
        ax1.set_yscale("log")
        # Calculate an appropriate upper limit for log scale
        max_val = max(np.max(total_mc_vals), np.max(data_vals))
        ax1.set_ylim(10, max_val * 100000)

    ax1.legend(fontsize=9, loc='upper right', ncol=2)
    ax1.grid(True, which="both", axis="y", alpha=0.3)

    # --- Bottom Panel ---
    ratio = np.divide(data_vals, total_mc_vals, out=np.zeros_like(data_vals), where=total_mc_vals != 0)

    ratio_err = ratio * np.sqrt(
        np.divide(1, data_vals, out=np.zeros_like(data_vals), where=data_vals != 0) +
        np.divide(1, total_mc_vals, out=np.zeros_like(total_mc_vals), where=total_mc_vals != 0)
    )

    ax2.errorbar(x, ratio, yerr=ratio_err, fmt='ko', capsize=3)
    ax2.axhline(1.0, color='#d63031', linestyle='--', linewidth=1.5)  # Reference line at 1.0

    ax2.set_ylabel("Data / MC", fontsize=10)
    ax2.set_ylim(0.0, 2.0)  # Adjust based on expected agreement
    ax2.set_xticks(x)
    ax2.set_xticklabels(cuts, rotation=45, ha="right")
    ax2.grid(True, axis='y', alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

    logger.info(f"Saved stacked ratio plot to {outfile}")


def plot_individual_efficiencies(cutflows_by_sample, outdir: Path, extension="png"):

    for sample_name, cf in cutflows_by_sample.items():
        cut_names = list(cf.keys())
        values = list(cf.values())

        # 1. Normalize to the first cut
        initial_val = values[0] if values[0] > 0 else 1.0
        efficiencies = [(v / initial_val) * 100 for v in values]

        # 2. Setup the Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        y_pos = np.arange(len(cut_names))

        # We reverse the order so the first cut is at the top
        ax.barh(y_pos, efficiencies, align='center', color='#5790fc', edgecolor='navy', alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(cut_names)
        ax.invert_yaxis()  # top-to-bottom

        # 3. Add text labels for the exact % on each bar
        for i, v in enumerate(efficiencies):
            ax.text(v + 1, i, f"{v:.2f}%", va='center', fontweight='bold', fontsize=9)

        ax.set_xlabel('Efficiency (%)')
        ax.set_title(f'Cutflow Efficiency: {sample_name}')
        ax.set_xlim(0, 115)

        fig.tight_layout()

        # 4. Save with a clean filename
        clean_name = "".join(x for x in sample_name if x.isalnum() or x in "._- ")
        output = outdir / f"efficiency_{clean_name}.{extension}"
        fig.savefig(output)
        plt.close(fig)

    logger.info(f"Generated efficiency plots in {outdir}")

# --------------------------------------------------
# Main
# --------------------------------------------------


def main():
    parser = ArgumentParser(description="Generic cutflow plotter")
    parser.add_argument("json_file", type=Path, help="Path to the cutflows.json")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("plots"), help="Output directory")
    parser.add_argument("--ext", choices=["pdf", "svg", "png"], help="Output file format", default="pdf")
    parser.add_argument("--level", default="all", help="Cutflow level to use (default: all)")
    parser.add_argument("--linear", action="store_true", help="Use linear scale instead of log")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable logging output")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(
            format="[%(name)s(%(levelname)s)] %(message)s",
            level=logging.INFO
        )
    logger.info(f"Reading {args.json_file}")

    with args.json_file.open("r") as f:
        data = json.load(f)

    cutflows_by_sample = OrderedDict()

    for sample, content in data.items():
        cf = extract_cutflow(content, level=args.level)
        if cf:
            cutflows_by_sample[sample] = cf

    logger.info(f"Using {len(cutflows_by_sample)} samples")

    outdir_eff = args.outdir / "efficiencies"
    args.outdir.mkdir(parents=True, exist_ok=True)
    outdir_eff.mkdir(parents=True, exist_ok=True)

    # Cutflow comparison between Data and Stacked MC
    plot_stacked_data_mc_ratio(
        cutflows_by_sample,
        outfile=args.outdir / f"cutflow_ratio.{args.ext}",
        log=not args.linear,
        data_keys=["EGamma", "Muon"]  # Define which samples correspond to the Data streams
    )

    # Cutflow efficiency for each sample
    plot_individual_efficiencies(
        cutflows_by_sample,
        outdir=outdir_eff,
        extension=args.ext
        )


if __name__ == "__main__":
    main()
