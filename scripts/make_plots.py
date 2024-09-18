"""
Plotting script for pepper histogram output.
Can read both ROOT and hist (.coffea) histograms.

Steered via a plot_config.json file
(see `config_documentation.md` for format.)

While it is quite flexible, it is not meant as a one-fits-all
solution, but instead should be extended and adapted to the needs
of your analysis.

Structure of the code:
- Actual plotting happens in the `plot` function. If you want to
    change the plot style, this is your first step.
- Rebinning & grouping of histograms is done via utility functions
    imported from pepper.hist_utils. Feel free to use those also
    in your scripts.
- File I/O and categorization is done in the latter half of the
    script.

"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep
from argparse import ArgumentParser
import tqdm
import hist
import itertools
import traceback
from warnings import filterwarnings

# The pepper class "HistCollection" is used to load hists.json files
from pepper import Config, HistCollection
# Utility methods for working with `hist` histograms
from pepper.hist_utils import group_histogram, rebin_histogram, \
    replace_missing_systematics

filterwarnings("ignore", "List indexing selection is experimental")
filterwarnings("ignore", "invalid value encountered in divide")
filterwarnings("ignore", "divide by zero encountered in divide")
filterwarnings("ignore", "invalid value encountered in add")
filterwarnings("ignore", "invalid value encountered in double_scalars")

# Use Agg (non-interactive) backend, which plays nicer with headless
# environments
matplotlib.use("Agg")
plt.set_loglevel("error")

# Use the Official CMS Plot Style (TM)
# Requires an up-to-date mplhep version
plt.style.use(mplhep.style.CMS)
cms_color_scheme = plt.rcParams['axes.prop_cycle'].by_key()['color']


def plot(h, config, dense_axis, outfile, log=False):
    """
    Plot a already grouped 1D histogram.

    This is the method you should modify if you want to change the
    plotting style.

    Parameters
    ----------
    h
        A histogram containing all MC samples, possible data and signals,
        and systematic variations if present, already grouped by process.
    config
        The plotting config.
    dense_axis
        The dense axis of `h` that should become the x axis of the plot.
    outfile
        Output file for the plot.
    log
        Whether to plot in log scale.
    """

    # Make a plot with two vertical subplots (total and ratio)
    fig, (ax1, ax2) = plt.subplots(
        nrows=2, sharex=True, gridspec_kw={"height_ratios": [2, 1]})

    # If we have systematics: retrieve the nominal histogram
    if "sys" in h.axes.name:
        h_nominal = h[{"sys": "nominal"}]
    else:
        h_nominal = h

    # x axis edges
    edges = dense_axis.edges

    # Build a stack plot for the background processes.
    # mplhep wants lists of yields, labels and colors for each BG process.
    stack_vals = []
    stack_labels = []
    stack_colors = []
    # Iterate over BG processes defined in the config
    for i, (process, bg_conf) in enumerate(config["backgrounds"].items()):
        # Check that the process is actually in the histogram - if not, ignore
        if process in h_nominal.axes["process"]:
            # Retrieve a hist for the specific process
            h_bg = h_nominal[{"process": process}]
            stack_vals.append(h_bg.values())
            if "color" in bg_conf:
                # If a color is defined in the config: use it
                stack_colors.append(bg_conf["color"])
            else:
                # If no color is defined: use one from the CMS color scheme
                stack_colors.append(
                    cms_color_scheme[i % len(cms_color_scheme)])
            # If no label is defined in the config, use the process key
            stack_labels.append(bg_conf["label"]
                                if "label" in bg_conf else process)

    # Plot the background as a stack plot on the upper panel
    mplhep.histplot(
        stack_vals,
        edges,
        stack=True,
        histtype="fill",
        label=stack_labels,
        color=stack_colors,
        sort='yield',
        ax=ax1)

    # Get summed yields and variances for all background processes
    all_backgrounds = list(config["backgrounds"].keys())
    h_nominal_mc_summed = h_nominal.integrate("process", all_backgrounds)
    mc_vals_summed = h_nominal_mc_summed.values()
    mc_vars_summed = h_nominal_mc_summed.variances()

    if "sys" in h.axes.name:
        # Compute systematic uncertainty using the helper function below
        sys_unc_up, sys_unc_down = \
            compute_systematic_uncertainty(h, config)

        # Sum systematic uncertainty in quadrature with MC stat uncertainty
        total_unc_up = np.sqrt(mc_vars_summed + sys_unc_up**2)
        total_unc_down = np.sqrt(mc_vars_summed + sys_unc_down**2)
    else:
        # No systematics present in histogram - only use MC stat uncertainty
        total_unc_up = np.sqrt(mc_vars_summed)
        total_unc_down = total_unc_up

    # Plot the uncertainty on the upper panel as a grey transparent band
    ax1.stairs(
        values=mc_vals_summed + total_unc_up,
        baseline=mc_vals_summed - total_unc_down,
        edges=edges,
        fill=True,
        facecolor="grey",
        alpha=0.3,
        label="Uncertainty"
    )

    # Uncertainty on the ratio = divide by sum of backgrounds
    relative_unc_up = total_unc_up / mc_vals_summed
    relative_unc_down = total_unc_down / mc_vals_summed
    # Plot the uncertainty on the lower panel as a grey transparent band
    ax2.stairs(
        values=1. + relative_unc_up,
        baseline=1. - relative_unc_down,
        edges=edges,
        fill=True,
        facecolor="grey",
        alpha=0.3,
    )

    # Minimum y axis limit in the ratio plot
    ratio_ylim = 0.01

    # Plot data if it is there in the histogram
    if "data" in h_nominal.axes["process"]:
        data_hist = h_nominal[{"process": "data"}]
        # Data markers in the upper panel
        mplhep.histplot(
            data_hist.values(),
            edges,
            yerr=np.sqrt(data_hist.variances()),
            histtype="errorbar",
            color="black",
            label="Data",
            markersize=12,
            ax=ax1
        )

        ratio_hist = data_hist / mc_vals_summed

        # Data markers in the lower panel
        mplhep.histplot(
            ratio_hist.values(),
            edges,
            yerr=np.sqrt(ratio_hist.variances()),
            histtype="errorbar",
            color="black",
            markersize=12,
            ax=ax2
        )

        # Increase the y axis limit in the ratio plot if necessary
        # so that all data markers are visible
        ratio_vals = np.nan_to_num(ratio_hist.values(),
                                   nan=1., posinf=1., neginf=1.)
        ratio_ylim = max(ratio_ylim, np.amax(abs(ratio_vals-1.)))

    # Plot signals if they are defined in the config
    if "signals" in config:
        # Iterate over all signals that are defined
        for process, sig_conf in config["signals"].items():
            # If the signal is not in the histogram, ignore it
            if process in h_nominal.axes["process"]:
                sig_hist = h_nominal[{"process": process}]

                # Get color and label fron config, use defaults if
                # not defined
                sig_color = sig_conf["color"] \
                    if "color" in sig_conf else None
                sig_label = sig_conf["label"] \
                    if "label" in sig_conf else process

                sig_vals = sig_hist.values()

                # If the flag to stack signals on top of the background is
                # given, do exactly that
                if "stack_signals" in config and config["stack_signals"]:
                    sig_vals = sig_vals + mc_vals_summed

                # Plot the signal in the upper panel
                mplhep.histplot(
                    sig_vals,
                    edges,
                    yerr=np.sqrt(sig_hist.variances()),
                    histtype="step",
                    edges=False,
                    label=sig_label,
                    color=sig_color,
                    linewidth=1,
                    ax=ax1
                )

                sig_ratio_hist = sig_hist / mc_vals_summed

                # Plot the signal in the lower panel
                # Here we always want the ratio (sig + bg) / bg
                mplhep.histplot(
                    sig_ratio_hist.values() + 1.,
                    edges,
                    yerr=np.sqrt(sig_ratio_hist.variances()),
                    histtype="step",
                    edges=False,
                    color=sig_color,
                    linewidth=1,
                    ax=ax2
                )

                # Increase the y axis limit in the ratio plot if
                # necessary so that the signal is fully visible
                ratio_vals = np.nan_to_num(sig_ratio_hist.values(),
                                           nan=0., posinf=0., neginf=0.)
                ratio_ylim = max(ratio_ylim, np.amax(abs(ratio_vals)))

    # Black dashed line at 1 for the ratio plot
    ax2.hlines(1, edges[0], edges[-1], color="black", linestyle="dashed",
               linewidth=1., zorder=-90)

    # Axis formatting for the upper panel
    ax1.set_xlabel("")
    ax1.set_ylabel("Event yield")
    # Make sure to use scientific notation
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0),
                         useMathText=True)
    # We need to push the "* 10^X" to the left because otherwise it
    # overlaps with the CMS logo -.-
    ax1.yaxis.get_offset_text().set_x(-0.09)

    # Log scale if requested
    if log:
        ax1.set_yscale("log")
        ax1.autoscale("y")

    # Axis formatting for the lower panel
    ax2.set_ylabel("Data / Pred.")
    ax2.set_xlim(edges[0], edges[-1])
    ax2.set_ylim(1.-ratio_ylim*1.05, 1.+ratio_ylim*1.05)
    ax2.set_xlabel(dense_axis.label)

    # Plot legend (with two columns)
    ax1.legend(ncol=2)

    # Add CMS logo, year and center-of-mass energy if given in config
    if "year" in config:
        year = config["year"]
        luminosity = config["luminosity"]
        com_energy = config["com_energy"] if "com_energy" in config else 13
        mplhep.cms.label(ax=ax1, data=True, year=year,
                         lumi=f"{luminosity:.1f}", com=com_energy)

    # Generic figure adjustments
    fig.align_ylabels()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0)

    # Save the plot to file
    fig.savefig(outfile)
    plt.close()


def compute_systematic_uncertainty(h, config):
    """
    Computes the total systematic uncertainty (up and down )for a histogram
    by taking the minimum & maximum for different sources, and then summing
    different sources in quadrature.

    Change this method if you want to treat uncertainties differently.

    Parameters
    ----------
    h
        A histogram containing all MC samples and systematic variations,
        already grouped by process.
    config
        The plotting config.

    Returns
    -------
        Tuple (up, down) containing the total systematic uncertainty.
    """

    all_backgrounds = list(config["backgrounds"].keys())

    # Integrate over all BG samples
    h_summed = h.integrate("process", all_backgrounds)

    # Nominal values
    values_nominal = h_summed[{"sys": "nominal"}].values()

    # Dictionary that will store the maximum and minimum deviations
    # to the nominal for each uncertainty source
    sys_diff_dict = {}

    for sys in h_summed.axes["sys"]:
        if sys == "nominal":
            continue

        # Get the uncertainty source from the variation name
        # E.g. "hdamp_down" & "hdamp_up" belong to source "hdamp"
        # "PDF_10" belongs to source "PDF"
        sys_name = sys
        if "_" in sys_name:
            sys_dir = sys_name.split("_")[-1]
            # If the variation name ends neither in "up", "down", or a
            # number, treat it as its own source (i.e. an one-sided unc)
            if sys_dir == "up" or sys_dir == "down" or sys_dir.isdigit():
                sys_name = "_".join(sys_name.split("_")[:-1])

        values_sys = h_summed[{"sys": sys}].values()
        diff = values_nominal - values_sys

        # We want to take the minimum and maximum of the deviations from
        # the nominal for a given source
        if sys_name in sys_diff_dict:
            # If other variations for the same source are already
            # in the dict, retrieve them
            previous_max, previous_min = sys_diff_dict[sys_name]
        else:
            # Otherwise start at zero deviation
            previous_max = 0.
            previous_min = 0.

        # Take the minimum and maximum with what is already there
        new_max = np.maximum(diff, previous_max)
        new_min = np.minimum(diff, previous_min)
        # Save it in the dictionary
        sys_diff_dict[sys_name] = (new_max, new_min)

    # Sum the different sources (= entries in the dict) in quadrature
    all_diffs = np.array(list(sys_diff_dict.values()))
    total_unc = np.sqrt(np.sum(all_diffs**2, axis=0))
    # tuple structure is (maximum, minimum) --> indices 0 and 1
    total_up = total_unc[0]
    total_down = total_unc[1]
    # return up & down uncertainty as a tuple
    return total_up, total_down


def get_category_dict(h, config):
    """
    Build a dictionary of containing definitions of all categories
    we want to plot for a specific histogram, as specified in the
    plot config.

    Modify this if you want different category behaviour.

    Parameters
    ----------
    h
        The histogram to categorize, already grouped by process.
    config
        The plotting config.

    Returns
    -------
        A dictionary containing the category defintions.
    """

    # List of coarse (i.e. category) axes in the histograms
    # except for systematics and processes, we treat those extra
    coarse_axes = [
        ax
        for ax in h.axes
        if (isinstance(ax, hist.axis.StrCategory)
            or isinstance(ax, hist.axis.IntCategory))
        and ax.name != "sys" and ax.name != "process"
    ]

    # Make a dict which, for each category we want to plot, contains a
    # list of which axis values should be included in this category
    category_dict = {}

    # If 'sum_categories' is true in the config, simply make one category
    # called "sum" which sums over all coarse axes in the histogram
    if "sum_categories" in config and config["sum_categories"]:
        category_dict["sum"] = {
            ax.name: list(ax) for ax in coarse_axes
        }

    # If categories are given in the config explicitly, use those
    elif "categories" in config:
        # Iterate over the categories defined in the config
        for category_label, category_definition \
                in config["categories"].items():

            # Check that at least one of the axes in the histogram
            # has a definition in the category. Otherwise, skip the
            # category
            if any(axname not in h.axes.name
                   for axname in category_definition.keys()):
                print(f"Warning: None of the coarse axes specified for "
                      f"category {category_label} was found in histogram. "
                      f"Ignoring this category")
                continue

            # Copy the category definition for modification
            category_definition_modified = category_definition.copy()

            # Check for axes that are present in the histogram but missing
            # in the category definition in the config
            for ax in coarse_axes:
                if ax.name not in category_definition:
                    # Axis is present in histogram but missing in config
                    # --> sum over this axis
                    category_definition_modified[ax.name] = list(ax)
                    print(f"Warning: Coarse axis {ax.name} found in histogram "
                          f"but not in definition for category "
                          f"{category_label}. Will sum over this axis")

            # Add to the dictionary of categories
            category_dict[category_label] = category_definition_modified

        # If no categories are valid for this histogram (e.g. because it
        # is early in the cutflow), fall back to plotting the sum
        if len(category_dict) == 0:
            category_dict["sum"] = {
                ax.name: list(ax) for ax in coarse_axes
            }

    # If nothing is defined in the config, build a dictionary with all possible
    # combinations of values present in the coarse axes of the histogram
    else:
        # For each coarse_axis, use all values in the histogram, as well as
        # a sum over that axis
        all_axis_values = [
            [*ax, "sum"] for ax in coarse_axes
        ]
        # Make a cartesian product using itertools
        for axis_values in itertools.product(*all_axis_values):
            # Build a dictionary for this category and a label
            category_definition = {}
            label_parts = []
            for ax, axis_value in zip(coarse_axes, axis_values):
                if axis_value == "sum":
                    category_definition[ax.name] = list(ax)
                else:
                    category_definition[ax.name] = axis_value
                label_parts.append(ax.name + "_" + axis_value)
            category_label = "_".join(label_parts)
            category_dict[category_label] = category_definition

    return category_dict


# Set up command line arguments
parser = ArgumentParser(
    description="Plot histograms from previously created histograms")
parser.add_argument(
    "plot_config", help="Path to a configuration file for plotting")
parser.add_argument(
    "histfile", help="Coffea file with a single histogram or a "
    "JSON file containing histogram info. See output of select_events.py")
parser.add_argument(
    "-o", "--outdir", help="Output directory. If not given, output to the "
    "directory where histfile is located")
parser.add_argument(
    "--log", action="store_true", help="Make logarithmic plots")
parser.add_argument(
    "--ext", choices=["pdf", "svg", "png"], help="Output file format",
    default="pdf")
parser.add_argument(
    "-c", "--cut", type=str, default=None, help="If specified, only "
    "plot a given cut (i.e. only 'Req MET'.)")
args = parser.parse_args()

# Read the config from json using the pepper Config class
config = Config(args.plot_config)

# From the config, build a dictionary for grouping the datasets in
# the histogram into processes
# This includes BG processes, signal processes, and observed data
groups = {}
for process, bg in config["backgrounds"].items():
    groups[process] = bg["datasets"]
if "signals" in config:
    for process, sig in config["signals"].items():
        groups[process] = sig["datasets"]
if "data" in config:
    groups["data"] = config["data"]

# Do we plot a collection of hists using a hists.json file, or
# just a single histogram ?
is_hist_json = args.histfile.endswith(".json")
if is_hist_json:
    # Read the info about all histograms from the hists.json
    # by constructing a HistCollection object
    with open(args.histfile) as f:
        hists = HistCollection.from_json(f)
else:
    # Create a HistCollection for only the single histogram
    hists = HistCollection.from_single_hist(args.histfile)

# Iterate over all histograms in our collection
for histkey in tqdm.tqdm(hists.keys()):

    # Retrieve the cut and the variable for the histogram
    if is_hist_json:
        # For hists.json, this is stored in the histkey
        cutname = histkey[0]
        variable = histkey[1]

    else:
        # Try to get cut and variable from the hist file name
        filename = os.path.basename(args.histfile)
        filename = filename.split(".")[0]

        if filename.startswith("Cut "):
            # Old file name syntax (with spaces)
            filename = " ".join(filename.split(" ")[2:])
            filename_parts = filename.split("_")
            cutname = filename_parts[0]
            variable = "_".join(filename_parts[1:])
        elif filename.startswith("Cut_"):
            # New file name syntax (with underscores)
            filename_parts = filename.split("_")
            cutname = filename_parts[2]
            variable = "_".join(filename_parts[3:])
        else:
            # Otherwise give up
            cutname = ""
            variable = ""

    # If "--cut" is given as a command line argument, and
    # the histogram is from a different cut, skip
    if args.cut is not None:
        if cutname != args.cut:
            continue

    print(f"Processing file {hists[histkey]}...")

    # Actually load the histogram into memory
    # This will internally convert ROOT histograms to hist
    h = hists.load(histkey)

    # Only if we have systematic uncertainties:
    # Sometimes certain systematics are missing for some datasets,
    # e.g. for systematics that come from alternate MC samples,
    # which will only be defined for that sample.
    # We need to fix this before grouping the datasets into processes.
    # The utility method called here simply replaces systematics
    # that are 0 everywhere by the nominal
    if "sys" in h.axes.name:
        h = replace_missing_systematics(h)

    # Group the histogram by datasets, using the dictionary defined
    # above from the config.
    # This will create a new coarse axis in the histogram, which we
    # call "process".
    h = group_histogram(h, old_axis="dataset", new_axis="process",
                        mapping=groups)

    # Get a list of dense axes in the histogram
    dense_axes = [
        ax
        for ax in h.axes
        if not (isinstance(ax, hist.axis.StrCategory)
                or isinstance(ax, hist.axis.IntCategory))
    ]

    # For now, we support only one dense axis per histogram
    if len(dense_axes) != 1:
        print(f"Skipping histogram {variable} for cut {cutname} "
              f"with {len(dense_axes)} dense axes")
        continue
    # Get the name of our single dense axis
    dense_axis_name = dense_axes[0].name

    # --- Rebinning ---
    if "rebin" in config and variable in config["rebin"]:
        # Rebin histogram as given in the config for this variable
        rebin_config = config["rebin"][variable]
        if "n_or_arr" not in rebin_config:
            # Invalid rebin config, raise error
            raise ValueError(f"Invalid rebin config for variable {variable}")
        if "lo" in rebin_config and "hi" in rebin_config:
            # Rebinning using a regularly spaced axis
            new_binning = np.linspace(rebin_config["lo"], rebin_config["hi"],
                                      rebin_config["n_or_arr"]+1)
        else:
            # Rebinning using a variable axis
            new_binning = np.array(rebin_config["n_or_arr"])
        # Use new label for the axis if given in config
        new_axis_label = rebin_config["label"] \
            if "label" in rebin_config else None
        try:
            # Rebin using a utility method from pepper.hist_utils
            h = rebin_histogram(h, dense_axis_name, new_binning,
                                new_axis_label)
        except ValueError:
            print(f"Could not rebin histogram {variable} for cut {cutname}")
            traceback.print_exc()
            continue

    # Get a dictionary of all categories we want for this histogram
    # using the helper method defined above
    category_dict = get_category_dict(h, config)

    # Iterate over all combinations of categories defined earlier
    for category_label, category_definition in category_dict.items():

        # Make a subfolder for this cut and category
        
        if args.outdir is not None:
            cat_outfolder = os.path.join(args.outdir, cutname, category_label)
        else:
            cat_outfolder = os.path.join(os.path.dirname(args.histfile), cutname, category_label)
        os.makedirs(cat_outfolder, exist_ok=True)

        # Build a output file name for the plot
        outfile = f"{variable}_{cutname}_{category_label}.{args.ext}"
        outfile = os.path.join(cat_outfolder, outfile)

        # Slice the histogram according to the category definition
        h_cat = h
        for axname, axvals in category_definition.items():
            # Make sure every entry is a list, even if it is just one
            # Otherwise hist.integrate will give nonsense
            if not isinstance(axvals, list):
                axvals = [axvals]
            h_cat = h_cat.integrate(axname, axvals)

        try:
            # Plot with the method defined above
            plot(h_cat, config, h_cat.axes[dense_axis_name], outfile,
                 log=args.log)
        except ValueError:
            print(f"Error for variable {variable}, cutname {cutname}, "
                  f"category {category_label}:")
            traceback.print_exc()
