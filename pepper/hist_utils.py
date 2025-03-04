import numpy as np
import hist


def group_histogram(h, old_axis, new_axis, mapping):
    """
    Group the categories present in the histogram for a specified axis
    by a mapping provided. Can e.g. be used to group datasets into
    processes.

    This function is similar to the deprecated `hist.group()` function
    from the coffea.hist package.

    Parameters
    ----------
    h
        A histogram to group.
    old_axis
        Name of the axis that should be grouped.
    new_axis
        Name of the new axis that will be created after grouping.
    mapping
        A dictionary mapping categories from the old axis to the new.

    Returns
    ----------
        The grouped histogram.
    """

    new_axes = []

    ax_index = None
    for i, ax in enumerate(h.axes):
        if ax.name == old_axis:
            new_axes.append(hist.axis.StrCategory(categories=mapping.keys(),
                                                  name=new_axis))
            ax_index = i
        else:
            new_axes.append(ax)

    if ax_index is None:
        raise ValueError(f"Axes {old_axis} not found in histogram!")

    h_new = hist.Hist(*new_axes, name=h.name, storage=h._storage_type())
    h_new_view = h_new.view(flow=True)

    categories_in_hist = list(h.axes[old_axis])

    for i, categories_for_group in enumerate(mapping.values()):
        categories_for_group = [cat for cat in categories_for_group
                                if cat in categories_in_hist]
        if len(categories_for_group) > 0:
            h_group = h.integrate(old_axis, categories_for_group)
            index_tuple = [i if j == ax_index else None
                           for j in range(len(new_axes))]
            h_new_view[tuple(index_tuple)] = h_group.view(flow=True)

    return h_new


def rebin_histogram(h, axis_name, edges, new_axis_label=None):
    """
    Rebin a dense axis of a histogram.

    Code was adapted from Kenneth Long, found here:
    https://gist.github.com/kdlong/d697ee691c696724fc656186c25f8814

    Parameters
    ----------
    h
        A histogram.
    axis_name
        Name of the axis that should be rebinned.
    edges
        List or numpy array of the new bin edges
    new_axis_label
        Optional, a new label for the rebinned axis. Default is None,
        i.e. keeping the old label.

    Returns
    ----------
        The rebinned histogram.
    """

    if isinstance(edges, int):
        return h[{axis_name: hist.rebin(edges)}]

    ax = h.axes[axis_name]
    ax_idx = [a.name for a in h.axes].index(axis_name)
    if not all([np.isclose(x, ax.edges).any() for x in edges]):
        raise ValueError(f"Cannot rebin histogram due to incompatible edges "
                         f"for axis '{ax.name}'\n"
                         f"Edges of histogram are {ax.edges}, requested "
                         f"rebinning to {edges}")

    # If you rebin to a subset of initial range, keep the overflow and
    # underflow
    overflow = ax.traits.overflow \
        or (edges[-1] < ax.edges[-1]
            and not np.isclose(edges[-1], ax.edges[-1]))
    underflow = ax.traits.underflow \
        or (edges[0] > ax.edges[0]
            and not np.isclose(edges[0], ax.edges[0]))
    flow = overflow or underflow
    new_axis_label = ax.label if new_axis_label is None else new_axis_label
    new_ax = hist.axis.Variable(edges, name=ax.name, label=new_axis_label,
                                overflow=overflow, underflow=underflow)
    axes = list(h.axes)
    axes[ax_idx] = new_ax

    hnew = hist.Hist(*axes, name=h.name, storage=h._storage_type())

    # Offset from bin edge to avoid numeric issues
    offset = 0.5*np.min(ax.edges[1:]-ax.edges[:-1])
    edges_eval = edges+offset
    edge_idx = ax.index(edges_eval)
    # Avoid going outside the range, reduceat will add the last index anyway
    if edge_idx[-1] == ax.size+ax.traits.overflow:
        edge_idx = edge_idx[:-1]

    if underflow:
        # Only if the original axis had an underflow should you offset
        if ax.traits.underflow:
            edge_idx += 1
        edge_idx = np.insert(edge_idx, 0, 0)

    # Take is used because reduceat sums i:len(array) for the last entry, in
    # the case where the final bin isn't the same between the initial and
    # rebinned histogram, you want to drop this value. Add tolerance of
    # 1/2 min bin width to avoid numeric issues
    hnew.values(flow=flow)[...] = np.add.reduceat(
            h.values(flow=flow),
            edge_idx,
            axis=ax_idx
        ).take(
            indices=range(new_ax.size+underflow+overflow),
            axis=ax_idx
        )
    if hnew._storage_type() == hist.storage.Weight():
        hnew.variances(flow=flow)[...] = np.add.reduceat(
                h.variances(flow=flow),
                edge_idx,
                axis=ax_idx
            ).take(
                indices=range(new_ax.size+underflow+overflow),
                axis=ax_idx
            )
    return hnew


def replace_missing_systematics(h):
    """
    Replaces systematic variations that are zero in all bins by the
    nominal values.

    In pepper output histograms, systematics computed using alternate
    datasets (e.g. top mass, hdamp, tune) are only filled into the
    histograms for the corresponding nominal dataset, while other datasets
    are zero. This needs to be fixed before grouping the histograms since
    some datasets in the group might have this problem while others dont.
    This function takes care of this.

    Note that this operation is in-place, i.e. modifies the original histogram.

    Parameters
    ----------
    h
        A histogram, containing systematics variations in the "sys" axis and
        datasets in the "dataset" axis.

    Returns
    ----------
        The modified histogram.
    """

    if "sys" not in h.axes.name:
        raise ValueError("Histogram has no systematic axis!")

    sys_ax_index = h.axes.name.index("sys")
    nominal_index = h.axes["sys"].index("nominal")

    view = h.view(flow=True)

    other_axes = tuple([i for i, ax in enumerate(h.axes) if ax.name != "sys"
                        and ax.name != "dataset"])

    is_zero = np.all(h.values(flow=True) == 0., axis=other_axes, keepdims=True)

    nominal = np.take(view, nominal_index, axis=sys_ax_index)
    nominal = np.expand_dims(nominal, sys_ax_index)

    np.copyto(view, nominal, where=is_zero)

    return h


def scale_histogram(h, axis_name, scales):
    ax = h.axes[axis_name]
    ax_ind = h.axes.name.index(axis_name)
    for axval, scale in scales.items():
        ind_tuple = tuple([ax.index(axval) if i == ax_ind else slice(None)
                           for i in range(len(h.axes))])
        h.view(flow=True)[ind_tuple] *= scale
