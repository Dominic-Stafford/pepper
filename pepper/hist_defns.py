import hist as hi
import awkward as ak
import numpy as np
from collections import defaultdict
import functools
import itertools


def not_arr(arr):
    """Element wise logical not"""
    return ~arr


def equal(arr, other):
    """Element wise equal, operats can both be arrays or scalars"""
    return arr == other


def leaddiff(quantity):
    """Difference in ``quantity`` of the two leading particles (at positiion
    [:, 0] and [:, 1])"""
    if isinstance(quantity, np.ndarray):
        if quantity.ndim == 1:
            quantity = quantity.reshape((-1, 1))
        quantity = ak.Array(quantity)
    enough = quantity[ak.broadcast_arrays(ak.num(quantity) >= 2, quantity)[0]]
    return enough[:, 0:1] - enough[:, 1:2]


def concatenate(*arr, axis=0):
    """Concatenate arrays along an axis"""
    arr_processed = []
    for a in arr:
        if isinstance(a, list):
            a = np.array(a)
        if isinstance(a, np.ndarray):
            if a.ndim == 1:
                a = a.reshape((-1, 1))
            a = ak.Array(a)
        arr_processed.append(a)
    return ak.concatenate(arr_processed, axis=axis)


func_dict = {
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "arcsin": np.arcsin,
    "arccos": np.arccos,
    "arctan": np.arctan,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
    "arcsinh": np.arcsinh,
    "arccosh": np.arccosh,
    "arctanh": np.arctanh,
    "exp": np.exp,
    "log": np.log,
    "log10": np.log10,
    "sqrt": np.sqrt,
    "abs": np.abs,
    "sign": np.sign,

    "not": not_arr,
    "equal": equal,
    "leaddiff": leaddiff,
    "concatenate": concatenate,

    "sum": lambda x: ak.sum(x, axis=1),
    "num": ak.num,
    "zeros_like": ak.zeros_like,
    "ones_like": ak.ones_like,
    "full_like": ak.full_like,
}


class HistDefinitionError(Exception):
    pass


class HistFillError(ValueError):
    pass


class HistDefinition:
    """Holds a definition of a histogram and creates histograms from it

    The definition is parsed from a dictionary that can easily be saved to
    JSON. The created histogram can have any number of additional
    categorizations (``StrCategory`` axes) on top of what is already in the
    definition."""

    def __init__(self, config):
        """
        Parameters
        ----------
        config
            The definition of the histogram in form of a dict.
            Please see the ``config_documentation.md`` of Pepper
        """
        self._label = config.get("label", None)
        self.dataset_axis = hi.axis.StrCategory(
            [], name="dataset", label="Dataset name", growth=True)
        bins = []
        if "bins" in config:
            for bin_config in config["bins"]:
                unit = bin_config.pop("unit", None)
                if unit is not None:
                    bin_config["label"] = \
                        bin_config.get("label", "") + f" ({unit})"
                hist_type = bin_config.pop("type", "float")
                if "lo" in bin_config:
                    if hist_type == "int":
                        bin_config.pop("n_or_arr")
                        bin = hi.axis.Integer(bin_config.pop("lo"),
                                              bin_config.pop("hi"),
                                              **bin_config)
                    elif hist_type == "float":
                        bin = hi.axis.Regular(bin_config.pop("n_or_arr"),
                                              bin_config.pop("lo"),
                                              bin_config.pop("hi"),
                                              **bin_config)
                    else:
                        raise HistDefinitionError(
                            f"Unknown histogram type: {hist_type}")
                else:
                    bin = hi.axis.Variable(bin_config.pop("n_or_arr"),
                                           **bin_config)
                bin.unit = unit
                bins.append(bin)
        cats = []
        if "cats" in config:
            for kwargs in config["cats"]:
                default_key = kwargs.pop("default_key", None)
                cat = hi.axis.StrCategory([], growth=True, **kwargs)
                cat.default_key = default_key
                cats.append(cat)
        self.axes = bins + cats
        self.bin_fills = {}
        self.cat_fills = {}
        if "weight" in config:
            self.weight = config["weight"]
        else:
            self.weight = None
        for axisname, method in config["fill"].items():
            if axisname in [b.name for b in bins]:
                self.bin_fills[axisname] = method
            elif axisname in [c.name for c in cats]:
                self.cat_fills[axisname] = method
            else:
                raise HistDefinitionError(
                    f"Fill provided for non-existing axis: {axisname}")
        missing_fills = (set(a.name for a in self.axes)
                         - set(config["fill"].keys()))
        if len(missing_fills) != 0:
            raise HistDefinitionError(
                "Missing fills for axes: " + ", ".join(missing_fills))
        if "step_requirement" in config:
            self.step_requirement = config["step_requirement"]
        else:
            self.step_requirement = None
        if "do_systs" in config:
            self.do_systs = config["do_systs"]
        else:
            self.do_systs = True

    @staticmethod
    def _prepare_fills(fill_vals):
        """Checks for length consistency across the fill_vals,
        removes events where counts do not agree (in case of 2 dims)
        and flattens everything into numpy arrays.

        Parameters
        ----------
        fill_vals
            Dict of awkward arrays, no more than 2 dimensions

        Returns
        -------
        prepared
            Dict of the prepared values, which can be used to fill a Hist
        """

        size = len(next(iter(fill_vals.values())))
        mask = np.full(size, True)  # Mask is purely 1D, i.e. event-level
        counts = None
        prepared = {}
        for key, data in fill_vals.items():
            if data is None:
                # Early return, no need to bother with other fill_vals
                return {key: None}
            if data.ndim > 2:
                raise ValueError(f"Got more than 2 dimensions for {key}")
            if size != len(data):
                raise ValueError(f"Got inconsistant filling size ({size} and"
                                 f"{len(data)} for {key})")
            # Make sure all fills have the same mask originating from ak.mask
            mask = mask & ~np.asarray(ak.is_none(data))

            # Make sure all counts agree
            if data.ndim == 2:
                # Apply ak.fill_none to catch mixed 1D-2D arrays where
                # some entries are None at event level
                num_data = np.asarray(ak.fill_none(ak.num(data), 0))
                if counts is None:
                    counts = num_data
                else:
                    mask = mask & (counts == num_data)
        for key, data in fill_vals.items():
            if counts is not None and data.ndim == 1:
                # Manual broadcast using np.repeat instead of awkward
                # for performance reasons
                data = np.repeat(np.asarray(data[mask]), counts[mask])
                prepared[key] = data
            else:
                prepared[key] = np.asarray(ak.flatten(data[mask], axis=None))
        return prepared

    def create_hist(self, categorizations, has_systematic=False):
        """Create a histogram according to the definition of this instance

        Parameters
        ----------
        categorizations
            Add ``StrCategory`` axes to the histogram. The axes will be named
            and labeled according to the keys of this dict. An exception is a
            category named "dataset", which will be used to overwrite the
            default dataset name.
        has_systematic
            If true, add an ``StrCategory`` axis called 'sys' for systematic
            uncertainties

        Returns
        -------
        hist
            A new, empty histogram according to the definition
        """
        axes = self.axes.copy()
        for cat in categorizations.keys():
            if cat != "dataset":
                axes.append(
                    hi.axis.StrCategory([], name=cat, label=cat, growth=True))
        if has_systematic:
            axes.insert(0, hi.axis.StrCategory(
                [], name="sys", label="Systematic", growth=True))
        hist = hi.Hist(self.dataset_axis, *axes, storage="Weight")
        hist.label = self.label
        return hist

    def __call__(
            self, data, categorizations, dsname, is_mc, weight):
        """Create and fill a histogram according to the data given.

        Parameters
        ----------
        data
            Record array, usually the full NanoEvents, from which the data
            that goes into the histogram is picked
        categorizations
            Dict of categories to add to the histogram on top of what is
            already in the definition. Its keys name the axes, while its values
            are lists of strings that name the inidivuals bins in the axes and
            the fields inside ``data`` to take masks from if an event belongs
            to a category. If a category called "dataset" is supplied, this
            will be used instead of the dsname.
        dsname
            Name of the data set from where the event data is
        is_mc
            Whether this is simulation
        weight
            Event weight as array or dict of arrays (for systematics)
        """

        has_systematic = self.do_systs and self.weight is None \
            and isinstance(weight, dict)

        if not isinstance(weight, dict):
            weight = {None: weight}
        elif not has_systematic:
            if "nominal" in weight:
                weight = {None: weight["nominal"]}
            else:
                raise HistFillError("Hist should not include systs, but "
                                    "no nominal was found")

        hist = self.create_hist(categorizations, has_systematic)

        # Get the fills for the axes from the event-level data
        fill_vals = {name: DataPicker(method)(data)
                     for name, method in self.bin_fills.items()}
        # If any fills are None, we cannot fill the histogram
        if any(val is None for val in fill_vals.values()):
            none_keys = [k for k, v in fill_vals.items() if v is None]
            raise HistFillError(f"No fill for axes: {', '.join(none_keys)}")

        # Option to give an explicit weight in the config
        if self.weight is not None:
            new_weight = DataPicker(self.weight)(data)
            if new_weight is None:
                raise HistFillError(
                    "Weight specified in hist config not available")
            weight = {None: new_weight}

        # Add categories explicitly defined in the histogram config
        # to those already present from the processor
        categorizations = {name: {cat: [cat] for cat in cats}
                           for name, cats in categorizations.items()}
        categorizations.update(self.cat_fills)
        cat_present = defaultdict(dict)

        # Get the bitmask for the different categories
        for name, val in categorizations.items():
            for cat, method in val.items():
                cat_present[name][cat] = DataPicker(method)(data)
        for name, val in cat_present.items():
            if any(mask is None for mask in val.values()):
                # If any category is undefined in the data, try to
                # fill everything in the default key for this category
                ax = next(ax for ax in self.axes if ax.name == name)
                if ax.default_key is None:
                    none_keys = [k for k, v in val.items() if v is None]
                    raise HistFillError(
                        f"No fill for axes: {', '.join(none_keys)}")
                else:
                    key = ax.default_key
                    cat_present[name] = {key: np.full(len(data), True)}
        if "dataset" in categorizations:
            non_array_fills = {}
        else:
            non_array_fills = {"dataset": dsname}

        # Add the weights to the fill_values to process
        for weightname, w in weight.items():
            if w is not None:
                # Use names starting with '__weight' for the weights
                # to not have overlap with the bin fills
                wname = "__weight/" + weightname \
                    if weightname is not None else "__weight"
                fill_vals[wname] = ak.fill_none(w, 0.)

        cat_keys = []
        # Add the category bitmasks to the fill_values
        # and get all present keys for each category
        for cat_name, cat_dict in cat_present.items():
            cat_keys.append(cat_dict.keys())
            for cat_key, cat_mask in cat_dict.items():
                cname = "__cat/" + cat_name + "/" + cat_key
                fill_vals[cname] = cat_mask

        # Actually do the processing & broadcasting, for fills, weights
        # and category bitmasks at once
        prepared = self._prepare_fills(fill_vals)
        # Filter out the category bitmasks again
        prepared_nocats = {k: v for k, v in prepared.items()
                           if not k.startswith("__cat")}

        # Only fill if there are no Nones
        if all(v is not None for k, v in prepared_nocats.items()
               if not k.startswith("__weight")):
            # Iterate over all possible combinations of categories
            for cat_combination in itertools.product(*cat_keys):
                if len(cat_combination) == 0:
                    # No categories - no masking required
                    prepared_masked = prepared_nocats
                    cat_fills = {}
                else:
                    # Make a dictionary for this combination to pass to
                    # hist.fill() later
                    cat_fills = dict(zip(cat_present.keys(), cat_combination))
                    # Get the masks that define this particular combination
                    prepared_masks = [
                        prepared["__cat/" + cat_name + "/" + cat_key]
                        for cat_name, cat_key in cat_fills.items()
                        ]

                    # Combine the masks into one with bitwise and
                    comb_mask = functools.reduce(
                        lambda a, b: a & b, prepared_masks)

                    # Apply the bitmask to the fill values
                    prepared_masked = {}
                    for k, v in prepared_nocats.items():
                        if v is not None:
                            v = v[comb_mask]
                        if v is not None and len(v) > 0:
                            prepared_masked[k] = v
                        else:
                            # Workaround when there are no events:
                            # put zero everywhere including weights
                            # so that the histogram is filled anyway
                            prepared_masked[k] = 0

                # Filter out the weights again for the pass to hist.fill()
                prepared_noweights = {k: v for k, v in prepared_masked.items()
                                      if not k.startswith("__weight")}

                # Iterate over all weights (i.e. systematics) to fill
                for weightname, w in weight.items():
                    if weightname is not None:
                        non_array_fills["sys"] = weightname
                    if w is not None:
                        # Get the prepared and masked weight array
                        wname = "__weight/" + weightname \
                            if weightname is not None else "__weight"
                        prepared_weight = prepared_masked[wname]
                        if prepared_weight is not None:
                            hist.fill(**non_array_fills, **prepared_noweights,
                                      **cat_fills, weight=prepared_weight)
                    else:
                        # No weight array present - fill without weights
                        hist.fill(**non_array_fills, **prepared_noweights,
                                  **cat_fills)

        return hist

    @property
    def label(self):
        """y axis label of the histogram in a format following the CMS
        style guidelines as close as possible"""
        # Make this a property so that if the axes change, label is updated
        if self._label is not None:
            return self._label
        bins = [ax for ax in self.axes if (isinstance(ax, hi.axis.Regular)
                or isinstance(ax, hi.axis.Variable)
                or isinstance(ax, hi.axis.Integer))]
        if len(bins) == 1 and not isinstance(bins[0], hi.axis.Variable):
            # This is a 1d, fixed-bin-width histogram.
            edges = bins[0].edges
            width = (edges[-1] - edges[0]) / bins[0].size
            unit = bins[0].unit
            if unit is None:
                unit = "unit" if width == 1 else "units"
            width = f"{width:.3g}"
            if "e" in width:
                width = width.replace("e", r"\times 10^{") + "}"
            return f"Events / ${width}$ {unit}"
        else:
            return "Events / bin"

    @label.setter
    def label(self, value):
        self._label = value


class DataPicker:
    """Extract and perform operations on specific data inside an Awkward array
    using only a very basic dict that can be specified in JSON"""

    def __init__(self, method):
        """
        Parameters
        ----------
        method
            Defines how to pick the data from an array and what operations to
            perform. Please see the ``config_documentation.md`` of Pepper
        """
        self._method = method

    def __call__(self, data):
        """Perform the picking on the data
        Parameters
        ----------
        data
            The array from which data is extracted

        Returns
        -------
        data
            The result of the data picking operation specified in the data
            picker using the data given
        """
        method = self._method
        orig_data = data
        if not isinstance(method, list):
            raise HistDefinitionError("Fill method must be list")
        for i, sel in enumerate(method):
            if isinstance(sel, str):
                try:
                    data = getattr(data, sel)
                except AttributeError:
                    try:
                        if (isinstance(data, ak.Array)
                                and sel not in ak.fields(data)):
                            # Workaround for awkward raising a ValueError
                            # instead of a KeyError
                            raise KeyError
                        else:
                            data = data[sel]
                    except KeyError:
                        break
                if callable(data):
                    data = data()
                if data is None:
                    break
            elif isinstance(sel, list):
                try:
                    data = data[sel]
                except (ValueError, KeyError):
                    try:
                        data_proc = {}
                        counts = None
                        can_zip = True
                        for attr in sel:
                            data_proc[attr] = getattr(data, attr)
                            if (can_zip
                                    and isinstance(data_proc[attr], ak.Array)
                                    and data_proc[attr].ndim > 1):
                                if counts is None:
                                    counts = ak.num(data_proc[attr])
                                else:
                                    can_zip = ak.all(
                                        counts == ak.num(data_proc[attr]))
                            else:
                                can_zip = False
                        if can_zip:
                            data = ak.zip(data_proc)
                        else:
                            data = ak.Array(data_proc)
                    except AttributeError:
                        break
            elif isinstance(sel, dict):
                if "function" in sel:
                    data = self._pick_data_from_function(
                        i, sel, data, orig_data)
                    if data is None:
                        break
                elif "key" in sel:
                    try:
                        data = data[sel["key"]]
                    except KeyError:
                        break
                elif "attribute" in sel:
                    try:
                        data = getattr(data, sel["attribute"])
                    except AttributeError:
                        break

                if "leading" in sel:
                    leading = sel["leading"]
                    if isinstance(leading, (list, tuple)):
                        start = leading[0] - 1
                        end = leading[1] - 1
                        if start >= end:
                            raise HistDefinitionError(
                                "First entry in 'leading' must be larger than "
                                "the second one")
                    else:
                        start = leading - 1
                        end = start + 1
                    if isinstance(data, np.ndarray):
                        data = ak.Array(data)
                    elif not isinstance(data, ak.Array):
                        raise HistDefinitionError(
                            "Can only do 'leading' on numpy or awkward arrays")
                    data = data[:, start:end]
            else:
                raise HistDefinitionError("Fill constains invalid type, must "
                                          f"be str, list or dict: {sel}")
        else:
            return data
        return None

    def _pick_data_from_function(self, i, sel, data, orig_data):
        funcname = sel["function"]
        if funcname not in func_dict.keys():
            raise HistDefinitionError(f"Unknown function {funcname}")
        func = func_dict[funcname]
        args = []
        kwargs = {}
        if "args" in sel:
            if not isinstance(sel["args"], list):
                raise HistDefinitionError("args for function must be list")
            for value in sel["args"]:
                if isinstance(value, list):
                    args.append(self.__class__(value)(orig_data))
                    if args[-1] is None:
                        return None
                else:
                    args.append(value)
        if "kwargs" in sel:
            if not isinstance(sel["kwargs"], dict):
                raise HistDefinitionError("kwargs for function must be dict")
            for key, value in sel["kwargs"].items():
                if isinstance(value, list):
                    kwargs[key] = self.__class__(value)(orig_data)
                    if kwargs[key] is None:
                        return None
                else:
                    kwargs[key] = value
        if i == 0:
            return func(*args, **kwargs)
        else:
            return func(data, *args, **kwargs)

    @property
    def name(self):
        """Principle string identifying the data picker"""
        name = ""
        for sel in self._method:
            if isinstance(sel, str):
                name += sel.replace("/", "")
            elif isinstance(sel, list):
                # Output will contain multiple attributes of an object, just
                # give the name of this object
                break
            elif isinstance(sel, dict):
                if "function" in sel:
                    name += sel["function"].replace("/", "")
                elif "key" in sel:
                    name += sel["key"].replace("/", "")
                elif "attribute" in sel:
                    name += sel["attribute"].replace("/", "")
                elif "leading" in sel:
                    name += "leading"
                    if isinstance(sel["leading"], tuple):
                        if (not isinstance(sel["leading"][0], int)
                                or not isinstance(sel["leading"][1], int)):
                            raise HistDefinitionError(
                                "Invalid value for leading tuple. Must be int")
                        name += (str(sel["leading"][0]) + "-"
                                 + str(sel["leading"][1]))
                    elif isinstance(sel["leading"], int):
                        name += str(sel["leading"])
                    else:
                        raise HistDefinitionError(
                            "Invalid value for leading. Must be int or tuple")
                else:
                    raise HistDefinitionError(
                        f"Dict without any recognized keys: {sel}")
            else:
                raise HistDefinitionError("Selection constains invalid type, "
                                          f"must be str, list or dict: {sel}")
            name += "_"
        name = name[:-1]  # Remove last underscore
        return name
