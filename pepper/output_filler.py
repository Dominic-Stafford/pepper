import logging
import numpy as np
import awkward as ak
import hist as hi
import itertools
from collections import defaultdict
from collections.abc import MutableMapping
import copy

import pepper


logger = logging.getLogger(__name__)


class DummyOutputFiller:
    """An output filler that does not do anything. Useful if no output filling
    should be performed"""
    def __init__(self, output):
        self.output = output

    def get_callbacks(self):
        return []


class AddableDict(dict):
    """A dictionary add can be added to other dictionaries"""
    # Workaround for Coffea shuffling dict keys when accumulating
    def __add__(self, other):
        if not isinstance(other, MutableMapping):
            raise ValueError("Cannot add object of incompatible type to dict: "
                             f"{type(other)}")
        out = copy.copy(self)
        out.clear()
        lhs, rhs = set(self), set(other)
        for key in self:
            if key in rhs:
                out[key] = self[key] + other[key]
            else:
                out[key] = copy.deepcopy(self[key])
        for key in other:
            if key not in lhs:
                out[key] = copy.deepcopy(other[key])
        return out


class OutputFiller:
    """Fills histograms and cutflows

    Attributes
    ----------
    output
        Dict holding all the cutflows and histograms filled by the filler
    """
    def __init__(self, hist_dict, is_mc, dsname, dsname_in_hist, sys_enabled,
                 sys_overwrite=None, cuts_to_histogram=None,
                 systs_to_histogram=None):
        """
        Parameters
        ----------
        hist_dict
            Dictionary mapping histogram names to functions. This class
            does not further assumptions about the histograms other than
            the functions are used to fill the histogram. The functions should
            return the histogram that was filled.
        is_mc
            Whether the data to be used is simulation
        dsname
            Data set name of the data to be used
        dsname_in_hist
            Data set name to be used for the data when filling the histogram.
            Usually this is the same as ``dsname``
        sys_enabled
            Whether systematic uncertainties are to be computed if available.
            For example, this might be False if the user disabled systematic
            computation in their config
        sys_overwrite
            If not None, the data that is being filled is treated as a
            systematic variation. The parameter then names the systematic
            variation. An example would be JER variations or events generated
            with different generator settings
        cuts_to_histogram
            Lists cuts for which histograms should be produced. If ``None``
            all cuts will create histograms
        systs_to_histogram
            Lists systematics which should be included in histograms.
            If ``None``, all systematics will be included.
        """
        self.output = {
            "hists": {},
            "cutflows": defaultdict(AddableDict)
        }
        if hist_dict is None:
            self.hist_dict = {}
        else:
            self.hist_dict = hist_dict
        self.is_mc = is_mc
        self.dsname = dsname
        self.dsname_in_hist = dsname_in_hist
        self.sys_enabled = sys_enabled
        self.sys_overwrite = sys_overwrite
        self.cuts_to_histogram = cuts_to_histogram
        self.systs_to_histogram = systs_to_histogram
        self.done_hists = set()

    def fill_cutflows(self, data, systematics, cut, done_steps, cats):
        """Fill the cutflows for a specific step or cut

        Parameters
        ----------
        data
            Event data, usually NanoEvents
        systematics
            Record array with systematic weights
        cut
            Name of the cut applied last
        done_steps
            List of steps that have been done. Each step is represented by a
            string
        cats
            Categorizations to split the numbers into. The keys of cats name
            the categorization, while the values are lists, of which each
            element names a category. If not all categories are present in
            ``data``, the particular categorization is ignored.
        """
        if self.sys_overwrite is not None:
            return
        accumulator = self.output["cutflows"]
        if cut in accumulator[self.dsname]:
            return
        if systematics is not None:
            weight = systematics["weight"]
        else:
            weight = ak.Array(np.ones(len(data)))
            if hasattr(data.layout, "bytemask"):
                weight = weight.mask[~ak.is_none(data)]
        data_fields = ak.fields(data)
        # Skip cats that are missing in data
        cats = {k: v for k, v in cats.items()
                if all(field in data_fields for field in v)}
        axes = [hi.axis.StrCategory([], name=cat, label=cat, growth=True)
                for cat in cats.keys()]
        hist = hi.Hist(hi.axis.Integer(0, 1), *axes, storage="Weight")
        if len(cats) > 0:
            for cat_position in itertools.product(*list(cats.values())):
                masks = []
                for pos in cat_position:
                    masks.append(np.asarray(ak.fill_none(data[pos], False)))
                mask = np.bitwise_and.reduce(masks, axis=0)
                count = ak.sum(weight[mask])
                args = {name: pos
                        for name, pos in zip(cats.keys(), cat_position)}
                hist.fill(0, **args, weight=count)
        else:
            hist.fill(0, weight=ak.sum(weight))
        count = hist.values().sum()
        if logger.getEffectiveLevel() <= logging.INFO:
            num_rows = len(data)
            num_masked = ak.sum(ak.is_none(data))
            logger.info(f"Filling cutflow. Current event count: {count} "
                        f"({num_rows} rows, {num_masked} masked)")
        accumulator[self.dsname][cut] = hist

    def _add_hist(self, cut, histname, sysname, dsname, hist):
        """Add a histogram to ``self.output``

        Parameters
        ----------
        cut
            Name of the cut applied last
        histname
            Name of the histogram
        sysname
            Name of the systematic variation if this histogram is for
            a systematic variation only
        dsname
            Name of the data set
        hist
            The histogram to add
        """
        acc = self.output["hists"]
        # Split histograms by data set name. Summing histograms of the same
        # data set is generally much faster than summing across data sets
        # because normally the category axes for one data set are always
        # the same. Thus not summing across data sets will increase speed
        # significantly.
        if dsname not in acc:
            acc[dsname] = {}
        if (cut, histname) in acc[dsname]:
            try:
                acc[dsname][(cut, histname)] += hist
            except ValueError as err:
                raise ValueError(
                    f"Error adding sys {sysname} to hist {histname} for cut"
                    f" {cut} due to incompatible axes. This most likely"
                    f" caused by setting a new column after setting a new "
                    f" channel and a new systematic, consider setting "
                    f"'no_callback=True' on all new columns set before the"
                    f" next cut") from err
        else:
            acc[dsname][(cut, histname)] = hist
        self.done_hists.add((cut, histname, sysname))

    def _is_sysname_allowed(self, sysname):
        """
        Helper method to check whether a systematic should be included in
        a histogram.

        Parameters
        ----------
        sysname
            The name of the systematic (e.g. "muonsf_down").
        """
        if self.systs_to_histogram is None:
            return True
        else:
            if "_" in sysname:
                sysname_split = sysname.split("_")
                if sysname_split[-1] == "down" or sysname_split[-1] == "up" \
                        or sysname_split[-1].isdigit():
                    sysname = "_".join(sysname_split[:-1])
            return sysname in self.systs_to_histogram

    def fill_hists(self, data, systematics, cut, done_steps, cats):
        """Fill the histograms for a specific step or cut

        This will call the functions found in ``self.hist_dict``.

        Parameters
        ----------
        data
            Event data, usually NanoEvents
        systematics
            Record array with systematic weights
        cut
            Name of the cut applied last
        done_steps
            List of steps that have been done. Each step is represented by a
            string
        cats
            Categorizations to split the numbers into. The keys of cats name
            the categorization, while the values are lists, of which each
            element names a category.
        """
        if self.cuts_to_histogram is not None:
            if cut not in self.cuts_to_histogram:
                return
        do_systematics = self.sys_enabled and systematics is not None
        if do_systematics and self.sys_overwrite is None:
            weight = {}
            for syscol in ak.fields(systematics):
                if syscol == "weight":
                    sysname = "nominal"
                    sysweight = systematics["weight"]
                elif not self._is_sysname_allowed(syscol):
                    continue
                else:
                    sysname = syscol
                    sysweight = systematics["weight"] * systematics[syscol]
                weight[sysname] = sysweight
        elif systematics is not None:
            weight = systematics["weight"]
        elif self.sys_enabled:
            weight = {"nominal": None}
        else:
            weight = None
        for histname, fill_func in self.hist_dict.items():
            if (fill_func.step_requirement is not None
                    and fill_func.step_requirement not in done_steps):
                continue
            try:
                if self.sys_overwrite is not None:
                    sysname = self.sys_overwrite
                    # But only if we want to compute systematics
                    if do_systematics and fill_func.do_systs:
                        if (cut, histname, sysname) in self.done_hists:
                            continue
                        elif not self._is_sysname_allowed(sysname):
                            continue

                        sys_hist = fill_func(
                            data=data, categorizations=cats,
                            dsname=self.dsname_in_hist, is_mc=self.is_mc,
                            weight={sysname: weight})
                        self._add_hist(cut, histname, sysname, self.dsname,
                                       sys_hist)
                else:
                    if (cut, histname, None) in self.done_hists:
                        continue
                    hist = fill_func(
                        data=data, categorizations=cats,
                        dsname=self.dsname_in_hist, is_mc=self.is_mc,
                        weight=weight)
                    self._add_hist(cut, histname, None, self.dsname, hist)
            except pepper.hist_defns.HistFillError:
                # Ignore if fill is missing in data
                continue

    def get_callbacks(self):
        """Get all functions that should be called after every step in a
        selection is done"""
        return [self.fill_cutflows, self.fill_hists]

    @property
    def channels(self):
        """Deprecated"""
        raise AttributeError("'channels' is not used anymore. Use "
                             "Selector.set_cat('channel', [...])")
