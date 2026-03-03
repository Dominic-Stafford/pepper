#!/usr/bin/env python3

import abc
from collections.abc import Mapping
import numpy as np
import coffea
from coffea.lookup_tools.extractor import file_converters
from coffea.btag_tools import BTagScaleFactor
import awkward as ak
import correctionlib
from functools import partial
from itertools import islice
import warnings
import logging

from pepper.misc import onedimeval, LHCRun, get_run_for_year

logger = logging.getLogger(__name__)


def get_evaluator(filename, fileform=None, filetype=None):
    """Get the Coffea evaluator for a specific type of scale factor

    Parameters
    ----------
    filename
        File name containing the scale factor. The files will be opened and
        read by a Coffea converter
    fileform
        The format of the file. If None, determined by the file extension in
        ``filename``
    filetype
        Specific type of scale factor contained in the file. For example "junc"
        for jet energy uncertainties

    Returns
    -------
        Coffea evaluator loaded with the scale factor
    """
    if fileform is None:
        fileform = filename.split(".")[-1]
    if filetype is None:
        filetype = "default"
    converter = file_converters[fileform][filetype]
    extractor = coffea.lookup_tools.extractor()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        for key, value in converter(filename).items():
            extractor.add_weight_set(key[0], key[1], value)
    extractor.finalize()
    return extractor.make_evaluator()


class ScaleFactors:
    """Scale factor from n-dimensional arrays, for example from histograms"""
    def __init__(self, factors, factors_up, factors_down, bins):
        """
        Paramters
        ---------
        factors
            Numpy array with the scale factors for the central variation
        factors_up
            Numpy array with the scale factors for the up variation
        factors_down
            Numpy array with the scale factors for the down variation
        bins
            A dict whose length is equal to the number of dimensions of
            factors. Its nth key gives a name to the nth dimension, so
            that it can be used as a kwarg in the __call__ method of this
            class. "variation" can not be a key. Its values are Numpy
            arrays and determine the edges of the binnding used for the
            factors.
        """
        if "variation" in bins:
            raise ValueError("'variation' must not be in bins")
        for i, (axis, edges) in enumerate(bins.items()):
            if factors.shape[i] + 1 != len(edges):
                raise ValueError(
                    f"axis '{axis}' should have {factors.shape[i] + 1} "
                    f"edges but has {len(edges)}")
        if (factors_up is not None and factors_down is not None and
            (factors.shape != factors_up.shape
                or factors.shape != factors_down.shape)):
            raise ValueError(
                "factors_up or factors_down have inconsistent shape")
        self._factors = factors
        self._factors_up = factors_up
        self._factors_down = factors_down
        self._bins = bins

    @staticmethod
    def _setoverflow(factors, value):
        """Set the overflow (value to use when events are outside the range)
        for a given ``factors`` array to a given ``value``"""
        for i in range(factors.ndim):
            factors[tuple([slice(None)] * i + [slice(0, 1)])] = value
            factors[tuple([slice(None)] * i + [slice(-1, None)])] = value

    @classmethod
    def from_hist(cls, hist, dimlabels=None):
        """Create a new instance from a histogram

        This uses the variances found in the histogram as up/down variations.

        Parameters
        ----------
        hist
            The uproot TH1 the scale factors are extracted from
        dimlabels
            The names of the variables the scale factor depends on. If
            ``None``, use the name of the axes of the histogram
        """
        edges = hist.to_numpy(flow=True)[1:]
        if dimlabels is None:
            dimlabels = []
            for member in ("fXaxis", "fYaxis", "fZaxis")[:len(edges)]:
                dimlabels.append(hist.all_members[member].all_members["fName"])
        factors = hist.values(flow=True)
        sigmas = np.sqrt(hist.variances(flow=True))
        factors_up = factors + sigmas
        factors_down = factors - sigmas
        if len(edges) != len(dimlabels):
            raise ValueError(
                f"Got {len(edges)} dimensions but {len(dimlabels)} labels")
        # Set overflow bins to 1 so that events outside the histogram
        # get scaled by 1
        cls._setoverflow(factors, 1)
        if factors_up is not None:
            cls._setoverflow(factors_up, 1)
        if factors_down is not None:
            cls._setoverflow(factors_down, 1)
        bins = dict(zip(dimlabels, edges))
        return cls(factors, factors_up, factors_down, bins)

    @classmethod
    def from_hist_per_bin(cls, hists, dimlabels=None, **kwargs):
        def get_factors(depth, hists):
            if depth == 0:
                factors = hists.values(flow=True)
                cls._setoverflow(factors, 1)
                return factors
            else:
                return [get_factors(depth-1, h) for h in hists]

        def get_sigmas(depth, hists):
            if depth == 0:
                sigmas = np.sqrt(hists.variances(flow=True))
                cls._setoverflow(sigmas, 0)
                return sigmas
            else:
                return [get_factors(depth-1, h) for h in hists]

        tmp_hists = hists
        for _i, _k in enumerate(kwargs):
            tmp_hists = tmp_hists[0]
        hist_edges = list(tmp_hists.to_numpy(flow=True)[1:])
        edges = list(kwargs.values()) + hist_edges

        if dimlabels is None:
            dimlabels = list(kwargs.keys())
            for member in ("fXaxis", "fYaxis", "fZaxis")[:len(hist_edges)]:
                dimlabels.append(
                    tmp_hists.all_members[member].all_members["fName"])
        else:
            dimlabels = list(kwargs.keys()) + dimlabels

        factors = ak.from_iter(get_factors(len(kwargs), hists))
        sigmas = ak.from_iter(get_sigmas(len(kwargs), hists))
        factors_up = factors + sigmas
        factors_down = factors - sigmas
        if len(edges) != len(dimlabels):
            raise ValueError(
                f"Got {len(edges)} dimenions but {len(dimlabels)} labels")
        bins = dict(zip(dimlabels, edges))
        return cls(factors, factors_up, factors_down, bins)

    def __call__(self, variation="central", **kwargs):
        """Evaluate the scale factor

        Parameters
        ----------
        variation
            Direction of the systematic variation. One of "central", "up"
            or "down"
        **kwargs
            The parameters the scale factor depends on. For example
            the pt, eta and so on.

        Returns
        -------
            Array of scale factors
        """
        if variation not in ("central", "up", "down"):
            raise ValueError("variation must be one of 'central', 'up', "
                             "'down'")
        factors = self._factors
        if variation == "up":
            factors = self._factors_up
        elif variation == "down":
            factors = self._factors_down

        bin_idxs = []
        counts = None
        for key, bins_for_key in self._bins.items():
            try:
                val = kwargs[key]
            except KeyError:
                raise ValueError(
                    f"Scale factor depends on \"{key}\" but no such "
                    f"argument was not provided")
            if isinstance(val, ak.Array):
                counts = []
                for _i in range(val.ndim - 1):
                    counts.append(ak.num(val))
                    val = ak.flatten(val)
                val = np.asarray(val)
            bin_idxs.append(np.digitize(val, bins_for_key) - 1)
        ret = factors[tuple(bin_idxs)]
        if counts is not None:
            for count in reversed(counts):
                ret = ak.unflatten(ret, count)
        return ret

    @property
    def dimlabels(self):
        """The names of the variables the scale factor depends on"""
        return self._bins.keys()


class CorrLibSFs:
    def __init__(self, sysnaming, jsonfile, corrname, add_args):
        """Create a thin wrapper around a single correctionlib correction.

        This object extracts the named correction from a correctionlib
        CorrectionSet (either provided directly or loaded from a JSON file)
        and records the expected input names so the correction can be
        evaluated later using the `__call__` method.

        Parameters
        ----------
        sysnaming
            Mapping from variation name (e.g. "central", "up", "down") to
            the value that should be supplied for correction inputs named
            "ValType", "scale_factors", "weights" or "systematic". The
            mapped value is inserted automatically based on the requested
            variation when evaluating the correction.
        jsonfile
            Either a path to a correctionlib JSON file or a
            ``correctionlib.CorrectionSet`` instance containing the
            corrections. If a filename is given the correction set will be
            loaded from that file.
        corrname
            The name of the correction within the correction set to use.
        add_args
            Optional dict of additional fixed arguments to pass to the
            correction. If the value for a given input name is a Mapping,
            it is interpreted as a set of "correction ranges" mapping a
            correction key to a (low, high) pair; these ranges are used by
            ``__call__`` to evaluate different correction keys over
            disjoint ranges (e.g. different pt bins). Otherwise the value
            is passed directly as the input for that parameter.
        """
        self.corrname = corrname
        if isinstance(jsonfile, correctionlib.CorrectionSet):
            try:
                self.corrset = jsonfile[corrname]
            except IndexError as e:
                raise ValueError(
                    f"Correction '{corrname}' not found in file "
                    f"'{jsonfile}'") from e
        elif isinstance(jsonfile, str):
            self.corrset = self._read_correctionfile(jsonfile, corrname)
        else:
            raise ValueError(f"Invalid correctionlib file: {jsonfile}")
        self.names = [i.name for i in self.corrset.inputs]
        self.add_args = add_args
        self.sysnaming = sysnaming

    @staticmethod
    def _read_correctionfile(jsonfile, corrname):
        try:
            corrset = correctionlib.CorrectionSet.from_file(jsonfile)[corrname]
        except IndexError as e:
            raise ValueError(
                f"Correction '{corrname}' not found in file "
                f"'{jsonfile}'") from e
        return corrset

    def __repr__(self) -> str:
        return f"CorrLibSFs(corrname={self.corrname})"

    def __call__(self, variation="central", **kwargs):
        """Evaluate the scale factor

        Parameters
        ----------
        variation
            Direction of the systematic variation. One of "central", "up"
            or "down"
        **kwargs
            The parameters the scale factor depends on. For example
            the pt, eta and so on.

        Returns
        -------
            Array of scale factors
        """
        if variation not in ("central", "up", "down"):
            raise ValueError("variation must be one of 'central', 'up', "
                             "'down'")
        arrays = dict()  # Since python 3.7 dicts are ordered
        output_like = None
        correction_ranges = []
        correction_ranges_arg_idx = -1
        for i, in_name in enumerate(self.names):
            if in_name in ["ValType", "scale_factors", "weights", "systematic"]:
                arrays[in_name] = self.sysnaming[variation]
            elif in_name in self.add_args:
                # Check if arg is a mapping
                if isinstance(self.add_args[in_name], Mapping):
                    # Save the argument index of the correction ranges. ie
                    # what index to put the argument when calling evaluate.
                    correction_ranges_arg_idx = i
                    # Insert all correction ranges into the arrays
                    for corr, (low, high) in self.add_args[in_name].items():
                        correction_ranges.append(
                                (corr, float(low), float(high)))
                else:
                    arrays[in_name] = self.add_args[in_name]
            elif in_name in kwargs:
                arrays[in_name] = kwargs[in_name]
                output_like = i
            else:
                raise ValueError(
                    f"Scale factor depends on \"{in_name}\" but no such "
                    f"argument was not provided")
        if correction_ranges_arg_idx == -1:
            # No correction ranges were provided, so we can just evaluate
            # the scale factor
            return onedimeval(self.corrset.evaluate, *arrays.values(),
                              output_like=output_like)
        counts = ak.num(arrays["pt"])
        flattened_arrays = dict()  # Since python 3.7 dicts are ordered
        for key, arr in arrays.items():
            if isinstance(arr, str):
                flattened_arrays[key] = arr
            else:
                flattened = arr
                for _i in range(flattened.ndim - 1):
                    flattened = ak.flatten(flattened)
                flattened_arrays[key] = np.asarray(flattened)
        pt = flattened_arrays["pt"]
        result = np.empty_like(pt, dtype=np.float64)
        for key, low, high in correction_ranges:
            mask = (low <= pt) & (pt < high)
            if not np.any(mask):
                continue
            eval_args = list(flattened_arrays.values())
            # Only use masked data for all input arrays i.e. mask non-strings
            for i, arr in enumerate(flattened_arrays.values()):
                if isinstance(arr, np.ndarray):
                    eval_args[i] = arr[mask]
            # Insert the correction range key at the correct position
            eval_args.insert(correction_ranges_arg_idx, key)
            chunk_result = self.corrset.evaluate(*eval_args)
            result[mask] = chunk_result
        return ak.unflatten(result, counts)

    @property
    def dimlabels(self):
        """The names of the variables the scale factor depends on"""
        return [in_name for in_name in self.names if (
                (in_name not in ["ValType", "scale_factors"])
                and in_name not in self.add_args)]


class MuonScaleFactor:
    """Scale factors for muons that come with split statistical and systematic
    uncertainty sources"""
    def __init__(self, nominal, stat, syst):
        """
        Parameters
        ----------
        nominal
            ScaleFactors for nominal weights
        stat
            ScaleFactors for deriving weights with statistical variation
        syst
            ScaleFactors for deriving weights with systematic variation
        """
        self.nominal = nominal
        self.stat = stat
        self.syst = syst

    def __call__(self, variation="central", **kwargs):
        """Evaluate the scale factor

        Parameters
        ----------
        variation
            Direction of the systematic variation. One of "central", "up",
            "down", "syst up", "syst down", "stat up" or "stat_down"
        **kwargs
            The parameters the scale factor depends on. For example
            the pt

        Returns
        -------
            Array of scale factors
        """
        if variation not in ("central", "up", "down", "syst up", "syst down",
                             "stat up", "stat down"):
            raise ValueError(f"Invalid variation '{variation}'")
        if variation in ("central", "up", "down"):
            sf = self.nominal
        else:
            unctype, variation = variation.split(" ")
            if unctype == "syst":
                sf = self.syst
            else:
                sf = self.stat
        return sf(variation=variation, **kwargs)

    @property
    def dimlabels(self):
        """The names of the variables the scale factor depends on"""
        return self.nominal.dimlabels


class BTagWPs:
    """A helper class to access b-tagging working points

    Attributes
    ----------
    name
        The name of the tagger (with appropriate capitilsation)
    """
    def __init__(self, wp_dict, tagger, year):
        """
        Parameters
        ----------
        wp_dict
            A dicitonary defining for each tagger the name (with appropriate
            capitilsation), the available working points, and the values
        tagger
            Name of the tagger, for example "deepCSV" or "deepJet".
            Case insensitive
        year
            Year of the data
        """
        self.wps = wp_dict[tagger]["wps"][year]
        # Slice the wpnames to the number of wps available for this year
        self.wpnames = dict(islice(
            wp_dict[tagger]["wpnames"].items(), len(self.wps)))
        self.name = wp_dict[tagger]["name"]

    def __call__(self, wp):
        wp = wp.lower()
        if wp not in self.wpnames.keys():
            raise ValueError(f"IInvalid value for b-tag WP. Expected one of "
                             f"{list(self.wpnames.keys())}")
        wp_idx = list(self.wpnames.keys()).index(wp)
        wp_val = self.wps[wp_idx]
        wp_letter = self.wpnames[wp].upper()
        return wp_letter, wp_idx, wp_val


class BTagWeighter:
    """Calculate the weights for b tagging

    This implements the methods described in
    https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods

    Attributes
    ----------
    sources
        List of uncertainty sources available
    """
    def __init__(self, sf_filename, eff_filename, tagger, year, wps,
                 method="fixedwp", meastype="mujets", ignore_missing=False):
        """
        Parameters
        ----------
        sf_filename
            File name of the scale factors made by the POG
        eff_filename
            Path to file with efficiencies of the analysis b tagging cut.
            In Pepper the file can be generated by the
            generate_btag_efficiencies.py script
        tagger
            Name of the tagger, for example "deepCSV" or "deepJet".
            Case insensitive
        year
            Year of the data
        wps
            A BTagWPs object defining the wps for this tagger and year, and
            the name of the tagger (with appropriate capitilsation)
        method
            Method to use, eithwe "fixedwp" or "iterativefit". "fixedwp" is
            used when the analysis is using one of the POG defined working
            points. "iterativefit" should be used when the b-tagging
            descriminant itself is used in the analysis.
        meastype
            Which type of scale factors provided by the POG should be used.
        ignore_missing
            If False, raise error if weight could not be computed, else
            return 1.
        """

        if isinstance(method, str):
            method = method.lower()
            possible_methods = ("iterativefit", "fixedwp")
            if method not in possible_methods:
                raise ValueError(
                    "Method must be one of: " + ", ".join(possible_methods)
                )
        else:
            raise TypeError(f"Expected str for method, got {method}")

        tagger = tagger.lower()

        if sf_filename.endswith(".csv"):
            self.filetype = "csv"
            self.sf = []
            if method == "fixedwp":
                self.eff_evaluator = get_evaluator(eff_filename)
                for i in range(3):
                    # Note that for light flavor normally only inclusive is
                    # available, thus we fix it to 'incl' here
                    self.sf.append(BTagScaleFactor(
                        sf_filename, i, f"{meastype},{meastype},incl"))
            else:
                self.sf.append(BTagScaleFactor(
                    sf_filename, BTagScaleFactor.RESHAPE,
                    "iterativefit,iterativefit,iterativefit"))
        elif sf_filename.endswith(".json") or sf_filename.endswith(".json.gz"):
            self.filetype = "json"
            if method == "fixedwp":
                self.eff_evaluator = get_evaluator(eff_filename)
                correset = correctionlib.CorrectionSet.from_file(sf_filename)
                if get_run_for_year(year) == LHCRun.Run3:
                    self.sf = [correset[f"{wps.name}_{meastype}"],
                               correset[f"{wps.name}_light"]]
                else:
                    self.sf = [correset[f"{wps.name}_{meastype}"],
                               correset[f"{wps.name}_incl"]]
            else:
                raise NotImplementedError(
                    "iterativefit with json format is not implemented")

        self.method = method
        if meastype == "mujets":
            self.sources = ["jes", "pileup", "statistic", "type3"]
        elif meastype == "comb":
            self.sources = ["isr", "fsr", "hdamp", "jes", "jer", "pileup",
                            "qcdscale", "statistic", "topmass", "type3"]
        else:
            self.sources = ["correlated", "uncorrelated"]
        self.wps = wps
        self.ignore_missing = ignore_missing

    def _fixedwp(self, wp, jf, eta, pt, discr, variation, efficiency):
        if isinstance(wp, str):
            _, wp, wp_val = self.wps(wp)
        elif not isinstance(wp, int):
            raise TypeError(f"Expected int or str for wp, got {wp}")
        elif wp > len(self.wps):
            raise ValueError(f"Expected value between 0 and {len(self.wps)} "
                             f"for wp, got {wp}")
        possible_variations = (
            "central", "light up", "light down", "heavy up", "heavy down",
            "light up_correlated", "light up_uncorrelated",
            "light down_correlated", "light down_uncorrelated",
            "heavy up_correlated", "heavy up_uncorrelated",
            "heavy down_correlated", "heavy down_uncorrelated")
        if variation not in possible_variations:
            raise ValueError(
                "variation must be one of: " + ", ".join(possible_variations))
        light_vari = "central"
        heavy_vari = "central"
        if variation != "central":
            vari_type, direction = variation.split(" ")
            if vari_type == "light":
                light_vari = direction
            else:
                heavy_vari = direction
        sf = self.sf[wp].eval(
            heavy_vari, jf, eta, pt, discr, self.ignore_missing)
        if light_vari != heavy_vari:
            sf = ak.where(
                jf >= 4, sf,
                self.sf[wp].eval(
                    light_vari, jf, eta, pt, discr, self.ignore_missing))

        eff = self.eff_evaluator[efficiency](jf, pt, abs(eta))
        sfeff = sf * eff
        is_tagged = discr > wp_val

        p_mc = ak.prod(eff[is_tagged], axis=1) * ak.prod(
            (1 - eff)[~is_tagged], axis=1)
        p_data = ak.prod(sfeff[is_tagged], axis=1) * ak.prod(
            (1 - sfeff)[~is_tagged], axis=1)

        # TODO: What if one runs into numerical problems here?
        return p_data / p_mc

    def _fixedwp_json(self, wp, jf, eta, pt, discr, variation, efficiency):
        def evaluate_two_flavor(jf, abseta, pt):
            is_light = jf == 0
            isnt_light = ~is_light
            sf = np.empty(len(is_light))
            sf[isnt_light] = self.sf[0].evaluate(
                heavy_vari, wp, jf[isnt_light], abseta[isnt_light],
                pt[isnt_light])
            sf[is_light] = self.sf[1].evaluate(
                light_vari, wp, jf[is_light], abseta[is_light], pt[is_light])
            return sf

        wp, _, wp_val = self.wps(wp)
        light_vari = "central"
        heavy_vari = "central"
        if variation != "central":
            vari_type, direction = variation.split(" ")
            if vari_type == "light":
                light_vari = direction
            elif vari_type == "heavy":
                heavy_vari = direction
            else:
                raise ValueError(f"Unknown variation type: {vari_type}")
        abseta = abs(eta)
        sf = onedimeval(evaluate_two_flavor, jf, abseta, pt)

        eff = self.eff_evaluator[efficiency](jf, pt, abseta)
        sfeff = sf * eff
        is_tagged = discr > wp_val

        p_mc = ak.prod(eff[is_tagged], axis=1) * ak.prod(
            (1 - eff)[~is_tagged], axis=1)
        p_data = ak.prod(sfeff[is_tagged], axis=1) * ak.prod(
            (1 - sfeff)[~is_tagged], axis=1)

        return p_data / p_mc

    def _iterativefit(self, wp, jf, eta, pt, discr, variation):
        possible_variations = (
            "central",
            "up_lf",
            "down_lf",
            "up_hf",
            "down_hf",
            "up_hfstats1",
            "down_hfstats1",
            "up_lfstats1",
            "down_lfstats1",
            "up_hfstats2",
            "down_hfstats2",
            "up_lfstats2",
            "down_lfstats2",
            "up_cferr1",
            "up_cferr2",
            "down_cferr1",
            "down_cferr2",
            "up_jes",
            "down_jes",
            "up_AbsoluteMPFBias",
            "down_AbsoluteMPFBias",
            "up_AbsoluteScale",
            "down_AbsoluteScale",
            "up_AbsoluteStat",
            "down_AbsoluteStat",
            "up_RelativeBal",
            "down_RelativeBal",
            "up_RelativeFSR",
            "down_RelativeFSR",
            "up_RelativeJEREC1",
            "down_RelativeJEREC1",
            "up_RelativeJEREC2",
            "down_RelativeJEREC2",
            "up_RelativeJERHF",
            "down_RelativeJERHF",
            "up_RelativePtBB",
            "down_RelativePtBB",
            "up_RelativePtEC1",
            "down_RelativePtEC1",
            "up_RelativePtEC2",
            "down_RelativePtEC2",
            "up_RelativePtHF",
            "down_RelativePtHF",
            "up_RelativeStatEC",
            "down_RelativeStatEC",
            "up_RelativeStatFSR",
            "down_RelativeStatFSR",
            "up_RelativeStatHF",
            "down_RelativeStatHF",
            "up_PileUpDataMC",
            "down_PileUpDataMC",
            "up_PileUpPtBB",
            "down_PileUpPtBB",
            "up_PileUpPtEC1",
            "down_PileUpPtEC1",
            "up_PileUpPtEC2",
            "down_PileUpPtEC2",
            "up_PileUpPtHF",
            "down_PileUpPtHF",
            "up_PileUpPtRef",
            "down_PileUpPtRef",
            "up_FlavorQCD",
            "down_FlavorQCD",
            "up_Fragmentation",
            "down_Fragmentation",
            "up_SinglePionECAL",
            "down_SinglePionECAL",
            "up_SinglePionHCAL",
            "down_SinglePionHCAL",
            "up_TimePtEta",
            "down_TimePtEta"
        )
        if variation not in possible_variations:
            raise ValueError(
                "variation must be one of: "
                + ", ".join(possible_variations))
        sf = self.sf[0].eval(
            variation, jf, eta, pt, discr, self.ignore_missing)
        # Apply a SF of 1. to c-flavored jets
        is_c = abs(jf) == 4
        if variation == "central":
            sf = ak.fill_none(ak.mask(sf, ~is_c), 1.0)
        # Apply cferr uncertainties only to c-flavored jets,
        # otherwise use central variation
        if "cferr" in variation:
            sf_central = self.sf[0].eval(
                "central", jf, eta, pt, discr, self.ignore_missing)
            sf = ak.where(is_c, sf, sf_central)
        # Apply lf, hf and jes uncertainties only to b- and
        # usdg-flavored jets, otherwise use central variation
        if "lf" in variation or "hf" in variation:
            sf_central = self.sf[0].eval(
                "central", jf, eta, pt, discr, self.ignore_missing)
            sf = ak.where(is_c, sf_central, sf)
        sf = ak.prod(sf, axis=1)
        return sf

    def __call__(
            self, wp, jf, eta, pt, discr, variation="central",
            efficiency="central"):
        """Evaluate the weight for the b tagging

        Parameters
        ----------
        wp
            Working point of the cut. Ignored if ``self.method``
            is ``iterativefit``
        jf
            Jet flavor, same numbering as the PDG ID, but without sign
        eta
            Eta of the jet four momentum
        pt
            pt of the jey four momentum
        discri
            b-tagging discriminator value of the jets
        variation
            Uncertainty variation to do. "central" do derive nominal weights
        efficiency
            Efficiency can depend on the systematic variation that is done.
            This specifies the name of the varition to use the efficiency from.
            The provided efficieny file must have a histogram named as such.

        Returns
        -------
            Array of b-tagging weights
        """
        if self.method == "fixedwp" and self.filetype == "csv":
            sf = self._fixedwp(
                wp, jf, eta, pt, discr, variation, efficiency)
        elif self.method == "fixedwp" and self.filetype == "json":
            sf = self._fixedwp_json(
                wp, jf, eta, pt, discr, variation, efficiency)
        elif self.method == "iterativefit":
            sf = self._iterativefit(wp, jf, eta, pt, discr, variation)
        else:
            raise ValueError("Unknown method or filetype")
        return sf

    @property
    def available_efficiencies(self):
        """Different available efficiencies as found in the provided efficiency
        file. Efficiencies can depend on the systematic variation"""
        return set(self.eff_evaluator.keys())


class JetPuIdWeighter:
    """Compute weights for the jet pileup ID"""
    def __init__(self, sf_filename, eff_filename=None):
        """
        Parameters
        ----------
        sf_filename
            Name of the scale factors probided by the POG in form of a JSON
        eff_filename
            File with efficiencies
        """
        self.sf_evaluator = correctionlib.CorrectionSet.from_file(sf_filename)
        if "PUJetID_mis" in [k for k in self.sf_evaluator.keys()]:
            self.has_mis_prob = True
        else:
            # Mistagging probabilties are not currently provided in UL,
            # recommendation from the POG is to just not consider mistags
            self.has_mis_prob = False
        if eff_filename is not None:
            self.eff_evaluator = get_evaluator(eff_filename)
        else:
            self.eff_evaluator = None

    def __call__(self, wp, eta, pt, pass_puid, has_gen_jet,
                 sf_type="eff", variation="nom"):
        """Compute the jet pileup ID weights for the specified data

        Parameters
        ----------
        wp
            The working point of the used ID cut
        eta
            Eta of the jet four vectors
        pt
            pt of the jet four vectors
        pass_puid
            Whether the jets pass the pileup ID cut
        has_gen_jet
            Whether a jet as a generator-level jet associated to it
        sf_type
            Type of scale factors to use. Either "eff" or "mis"
        variation
            Systematic variation to do. "nom" for nominal weights
        """
        if sf_type == "eff":
            eta = eta[has_gen_jet]
            pt = pt[has_gen_jet]
            pass_puid = pass_puid[has_gen_jet]
        elif sf_type == "mis":
            if not self.has_mis_prob:
                raise ValueError(
                    "Cannot compute SFs for mistagged PU jets as these are "
                    "not provided in this PU ID SF file")
            eta = eta[~has_gen_jet]
            pt = pt[~has_gen_jet]
            pass_puid = pass_puid[~has_gen_jet]
        wp = wp.lower()
        if wp == "loose":
            wp = "L"
        elif wp == "medium":
            wp = "M"
        elif wp == "tight":
            wp = "T"
        else:
            raise ValueError("Invalid value for wp. Expected 'loose', "
                             "'medium' or 'tight'")

        def sf_evaluate(var, eta, pt):
            return self.sf_evaluator["PUJetID_" + sf_type].evaluate(
                eta, pt, var, wp)

        sf = onedimeval(partial(sf_evaluate, variation), eta, pt)
        # Clip SF if greater than 5, as recommended by JetMET POG (particularly
        # relevant for mistagging SFs)
        sf = ak.where(sf > 5, 5, sf)
        if self.eff_evaluator is not None:
            eff = self.eff_evaluator[sf_type](pt, eta)
        else:
            try:
                eff = onedimeval(partial(sf_evaluate, "MCEff"), eta, pt)
            except IndexError:
                raise KeyError("No MC efficiencies in this Jet PU ID SF "
                               "file, please provide separate ones")

        sfeff = sf * eff
        p_mc = ak.prod(eff[pass_puid], axis=1) * ak.prod(
            (1 - eff)[~pass_puid], axis=1)
        p_data = ak.prod(sfeff[pass_puid], axis=1) * ak.prod(
            (1 - sfeff)[~pass_puid], axis=1)
        return p_data/p_mc


class PileupWeighter:
    """Compute weights for the pileup reweighting"""
    def __init__(self, rootfile):
        """
        Paramteres
        ----------
        rootfile
            Path to a Root file containing the weights. In pepper this file
            can be generated by the generate_pileup_weights.py script
        """
        self.central = {}
        self.up = {}
        self.down = {}

        self.upsuffix = "_up"
        self.downsuffix = "_down"
        for key, hist in rootfile.items():
            key = key.rsplit(";", 1)[0]
            sf = ScaleFactors.from_hist(hist, ["ntrueint"])
            if key.endswith(self.upsuffix):
                self.up[key[:-len(self.upsuffix)]] = sf
            elif key.endswith(self.downsuffix):
                self.down[key[:-len(self.downsuffix)]] = sf
            else:
                self.central[key] = sf
        if (self.central.keys() != self.up.keys()
                or self.central.keys() != self.down.keys()):
            raise ValueError(
                "Missing up/down or central weights for some datasets")

    def __call__(self, dsname, NumTrueInteractions, variation="central"):  # noqa: N803
        """Compute the weights for given events

        Parameters
        ----------
        dsname
            Name of the data set
        NumTrueInteractions
            Number of true pileup per event. Pileup_ntrueint in NanoAOD

        Returns
        -------
            Array of weights
        """
        # If all_datasets is present, use that instead of per-dataset weights
        if "all_datasets" in self.central:
            key = "all_datasets"
        else:
            key = dsname
        if variation == "up":
            return self.up[key](ntrueint=NumTrueInteractions)
        elif variation == "down":
            return self.down[key](ntrueint=NumTrueInteractions)
        elif variation == "central":
            return self.central[key](ntrueint=NumTrueInteractions)
        else:
            raise ValueError("variation must be either 'up', 'down' or "
                             f"'central', not {variation}")


class TopPtWeigter:
    """Top pt reweighting according to
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
    """
    def __init__(self, method, scale=1.0, sys_only=False, **kwargs):
        """
        Parameters
        ----------
        method
            The method to use. Either "datanlo" or "theory"
        scale
            An overall scale to multiple the weights by
        sys_only
            Whether this reweighting should be used only as a systematic
            variation and not to scale nominal event weights
        **kwargs
            Parameters in the reweighting formulas. See ``datanlo_sf()`` and
            ``theory_sf``
        """
        if method.lower() == "datanlo":
            self.sffunc = self.datanlo_sf
        elif method.lower() == "theory":
            self.sffunc = self.theory_sf
        else:
            raise ValueError(f"Invalid method: {method}")
        self.scale = scale
        self.sys_only = sys_only
        self.kwargs = kwargs

    def datanlo_sf(self, pt):
        """Data-NLO method

        Uses the formula ``exp(a + b * pt)``. ``a`` and ``b`` are obtained
        from the ``kwargs`` attribute.

        Parameters
        ----------
        pt
            pt of the top four momentum
        """
        return np.exp(self.kwargs["a"] + self.kwargs["b"] * pt)

    def theory_sf(self, pt):
        """Theory method

        Uses the formula ``a * exp(b * pt) + c * pt + d``.
        ``a``, ``b``, ``c`` and ``d`` are obtained from the ``kwargs``
        attribute.
        """
        arg = self.kwargs
        return arg["a"] * np.exp(arg["b"] * pt) + arg["c"] * pt + arg["d"]

    def __call__(self, toppt, antitoppt):
        """Compute the weights for the top pt reweighting

        Weights for top quark and top antiquark are multiplied under square
        root.

        Parameters
        ----------
        toppt
            pt of the four momentum of the top quarks
        antitoppt
            pt of the four momentum of the top antiquarks

        Returns
        -------
            Array of one weight per event
        """
        sf = self.sffunc(toppt)
        antisf = self.sffunc(antitoppt)
        return np.sqrt(sf * antisf) * self.scale


class JetCorrectionAdapterBase(abc.ABC):
    """Base adapter for jet corrections"""

    @abc.abstractmethod
    def __call__(self, **kwargs):
        pass


class LegacyJetAdapter(JetCorrectionAdapterBase):
    """Generic adapter for legacy coffea jetmet_tools objects"""
    def __init__(self, legacy_obj, method_name):
        """
        Args:
            legacy_obj: The coffea jetmet_tools object
            method_name: Name of method to call (e.g., 'getCorrection', 'getUncertainty')
        """
        self.legacy_obj = legacy_obj
        self.method = getattr(legacy_obj, method_name)

    @property
    def levels(self):
        return self.legacy_obj.levels

    def __call__(self, **kwargs):
        return self.method(**kwargs)


class CorrectionlibAdapter(JetCorrectionAdapterBase):
    """
    Generic adapter for correctionlib corrections.
    Passes kwargs directly to correctionlib evaluate().
    """

    def __init__(self, path, correction_name: str, sysnaming={
            "central": "nom", "up": "up", "down": "down"},
            is_uncertainty: bool = False):
        """
        Args:
            cset: CorrectionSet from correctionlib
            correction_name: Correction name (str) or empty string/None
        """
        self.path = path
        self._is_uncertainty = is_uncertainty
        # Only allow exactly one correction_name
        if not isinstance(correction_name, str):
            msg = f"'correction_name' must be a string, got {type(correction_name)}."
            raise ValueError(msg)
        self.correction_name = correction_name
        # Create CorrLibSFs wrapper
        self.corr = CorrLibSFs(sysnaming, path, self.correction_name, {})

    def __call__(self, **kwargs):
        """Pass kwargs directly to correctionlib evaluate()"""
        if self.corr is None:
            # Return trivial result for empty correction
            if self._is_uncertainty:
                return 0.0
            else:
                return 1.0
        # Broadcast to match all array types.
        # Separate arrays from non-arrays while preserving order
        # Slightly involved, but we want to maintain original ordering of kwargs
        array_keys = [k for k, v in kwargs.items() if isinstance(v, (ak.Array, np.ndarray))]
        array_values = [kwargs[k] for k in array_keys]

        # Broadcast only the arrays
        if array_values:
            broadcasted = dict(zip(array_keys, ak.broadcast_arrays(*array_values)))

            # Rebuild kwargs in original order
            kwargs = {k: broadcasted.get(k, v) for k, v in kwargs.items()}
        return self.corr(**kwargs)


class CorrectionlibCompoundAdapter(JetCorrectionAdapterBase):
    """
    Adapter for compound corrections (JEC only).
    Applies multiple corrections sequentially: L1 * L2 * L3 * ...
    """

    def __init__(self, path, correction_names: str | list[str], sysnaming={
            "central": "sf", "up": "sfup", "down": "sfdown"}):
        """
        Args:
            cset: CorrectionSet from correctionlib
            correction_names: Single name (str) or list of names to apply sequentially
        """
        self.path = path

        # Normalize to list
        if isinstance(correction_names, str):
            self.correction_names = [correction_names] if correction_names else []
        else:
            self.correction_names = correction_names if correction_names else []

        # Create CorrLibSFs wrappers for each correction
        self.corrections = {}
        for name in self.correction_names:
            if name:
                try:
                    corr = CorrLibSFs(sysnaming, path, name, {})
                    self.corrections[name] = corr
                except Exception as e:
                    raise ValueError(
                        f"Correction '{name}' not found in CorrectionSet at {path}. "
                        f"Error: {e}"
                    )

    def __call__(self, **kwargs):
        """
        Apply corrections sequentially, updating JetPt after each correction.

        Args:
            **kwargs: Must include JetPt, and typically JetEta, JetA, Rho

        Returns:
            Combined correction factor
        """
        if "JetPt" not in kwargs:
            msg = "'JetPt' must be provided in kwargs for compound correction."
            raise ValueError(msg)

        if not self.corrections:
            return 1.0

        correction_factor = 1.0
        current_kwargs = kwargs.copy()
        for corr in self.corrections.values():
            # Check if we need to broadcast anything
            current_kwargs = dict(zip(
                current_kwargs.keys(),
                ak.broadcast_arrays(*current_kwargs.values())
            ))
            factor = corr(**current_kwargs)
            correction_factor = correction_factor * factor

            # Update pt for next correction in chain
            current_kwargs['JetPt'] = current_kwargs['JetPt'] * factor
        return correction_factor


class JetIdProducer:
    def __init__(self, jettype, jsonfile):
        """Module to determine jetID variables (passTight, passTightLepVeto),
        packed in Jet_jetId as in nanoAODv12

        Parameters
        ----------
        json
            jetID json file provided by the JetMET POG
        jettype
            jet type used in the analysis ("AK4PUPPI" or "AK4CHS")
        """
        self.corrset = correctionlib.CorrectionSet.from_file(jsonfile)
        self.tight = self.corrset[f"{jettype}_Tight"]
        self.tightLepVeto = self.corrset[f"{jettype}_TightLeptonVeto"]

    def evaluate(self, jet):

        multiplicity = jet.chMultiplicity + jet.neMultiplicity

        passTight = onedimeval(
            self.tight.evaluate,
            jet.eta,
            jet.chHEF,
            jet.neHEF,
            jet.chEmEF,
            jet.neEmEF,
            jet.muEF,
            jet.chMultiplicity,
            jet.neMultiplicity,
            multiplicity,
        )

        passTightLepVeto = onedimeval(
            self.tightLepVeto.evaluate,
            jet.eta,
            jet.chHEF,
            jet.neHEF,
            jet.chEmEF,
            jet.neEmEF,
            jet.muEF,
            jet.chMultiplicity,
            jet.neMultiplicity,
            multiplicity,
        )

        passTight = ak.values_astype(passTight, bool)
        passTightLepVeto = ak.values_astype(passTightLepVeto, bool)

        return passTight, passTightLepVeto
