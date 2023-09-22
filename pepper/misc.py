import os
from glob import glob
import inspect
import gc
from itertools import product
from functools import wraps, partial
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import awkward as ak
import hist as hi


"""
Attributes
----------
XROOTDTIMEOUT
    Seconds to wait for XRootD request to be replied by the server before
    raising an error
"""


XROOTDTIMEOUT = 10  # 10 s, no need to bother with slow sites


def normalize_trigger_path(path):
    """Remove trigger prefixes, such as HLT and L1, as well as suffixes, such
    as _v, from trigger paths"""
    # Remove HLT and L1 prefixes
    if path.startswith("HLT_"):
        path = path[4:]
    elif path.startswith("L1_"):
        path = path[3:]
    # Remove _v suffix
    if path.endswith("_v"):
        path = path[:-2]
    return path


def get_trigger_paths_for(dataset, is_mc, trigger_paths, trigger_order=None,
                          normalize=True, era=None):
    """Get trigger paths needed for the specific dataset.

    Parameters
    ----------
    dataset
        Name of the dataset
    trigger_paths
        dict mapping dataset names to their triggers
    trigger_order
        Datasets to define the order in which the triggers are applied.
    normalize
        Whether to remove HLT_ from the beginning
    era
        If not None and if <name>_era, where name is any dataset name, is
        present in `trigger_paths`, it will be used over just <name>. This can
        be used to define per era triggers.

    Returns
    -------
    pos_triggers
        Triggers that events in the data set need to pass
    neg_triggers
        Trigger that events in the data set must not pass (to avoid double
        counting)
    """
    if isinstance(trigger_order, dict):
        if era in trigger_order.keys():
            trigger_order = trigger_order[era]
        else:
            trigger_order = trigger_order["other"]
    pos_triggers = []
    neg_triggers = []
    if is_mc:
        for key in trigger_order:
            pos_triggers.extend(trigger_paths[key])
    else:
        for key in trigger_order:
            if (key == dataset) or (
                era is not None and key == dataset + "_" + era
            ):
                break
            neg_triggers.extend(trigger_paths[key])
        else:
            raise ValueError(f"Dataset {dataset} not in trigger_order")
        pos_triggers = trigger_paths[key]
    pos_triggers = list(dict.fromkeys(pos_triggers))
    neg_triggers = list(dict.fromkeys(neg_triggers))
    if normalize:
        pos_triggers = [normalize_trigger_path(t) for t in pos_triggers]
        neg_triggers = [normalize_trigger_path(t) for t in neg_triggers]
    return pos_triggers, neg_triggers


def get_event_files(eventdir, eventext, datasets):
    out = {}
    for dsname in datasets:
        out[dsname] = glob(os.path.join(eventdir, dsname, "*" + eventext))
    return out


def hist_split_strcat(hist):
    ret = {}
    cats = {}
    for ax in hist.axes:
        if isinstance(ax, hi.axis.StrCategory):
            cats[ax.name] = tuple(ax)
    for idx in product(*cats.values()):
        hist_idx = {ax: pos for ax, pos in zip(cats.keys(), idx)}
        ret[idx] = hist[hist_idx]
    return ret


def get_hist_cat_values(hist):
    """Return a map from the different categories of a hist histogram
    to the values (in the same way hist.values did for a coffea hists).
    Will hopefully be superseded in the near future by a dedicated hist
    function."""
    axs = [ax for ax in hist.axes if isinstance(ax, hi.axis.StrCategory)]
    if len(axs) == 0:
        return {(): hist.values()}
    cat_combs = None
    for ax in axs:
        if cat_combs is None:
            cat_combs = [((ax.name, cat),) for cat in ax]
        else:
            cat_combs = [cc + ((ax.name, cat),)
                         for cat in ax for cc in cat_combs]
    return {tuple([c[1] for c in cc]): hist[{c[0]: c[1] for c in cc}].values()
            for cc in cat_combs}


def hist_divide(num, denom):
    """Return a histogram with bin heights = num / denom and errors set
    accordinging to error propagation. The result will have zeros where num and
    denom both have zero bin heights."""
    hout = num.copy()

    num_val = num.values()
    denom_val = denom.values()
    both_zero = (num_val == 0) & (denom_val == 0)
    denom_val = np.where(both_zero, 1, denom_val)
    ratio = num_val / denom_val
    if issubclass(hout.storage_type, hi.storage.Weight):
        num_var = num.variances(flow=True)
        denom_var = denom.variances(flow=True)
        var = (num_var * denom_val ** 2
               + denom_var * num_val ** 2) / denom_val ** 4
        var[both_zero] = np.nan
        hout[:] = np.stack([ratio, var], axis=-1)
    else:
        hout[:] = ratio

    hout.label = "Ratio"

    return hout


def chunked_calls(array_param, returns_multiple=False, chunksize=10000,
                  num_threads=1):
    """A decorator that will split a function call into multiple calls on
    smaller chunks of data. This can be used to reduce peak memory usage or
    to paralellzie the call.
    Only arguments that have the attribute shape and have the same size in
    their first dimension as the value of the parameter named by array_param
    will get chunked.
    The return values of each call will be concatenated.
    The resulting functions will have two additional parameters, chunksize and
    num_threads. For a description see below.

    Parameters
    ----------
    array_param
        Parameter that will defninitely be a chunkable argument
    returns_multiple
        Needs to be set to true if the function returns more than one
        variable, e.g. as a tuple or list
    chunksize
        Default maximum chunk size to call the function on. The chunksize can
        be adjusted by using the keyword argument `chunksize` of the resulting
        function.
    num_threads
        Number of simultaneous threads. Each thread processes one chunk at a
        time, allowing to process multiple chunks in parallel and on multiple
        cores. The number of threads can be adjusted by using the keyword
        argument `num_threads` of the resulting function.

    Returns
    -------
    decorator
        Decorator function
    """

    def concatenate(arrays):
        if isinstance(arrays[0], np.ndarray):
            return np.concatenate(arrays)
        else:
            return ak.concatenate(arrays)

    def decorator(func):
        sig = inspect.signature(func)

        def do_work(kwargs, array_parameters, start, stop):
            chunked_kwargs = kwargs.copy()
            for param in array_parameters:
                chunked_kwargs[param] = kwargs[param][start:stop]
            return func(**chunked_kwargs)

        @wraps(func)
        def wrapper(*args, chunksize=chunksize, num_threads=num_threads,
                    **kwargs):
            kwargs = sig.bind(*args, **kwargs).arguments
            rows = len(kwargs[array_param])
            if rows <= chunksize:
                # Nothing to chunk, just return whatever func returns
                return func(**kwargs)
            if num_threads > 1:
                pool = ThreadPoolExecutor(max_workers=num_threads)
            array_parameters = {array_param}
            for param, arg in kwargs.items():
                if hasattr(arg, "__len__") and len(arg) == rows:
                    array_parameters.add(param)
            starts = np.arange(0, rows, chunksize)
            stops = np.r_[starts[1:], rows]
            result_funcs = []
            for start, stop in zip(starts, stops):
                if num_threads > 1:
                    result_funcs.append(pool.submit(
                        do_work, kwargs, array_parameters, start, stop).result)
                else:
                    result_funcs.append(partial(
                        do_work, kwargs, array_parameters, start, stop))
            ret_chunks = []
            for result_func in result_funcs:
                ret_chunk = result_func()
                if ret_chunk is None:
                    if num_threads > 1:
                        pool.shutdown(cancel_futures=True)
                    return None
                ret_chunks.append(ret_chunk)
                # Force clean up of memory to keep usage low
                gc.collect()
            if num_threads > 1:
                pool.shutdown()
            if len(ret_chunks) == 1:
                concated = ret_chunks[0]
            elif returns_multiple:
                # Have to transpose ret_chunks
                ret_chunks_t = [[] for _ in range(len(ret_chunks[0]))]
                for ret_chunk in ret_chunks:
                    for ret_split, chunk_val in zip(ret_chunk, ret_chunks_t):
                        chunk_val.append(ret_split)
                concated = tuple(concatenate(v) for v in ret_chunks_t)
            else:
                concated = concatenate(ret_chunks)
            return concated

        return wrapper
    return decorator


def onedimeval(func, *arrays, tonumpy=True, output_like=0):
    """Evaluate a function on the flattened one dimensional version of arrays

    Parameters
    ----------
    func
        Function to execute
    *arrays
        Arrays that will be planned and fed into ``func``
    tonumpy
        Whether to convert the flattened arrays to numpy arrays
    output_like
        Position of the array in ``arrays`` to take parameters and behavior
        from

    Returns
    -------
        An unflattened version of the array that was returned by ``func``.
        Its behavior and parameters are set according to the array pointed to
        by ``output_like``
    """
    counts_all_arrays = []
    flattened_arrays = []
    for array in arrays:
        flattened = array
        counts = []
        for i in range(flattened.ndim - 1):
            if isinstance(flattened.type.type, ak.types.RegularType):
                counts.append(flattened.type.type.size)
            else:
                counts.append(ak.num(flattened))
            flattened = ak.flatten(flattened)
        if tonumpy:
            flattened = np.asarray(flattened)
        counts_all_arrays.append(counts)
        flattened_arrays.append(flattened)
    res = func(*flattened_arrays)
    for count in reversed(counts_all_arrays[output_like]):
        res = ak.unflatten(res, count)
    for name, val in ak.parameters(arrays[output_like]).items():
        res = ak.with_parameter(res, name, val)
    res.behavior = arrays[output_like].behavior
    return res


def akremask(array, mask):
    """Make an array of length of ``mask``, where values inside ``array`` are
    used where mask is ``True`` and masked values where it is ``False``.

    Parameters
    ----------
    array
        Array to obtain the values from. Length is equal to the number of
        ``True`` occurances in ``mask``
    mask
        Array identifying masked values

    Returns
    -------
    array
        The produced array
    """
    if ak.sum(mask) != len(array):
        raise ValueError(f"Got array of length {len(array)} but mask needs "
                         f"{ak.sum(mask)}")
    if len(array) == 0:
        return ak.pad_none(array, len(mask), axis=0)
    offsets = np.cumsum(np.asarray(mask)) - 1
    return ak.mask(array[offsets], mask)


class VirtualArrayCopier:
    """Create a shallow copy of the an awkward Record Array such as NanoEvents
    while trying to not make virtual subarrays load their contents.

    Setting fields in record arrays containing virtual arrays often leads to
    the virtual arrays to be loaded, which should be avoided. The class makes
    it possible to achieve this.

    Notes
    -----
        With the removal of Virtual Arrays in Awkward version 2, this class
        will also lose its function.
    """
    def __init__(self, array, attrs=[]):
        """
        Parameters
        ----------
        array
            Record array containing virtual arrays that should not be touched
        attrs
            Attributes of the array to keep
        """
        self.data = {f: array[f] for f in ak.fields(array)}
        self.behavior = array.behavior
        self.attrs = {attr: getattr(array, attr) for attr in attrs}

    def __setitem__(self, key, value):
        """Add a field to the array"""
        self.data[key] = value

    def __getitem__(self, key):
        """Get a field from the array"""
        return self.data[key]

    def __delitem__(self, key):
        """Remove a field from the array"""
        del self.data[key]

    def get(self):
        """Get an awkward Array version of the copy"""
        array = ak.Array(self.data)
        array.behavior = self.behavior
        for attr, value in self.attrs.items():
            setattr(array, attr, value)
        return array

    def wrap_with_copy(self, func):
        """A decorator that will bind the result from ``get()`` to the first
        parameter of the function"""
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(self.get(), *args, **kwargs)
        return wrapper


def akismasked(arr):
    """Return ``True`` if arr is masked on any axis"""
    t = arr
    while hasattr(t, "type"):
        if isinstance(t.type, ak.types.OptionType):
            return True
        t = t.type
    return False
