import numpy as np
import awkward as ak
from collections import defaultdict
import logging
from copy import copy
import pepper.misc


logger = logging.getLogger(__name__)


class Selection:
    """Keep track of weighted cuts

    Attributes
    ----------
    names
        Name of the cuts
    cuts
        Record array of the cut weights
    """
    def __init__(self):
        self.names = []
        self.cuts = ak.Array({})

    def all(self, names=None):
        """Get the product of all cuts' weights

        Paramters
        ---------
        names
            Only use the cuts named in this list

        Returns
        -------
        total
            Product of all cuts
        """
        if names is None:
            names = self.names
        total = None
        for name in names:
            if total is None:
                total = np.array(self.cuts[name])
            else:
                total = total * self.cuts[name]
        return total

    def add_cut(self, name, accept):
        """Add a cut to the selection

        Parameters
        ----------
        name
            Name of the cut
        accept
            Array whether an event passes a cut. Unweighted cuts can be
            specified as bools and weighted cuts as floats
        """
        self.cuts[name] = accept
        self.names.append(name)

    def clear(self):
        """Remove all cuts from the selection"""
        self.names = []
        self.cuts = ak.Array({})

    def __len__(self):
        """Number of cuts"""
        return len(self.names)

    def __copy__(self):
        """Copy, also creating copies of the class attributes"""
        s = self.__class__.__new__(self.__class__)
        s.__dict__.update(self.__dict__)
        s.names = copy(self.names)
        s.cuts = copy(self.cuts)
        return s


class Selector:
    """Keeps track of the current event selection and data

    Attributes
    ----------
    data
        Event data for events passing the applied cuts
    metadata
        Event meta data such as file name
    systematics
        Record array of systematic weights
    cutnames
        Names of applied and possibly unapplied cuts
    cut_systematic_map
        Tells what systematic weights were added at which cut
    done_steps
        Tells which steps (cuts and column additions) where performed
    rng
        A random number generator that was seeded according to ``rng_seed``
        in the class constructor and that should be used throughout the
        selection
    cats
        Currently known categorizations. Maps names of categorizations to their
        categories.
    """

    def __init__(self, data, weight=None, on_update=None, applying_cuts=True,
                 rng_seed=None, output_filler=None):
        """
        Parameters
        ----------
        data
            Array holding the events' data
        weight
            One dimensional array of size equal to ``data`` size, describing
            the events' weight. ``None`` if unweighted events
        on_update
            callable or list of callables that get called after a
            call to ``add_cut`` or ``set_column``. The callable should accept
            the keyword argument data, systematics and cut.
        applying_cuts
            Wether to apply cuts added with ``add_cut``. If ``False`, cuts will
            be kept at ``unapplied_cuts``
        rng_seed
            int or tuple of ints to seed the random number generator
            with. If ``None``, a random seed will be used. For details
            see the parameter of ``numpy.random.default_rng()``.
        output_filler
            A direct link to the output_filler, currently only used to
            access the accumulators for normalisation of systematics
        """
        self.data = data
        if hasattr(self.data, "metadata"):
            self.metadata = self.data.metadata
        if on_update is None:
            self.on_update = []
        elif isinstance(on_update, list):
            self.on_update = on_update.copy()
        else:
            self.on_update = [on_update]

        if weight is not None:
            self.systematics = ak.Array({"weight": weight})
        else:
            self.systematics = None

        self.cutnames = []
        self.unapplied_cuts = Selection()
        self.cut_systematic_map = defaultdict(list)
        self.done_steps = set()

        self.rng = np.random.default_rng(rng_seed)
        self.cats = {}
        self.filler = output_filler

        self._applying_cuts = True
        self.add_cut("BeforeCuts", np.full(self.num, True))
        self._applying_cuts = applying_cuts

    @property
    def applying_cuts(self):
        return self._applying_cuts

    @applying_cuts.setter
    def applying_cuts(self, val):
        if val and not self._applying_cuts and len(self.unapplied_cuts) > 0:
            self.apply_all_cuts()
        self._applying_cuts = val

    @property
    def unapplied_product(self):
        """An array giving the effect on the event weight of all cuts added
        after the last cut was applied."""
        if len(self.unapplied_cuts) == 0:
            return np.full(self.num, 1.)
        else:
            return self.unapplied_cuts.all()

    @property
    def final(self):
        """Data of events which have passed all cuts (including unapplied ones)
        """
        return ak.mask(
            self.data, ak.values_astype(self.unapplied_product, bool))

    @property
    def final_systematics(self):
        """Systematics of events which have passed all cuts
        (including unapplied ones)"""
        if self.systematics is None:
            return None
        unapplied = self.unapplied_product
        masked = ak.mask(
            self.systematics, ak.values_astype(unapplied, bool))
        masked["weight"] = masked["weight"] * unapplied
        return masked

    @property
    def num(self):
        "Number of events passing the applied cuts"
        return len(self.data)

    @property
    def num_final(self):
        "Number of selected events passing all cuts (including unapplied ones)"
        if len(self.unapplied_cuts) == 0:
            return self.num
        else:
            return ak.sum(self.unapplied_product != 0)

    def set_cat(self, name, categories, safe=True):
        """Set a categorisation

        A categorisation is for example the decay channel of an event.
        Categories inside this categorisation would be the ee-decay or the
        µµ-decay.

        Parameters
        ----------
        name
            Name of the categorization
        categories
            Categories inside the categorisation. Each category must also be
            a field in ``data``
        safe --  Boolean on whether to check if all events belong to exactly
                 one category. Note that if this is not satisfied, the
                 cutflows will be incorrect, as will histograms when
                 integrating over the category axis, as events will either be
                 silently cut away (if not included in any catgories) or
                 double counted (if appearing in multiple categories)
        """
        if safe:
            selected = np.zeros(len(self.data), dtype=int)
            for cat in categories:
                selected = selected + self.data[cat]
            if not ak.all(selected == 1):
                raise ValueError(
                    "Not all events fit into exactly one category")
        for exist_categories in self.cats.values():
            if exist_categories & categories:
                raise ValueError(
                    f"The following category(s) are already used in "
                    f"a previous categorisation, please rename: "
                    f"{exist_categories & categories}")
        self.cats[name] = categories

    def _invoke_callbacks(self):
        """Execute the callbacks in ``on_update`` on the final data"""
        data = self.final
        systematics = self.final_systematics
        for cb in self.on_update:
            cb(data=data, systematics=systematics, cut=self.cutnames[-1],
               done_steps=self.done_steps, cats=self.cats)

    def _get_category_mask(self, categories):
        """Get a mask that is only true if an event is contained in all
        categories named by ``categories``"""
        num = self.num
        mask = np.full(num, True)
        for cat, regs in categories.items():
            if cat not in self.cats:
                raise ValueError(f"Unknown categorization {cat}")
            cat_mask = np.full(num, False)
            for reg in regs:
                cat_mask = cat_mask | self.data[reg]
            mask = mask & cat_mask
        return mask

    def add_cut(self, name, accept, systematics=None, no_callback=False,
                categories=None):
        """Adds a cut and applies it if `self.applying_cuts` is True, otherwise
        the cut will be stored in `self.unapplied_cuts`. Applying in this
        context means that rows of `self.data` are discarded accordingly.

        Parameters
        ----------
        name
            Name of the cut
        accept
            An array of bools or floats or a tuple of the former and a
            systematics dict or a callable returning any of the former.
            In the array a value of ``0`` or ``False`` means that the
            event corresponding to the row is discarded. Any other value
            will get multiplied into the event weight. Here a value of
            ``True`` corresponds to a ``1``.
            The systematics dict is a mapping of systematics name ->
            values, where name and values have the same meaning as in
            ``self.set_systematic``. The values arrays of lengths equal
            ``self.num`` either before or after the cut is applied.
            Systematics given for cut evets are ignored.
            In the case this is a callable, the callable will be called and
            its return value will be used as the new value for this
            parameter.
        systematics
            A dict of name and values and has the same effect has
            calling ``self.set_systematic`` on every item. Will be
            ignored if a systematics dict is given with ``accept``.
        no_callback
            A bool whether not to call the callbacks, which usually
            fill histograms etc.
        categories
            If not None, ignore events that are not part of any of
            the specified categorizations. This is done by specifying
            a dict of lists. A key gives the name of the
            categorization, while the list contains field names
            of ``self.data``. These fields needs to be flat bool
            arrays. An event is considered to be part of a
            categorization if any of the fields are True.
        """
        def pad_cats(mask, arr):
            full_arr = np.full(self.num, 1, dtype=arr.dtype)
            full_arr[mask] = arr
            return full_arr
        if name in self.cutnames:
            raise ValueError(f"A cut with the name '{name}' "
                             "already exists")
        logger.info(f"Adding cut '{name}'"
                    + (" (no callback)" if no_callback else ""))
        if categories is not None:
            cats_mask = self._get_category_mask(categories)
            if callable(accept):
                accept = accept(self.data[cats_mask])
        elif callable(accept):
            accept = accept(self.data)
        if isinstance(accept, tuple):
            accept, systematics = accept
        # Allow accept to be masked. Treat masked values as False
        accept = ak.fill_none(accept, 0)
        self.cutnames.append(name)
        if not isinstance(accept, np.ndarray):
            accept = np.array(accept)
        if categories is not None:
            accept = pad_cats(cats_mask, accept)
        if accept.dtype == bool:
            accept_weighted = accept.astype(float)
        else:
            accept_weighted = accept
            accept = accept != 0
        self.unapplied_cuts.add_cut(name, accept_weighted)
        if self.applying_cuts:
            self.apply_all_cuts()
        if systematics is not None:
            for sysname, values in systematics.items():
                if isinstance(values, list):
                    raise ValueError("Multiple systematic variations need to "
                                     "be given as tuple, not list")
                if not isinstance(values, tuple):
                    values = (values,)
                values_old = values
                values = []
                n = self.num
                for value_old in values_old:
                    if not isinstance(value_old, np.ndarray):
                        value_old = np.array(value_old)
                    if categories is not None:
                        value_old = pad_cats(cats_mask, value_old)
                    if len(value_old) != n:
                        if self.applying_cuts:
                            value = value_old[accept]
                        else:
                            value = pad_cats(accept, value_old)
                    elif not self.applying_cuts:
                        value = ak.mask(value_old, accept)
                    else:
                        value = value_old
                    values.append(value)
                self.set_systematic(sysname, *values, cut=name)
        self.done_steps.add("cut:" + name)
        if not no_callback:
            self._invoke_callbacks()

    def apply_all_cuts(self):
        """Apply all unapplied cuts, discarding rows of ``self.data`` where
        the resulting weight is ``0`` and modify the event weight
        accordingly."""
        weighted = self.unapplied_product
        mask = weighted != 0
        self.data = self.data[mask]
        if self.systematics is not None:
            self.systematics["weight"] = self.systematics["weight"] * weighted
            self.systematics = self.systematics[mask]
        self.unapplied_cuts.clear()

    def set_systematic(self, name, *values, scheme=None, cut=None,
                       norm_post=False):
        """Set the systematic weights for an uncertainty. These will be
        found in the `self.systematics`.

        Notes
        -----
        Please make sure to provide your systematic weights in form of a ratio
        over the nominal weight.

        Parameters
        ----------
        name
            Name of the systematic
        values
            Each of these arrays gives the ratio of a systematic
            variation and the central value of the event weight
        scheme
            One of 'updown', 'numeric', 'single' or None. Determines the
            column (of the systematics table) the values will appear in.
            'updown': Requires ``values`` to have length 2. Column will
            be ``name`` + _up and _down
            'numeric' column will be ``name`` + _i where i is determined by
            enumerating ``values``
            'single': Requires ``values`` to have length 1. Column will be
            ``name``
            None: scheme will be decided based on ``values`` length, where
            numeric will be used for lengths > 2
        cut
            Name of the cut after which the systematic needs to be accounted
            for. If not None, a corresponding item will be found in
            ``self.cut_systematic_map``
        norm_post
            Boolean for whether this systematic should be normailsed at the end
            of the run, will fill the sumw into the output filler, triggering
            normalisation when saving histograms
        """

        if name == "weight":
            raise ValueError("The name of a systematic can't be 'weight'")
        if len(values) == 0 and callable(values[0]):
            values = values(self.data)
        if scheme is None:
            if len(values) == 1:
                scheme = "single"
            elif len(values) == 2:
                scheme = "updown"
            else:
                scheme = "numeric"
        if scheme == "single":
            names = [name]
        elif scheme == "updown":
            names = [f"{name}_up", f"{name}_down"]
        elif scheme == "numeric":
            names = [f"{name}_{i}" for i in range(len(values))]
        else:
            raise ValueError("scheme needs to be either 'updown', 'numeric',"
                             f"'single' or None, got {scheme}")
        for name, value in zip(names, values):
            if not isinstance(value, np.ndarray):
                value = np.array(value)
            if value.ndim != 1:
                raise ValueError("Systematics must be flat 1-d arrays!")
            self.systematics[name] = value
            if cut is not None:
                self.cut_systematic_map[cut].append(name)
        if norm_post:
            if len(self.cutnames) > 1:
                raise RuntimeError(
                    "Systematics to be normed should be set before any cuts "
                    "to account for acceptance cuts. Please move the call "
                    "to selector.systematic earlier in your processor")
            if self.filler is None:
                raise UnboundLocalError(
                    "A call was made to selctor.set_systematic with "
                    "norm_post=True, but no filler was set for this selector")
            self.filler.add_sumws_for_norm(self.systematics, names)

    def _mask(self, column, mask):
        num_events = ak.sum(mask)
        if len(column) != num_events:
            raise ValueError(f"Column length exptected to be {num_events} but "
                             f"got {len(column)}")
        if len(column) > 0:
            return ak.mask(column[np.cumsum(np.asarray(mask)) - 1], mask)
        else:
            return ak.pad_none(column, len(mask), axis=0)

    def set_column(self, column_name, column, all_cuts=False,
                   no_callback=False, lazy=False, categories=None):
        """Sets a column of ``self.data``

        Parameters
        ----------
        column_name
            Name of the column to set
        column
            Column data or a callable that returns it.
            The callable that will be called with ``self.data`` as argument
        all_cuts
            The callable from ``column`` will be called only on events
            passing all cuts (including unapplied cuts).
        no_callback
            Whether not to call the callbacks. These callbacks usually fill
            histograms etc.
        lazy
            If True, column must be a callable and the column will be
            inserted as a virtual array, making the callable only called
            when the data array determines it has to.
        categories
            If not ``None``, ignore events that are not part of any of
            the specified categorizations. This is done by specifying
            a dict of lists. A key gives the name of the
            categorization, while the list contains field names
            of ``self.data``. These fields needs to be flat bool
            arrays. An event is considered to be part of a
            categorization if any of the fields are ``True``.
        """

        logger.info(
            f"Setting column {column_name}" + (" lazily" if lazy else ""))
        if all_cuts and not self.applying_cuts:
            data = self.final
            mask = ~ak.is_none(data)
            data = ak.flatten(self.final, axis=0)
        else:
            data = self.data
            mask = None

        if categories is not None:
            cats_mask = self._get_category_mask(categories)
            data = data[cats_mask]

        if callable(column):
            if lazy:
                column = pepper.misc.VirtualArrayCopier(data).wrap_with_copy(
                    column)
                column = ak.virtual(column, cache={}, length=len(data))
            else:
                column = column(data)
        if column is None:
            raise ValueError("column content must be an array, not None")
        if categories is not None:
            column = self._mask(column, cats_mask)
        if mask is not None:
            column = self._mask(column, mask)
        if lazy:
            data_copier = pepper.misc.VirtualArrayCopier(self.data)
            data_copier[column_name] = column
            self.data = data_copier.get()
        else:
            self.data[column_name] = column
        self.done_steps.add("column:" + column_name)
        if not no_callback:
            self._invoke_callbacks()

    def set_multiple_columns(self, columns, all_cuts=False, no_callback=False,
                             categories=None):
        """Sets multiple columns of ``self.data`` at once.

        Parameters
        ----------
        columns
            A dict of column names and data or a callable returning
            the former. The callable will be called with `self.data`
            as argument.
        all_cuts
            The callables from `columns` will be called only on events
            passing all cuts (including unapplied cuts).
        no_callback
            A bool whether not to call the callbacks, which usually
            fill histograms etc.
        categories
            If not ``None``, ignore events that are not part of any of
            the specified categorizations. This is done by specifying
            a dict of lists. A key gives the name of the
            categorization, while the list contains field names
            of ``self.data``. These fields needs to be flat bool
            arrays. An event is considered to be part of a
            categorization if any of the fields are ``True``.
        """
        if callable(columns):
            if all_cuts and not self.applying_cuts:
                data = self.final
                data = ak.flatten(self.final, axis=0)
            else:
                data = self.data
            if categories is not None:
                cats_mask = self._get_category_mask(categories)
                data = data[cats_mask]
            columns = columns(data)
        if isinstance(columns, ak.Array):
            columns = {k: columns[k] for k in ak.fields(columns)}
        for name, column in columns.items():
            self.set_column(name, column, all_cuts, no_callback=True,
                            categories=categories)
        if not no_callback:
            self._invoke_callbacks()

    def get_cuts(self):
        """Returns a tuple a list of all cut names and an array containing the
           currently unapplied cuts"""
        return self.cutnames, self.unapplied_cuts.cuts

    def __copy__(self):
        """Makes a shallow copy excluding some constituents. This allows to
        work with the copy without modifying the original"""
        s = self.__class__.__new__(self.__class__)
        s.__dict__.update(self.__dict__)
        # copy(self.data) loads all fields. Workaround
        s.data = pepper.misc.VirtualArrayCopier(self.data).get()
        if self.systematics is not None:
            s.systematics = copy(self.systematics)
        s.cutnames = copy(self.cutnames)
        s.unapplied_cuts = copy(self.unapplied_cuts)
        s.cut_systematic_map = copy(self.cut_systematic_map)
        s.done_steps = copy(self.done_steps)
        s.cats = {k: v.copy() for k, v in self.cats.items()}
        return s

    def copy(self):
        """Shorthard for `copy.copy(self)``, see `__copy__`"""
        return copy(self)
