#!/usr/bin/env python3

import os
from collections.abc import Mapping
import json
from itertools import product
import numpy as np
import uproot
import coffea.util
import hist as hi


class _JSONEncoderWithSets(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


class HistCollection(Mapping):
    """Class that provides access to and writing of a multiple histograms. The
       histograms can be in stored in different formats (hist and root)."""

    @staticmethod
    def root_key(idx):
        return "/".join(idx) + "/hist"

    @classmethod
    def rootdir_to_hist(cls, rootdir, infokey="pepper_hist_info"):
        """Convert an uproot directory containing histograms to a single hist
        histogram using the information given in infokey to merge the
        histograms along category axes"""
        info = json.loads(rootdir[infokey])
        order = info["axis_order"]
        axes = info["cat_axes"]
        cat_axes = [hi.axis.StrCategory(**axes[k], growth=True) for k in order]
        # Obtain dense axes from an arbitrary histogram in the rootdir
        for key, obj in rootdir.iteritems(cycle=False):
            if ("/" + key).endswith("/hist"):
                dense_axes = obj.to_hist().axes
                break
        else:
            raise RuntimeError("No hists in root directory")
        axes = cat_axes + list(dense_axes)
        val = np.zeros(tuple(ax.extent for ax in axes))
        var = np.zeros_like(val)
        key_to_idx = [dict(zip(ax, range(len(ax)))) for ax in cat_axes]
        cats_present = set()
        for key, obj in rootdir.iteritems(cycle=False):
            key = "/" + key
            if not key.endswith("/hist"):
                continue
            key = key.split("/")[1:-1]
            idx = tuple(m[k] for m, k in zip(key_to_idx, key))
            val[idx] = obj.values(flow=True)
            var[idx] = obj.variances(flow=True)
            cats_present.add(tuple(key))
        data = np.stack([val, var], axis=-1)
        hist = hi.Hist(
            *axes,
            name=info["name"],
            label=info["label"],
            storage=hi.storage.Weight(),
            data=data
        )
        hist.pepper_cats_present = cats_present
        return hist

    @staticmethod
    def get_cats_present(hist):
        cataxes = [
            ax for ax in hist.axes if isinstance(ax, hi.axis.StrCategory)
        ]
        return set(product(*cataxes))

    @staticmethod
    def hist_split_strcat(hist, cats_present=None):
        """Split hist along its StrCategory axes, if any, and return a dict of
        the splits."""
        ret = {}
        cats = {}
        for ax in hist.axes:
            if isinstance(ax, hi.axis.StrCategory):
                cats[ax.name] = tuple(ax)
        for idx in product(*cats.values()):
            hist_idx = {ax: pos for ax, pos in zip(cats.keys(), idx)}
            if tuple(hist_idx.values()) not in cats_present:
                continue
            ret[idx] = hist[hist_idx]
        return ret

    @staticmethod
    def get_hist_info(hist):
        """Return a dictionary containing information about the hist so that
        after saving all its values into a Root file, it can still be
        reconstructed as a hist histogram"""
        cat_axes = {}
        axis_order = []
        for ax in hist.axes:
            if not isinstance(ax, hi.axis.StrCategory):
                continue
            if ax.name in cat_axes:
                cat_axes[ax.name]["categories"] |= set(ax)
            else:
                cat_axes[ax.name] = {
                    "name": ax.name,
                    "label": ax.label,
                    "categories": set(ax)
                }
                axis_order.append(ax.name)
        info = {
            "name": hist.name,
            "label": hist.label,
            "cat_axes": cat_axes,
            "axis_order": axis_order
        }
        return info

    def __init__(self, path, key_fields, content=None, userdata=None):
        self.path = path
        self.key_fields = key_fields
        self._content = {} if content is None else content.copy()
        self.userdata = userdata

    @classmethod
    def from_json(cls, fileobj, path=None):
        if path is None:
            path = os.path.dirname(os.path.realpath(fileobj.name))
        data = json.load(fileobj)
        key_fields = data["key_fields"]
        content = {tuple(k): v for k, v in zip(*data["content"])}
        userdata = data["userdata"]
        return cls(
            path=path,
            key_fields=key_fields,
            content=content,
            userdata=userdata
        )

    @classmethod
    def from_single_hist(cls, histpath):
        path = os.path.dirname(os.path.realpath(histpath))
        content = {(): os.path.basename(histpath)}
        return cls(
            path=path,
            key_fields=[],
            content=content
        )

    def __len__(self):
        return len(self._content)

    def __iter__(self):
        return iter(self._content)

    def __getitem_tuple__(self, key):
        if len(key) > len(self.key_fields):
            raise KeyError(key)
        key = key + (None,) * (len(self.key_fields) - len(key))
        NoneType = type(None)
        unspec = []
        for i, field in enumerate(self.key_fields):
            if isinstance(key[i], (NoneType, list)):
                unspec.append((i, field))
        if len(unspec) == 0:
            return self._content[key]
        keys = self._content.keys()
        for i, keypart in enumerate(key):
            if len(keys) == 0:
                break
            testkeys = keys
            keys = []
            if keypart is None:
                keys = testkeys
                continue
            if isinstance(keypart, list):
                for testkey in testkeys:
                    if testkey[i] in keypart:
                        keys.append(testkey)
                continue
            for testkey in keys:
                if testkey[i] == keypart:
                    keys.append(testkey)

        if len(keys) == 0:
            raise KeyError(key)
        content = {}
        for k in keys:
            content[tuple(k[i] for i, f in unspec)] = self._content[k]
        fields = [f for i, f in unspec]
        ret = self.__class__(
            path=self.path,
            key_fields=fields,
            content=content
        )
        return ret

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self.__getitem_tuple__(key)
        elif isinstance(key, dict):
            bad_fields = key.keys() - set(self.key_fields)
            if len(bad_fields) > 0:
                raise ValueError(
                    f"Invalid key field. Expected any of {self.key_fields}, "
                    f"but key included: {', '.join(bad_fields)}"
                )
            tuplekey = tuple(key.get(f, None) for f in self.key_fields)
            try:
                return self[tuplekey]
            except KeyError as e:
                raise KeyError(key) from e
        else:
            raise TypeError(
                "HistCollection indices are either tuples or dicts, not "
                f"{type(key)}")

    def __add__(self, other):
        if isinstance(other, HistCollection):
            content = self._content.copy()
            content.update(other._content)
            return self.__class__(self.path, self.key_fields, content)
        else:
            raise ValueError("Can only add HistCollection")

    def __iadd__(self, other):
        if isinstance(other, HistCollection):
            self._content.update(other._content)
            return self
        else:
            raise ValueError("Can only add HistCollection")

    def copy(self):
        return self.__class__(
            path=self.path,
            key_fields=self.key_fields,
            content=self._content.copy()
        )

    def load(self, key):
        path = self[key]
        if isinstance(path, self.__class__):
            return {k: path.load(k) for k, v in path.keys()}
        if not os.path.isabs(path):
            path = os.path.join(self.path, path)
        if path.endswith(".root"):
            with uproot.open(path) as f:
                return self.rootdir_to_hist(f)
        else:
            return coffea.util.load(path)

    def items_loaded(self):
        for key in self.keys():
            yield key, self.load(key)

    def save(
        self, key, hist, filename, format, cats_present,
        infokey="pepper_hist_info"
    ):
        if len(key) != len(self.key_fields):
            raise ValueError(
                f"Invalid key length, expected {len(self.key_fields)}, "
                f"but got {len(key)}")
        if (
            not isinstance(key, tuple)
            or not all(isinstance(k, str) for k in key)
        ):
            raise ValueError("key must be a tuple of str")
        if os.path.isabs(filename):
            filepath = filename
        else:
            filepath = os.path.join(self.path, filename)
        if format == "hist":
            hist.pepper_cats_present = cats_present
            coffea.util.save(hist, filepath)
        elif format == "root":
            with uproot.recreate(filepath) as f:
                for idx, hist_split in self.hist_split_strcat(
                        hist, cats_present).items():
                    f[self.root_key(idx)] = hist_split
                info = self.get_hist_info(hist)
                f[infokey] = json.dumps(info, cls=_JSONEncoderWithSets)
        else:
            raise ValueError(
                f"format needs to be 'root' or 'hist', not {format}")
        self._content[key] = filename

    def save_metadata_json(self, fileobj):
        keys = [list(k) for k in self._content.keys()]
        vals = [v for v in self._content.values()]
        data = {
            "key_fields": list(self.key_fields),
            "content": [keys, vals],
            "userdata": self.userdata,
        }
        json.dump(data, fileobj, indent=4)
