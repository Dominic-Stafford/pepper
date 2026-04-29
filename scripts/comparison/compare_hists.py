import argparse
import numpy as np
import hist
import itertools
import sys

from pepper import HistCollection

parser = argparse.ArgumentParser()
parser.add_argument("json1")
parser.add_argument("json2")
parser.add_argument("--rtol", type=float, default=1e-5)
parser.add_argument("--atol", type=float, default=1e-3)
args = parser.parse_args()

with open(args.json1) as f:
    coll1 = HistCollection.from_json(f)
with open(args.json2) as f:
    coll2 = HistCollection.from_json(f)

assert all(k in coll2 for k in coll1.keys())
assert all(k in coll1 for k in coll2.keys())

for k in coll1.keys():
    h1 = coll1.load(k)
    h2 = coll2.load(k)

    assert all(ax.name in h2.axes.name for ax in h1.axes)
    assert all(ax.name in h1.axes.name for ax in h2.axes)

    cat_axes = [ax for ax in h1.axes
                if isinstance(ax, hist.axis.StrCategory)
                or isinstance(ax, hist.axis.IntCategory)]
    for ax1 in cat_axes:
        ax2 = h2.axes[ax1.name]
        assert all(key in ax1 for key in ax2)
        assert all(key in ax2 for key in ax1)

    for cat_vals in itertools.product(*cat_axes):
        cat_dict = dict(zip([ax.name for ax in cat_axes], cat_vals))

        h1cat = h1[cat_dict]
        h2cat = h2[cat_dict]

        if not np.all(np.isclose(h1cat.values(flow=True), h2cat.values(flow=True),
                                 rtol=args.rtol, atol=args.atol)):
            print(f"Values differ for {k[0]} {k[1]} ({','.join(cat_vals)}):")
            print(">> ", h1cat.values(flow=True))
            print("<< ", h2cat.values(flow=True))
            print("reldiff ", np.abs(h1cat.values(flow=True)/h2cat.values(flow=True) - 1.))
            sys.exit(-1)
        if not np.all(np.isclose(h1cat.variances(flow=True), h2cat.variances(flow=True),
                                 rtol=2*args.rtol, atol=args.atol**2)):
            print(f"Variances differ for {k[0]} {k[1]} ({','.join(cat_vals)}):")
            print(">> ", h1cat.variances(flow=True))
            print("<< ", h2cat.variances(flow=True))
            print("reldiff ", np.abs(h1cat.variances(flow=True)/h2cat.variances(flow=True) - 1.))
            sys.exit(-1)
