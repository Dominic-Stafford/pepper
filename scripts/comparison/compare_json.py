#!/usr/bin/env python3
import argparse
import json
import math
import sys


def compare(a, b, path, atol, rtol, diffs):
    if isinstance(a, bool) or isinstance(b, bool):
        if a != b:
            diffs.append(f"{path}: {a!r} != {b!r}")
        return

    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not (math.isclose(a, b, rel_tol=rtol, abs_tol=atol)):
            diffs.append(f"{path}: {a} != {b} (atol={atol}, rtol={rtol})")
        return

    if type(a) is not type(b):
        diffs.append(f"{path}: type mismatch ({type(a).__name__} vs {type(b).__name__})")
        return

    if isinstance(a, dict):
        keys_a, keys_b = set(a.keys()), set(b.keys())
        for k in sorted(keys_a - keys_b):
            diffs.append(f"{path}.{k}: key only in first file")
        for k in sorted(keys_b - keys_a):
            diffs.append(f"{path}.{k}: key only in second file")
        for k in sorted(keys_a & keys_b):
            compare(a[k], b[k], f"{path}.{k}", atol, rtol, diffs)

    elif isinstance(a, list):
        if len(a) != len(b):
            diffs.append(f"{path}: list length differs ({len(a)} vs {len(b)})")
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{path}[{i}]", atol, rtol, diffs)

    elif a != b:
        diffs.append(f"{path}: {a!r} != {b!r}")


def main():
    parser = argparse.ArgumentParser(description="Compare two JSON files.")
    parser.add_argument("file1", help="First JSON file")
    parser.add_argument("file2", help="Second JSON file")
    parser.add_argument("--atol", type=float, default=0.0,
                        help="Absolute tolerance for numeric comparisons (default: 0)")
    parser.add_argument("--rtol", type=float, default=1e-9,
                        help="Relative tolerance for numeric comparisons (default: 1e-9)")
    args = parser.parse_args()

    with open(args.file1) as f:
        data1 = json.load(f)
    with open(args.file2) as f:
        data2 = json.load(f)

    diffs = []
    compare(data1, data2, "root", args.atol, args.rtol, diffs)

    if diffs:
        print(f"Files differ ({len(diffs)} difference(s)):")
        for d in diffs:
            print(f"  {d}")
        sys.exit(1)


if __name__ == "__main__":
    main()
