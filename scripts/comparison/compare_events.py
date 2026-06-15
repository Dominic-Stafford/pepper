#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import awkward as ak


def get_subfolders(folder):
    return {f.name: f.path for f in os.scandir(folder) if f.is_dir()}


def get_event_files(folder):
    return [
        f.path for f in os.scandir(folder)
        if f.is_file() and f.name.endswith((".h5", ".hdf5", ".root"))
    ]


def compare_ak(a, b, path, rtol, atol, ignore):
    if path in ignore:
        return []

    errors = []

    if len(a) != len(b):
        errors.append(f"{path}: length {len(a)} vs {len(b)}")
        return errors

    fields_a = ak.fields(a)
    fields_b = ak.fields(b)

    if fields_a or fields_b:
        def field_path(f):
            return f"{path}/{f}" if path else f

        non_ignored_a = {f for f in fields_a if field_path(f) not in ignore}
        non_ignored_b = {f for f in fields_b if field_path(f) not in ignore}
        if non_ignored_a != non_ignored_b:
            errors.append(
                f"{path}: fields {sorted(non_ignored_a)} vs {sorted(non_ignored_b)}"
            )
            return errors

        for field in sorted(non_ignored_a):
            errors.extend(compare_ak(a[field], b[field], field_path(field), rtol, atol, ignore))
        return errors

    # Jagged array: check inner structure before flattening
    if a.ndim > 1:
        counts_a = ak.num(a, axis=1)
        counts_b = ak.num(b, axis=1)
        if not ak.all(counts_a == counts_b):
            errors.append(f"{path}: jagged structure mismatch")
            return errors
        errors.extend(compare_ak(
            ak.flatten(a, axis=1), ak.flatten(b, axis=1),
            f"{path}[*]", rtol, atol, ignore,
        ))
        return errors

    # 1D leaf array
    try:
        a_np = ak.to_numpy(a)
        b_np = ak.to_numpy(b)
    except Exception:
        # String or other non-numpy-convertible type
        try:
            if not ak.all(a == b):
                errors.append(f"{path}: value mismatch")
        except Exception as e:
            errors.append(f"{path}: comparison failed: {e}")
        return errors

    if np.issubdtype(a_np.dtype, np.floating):
        bad = ~np.isclose(a_np, b_np, rtol=rtol, atol=atol, equal_nan=True)
    else:
        bad = a_np != b_np

    if bad.any():
        i = int(np.where(bad)[0][0])
        errors.append(
            f"{path}: {bad.sum()} mismatch(es), first at [{i}]: "
            f"{a_np[i]} vs {b_np[i]}"
        )

    return errors


def compare_hdf5(path1, path2, rtol, atol, ignore):
    from pepper import HDF5File
    errors = []
    with HDF5File(path1, "r") as f1, HDF5File(path2, "r") as f2:
        keys1 = {k for k in f1.keys() if k not in ignore}
        keys2 = {k for k in f2.keys() if k not in ignore}
        if keys1 != keys2:
            errors.append(f"key mismatch: {sorted(keys1)} vs {sorted(keys2)}")
            return errors
        for key in sorted(keys1):
            if isinstance(f1[key], ak.Array):
                errors.extend(compare_ak(f1[key], f2[key], key, rtol, atol, ignore))
    return errors


def compare_root(path1, path2, rtol, atol, ignore):
    import uproot
    errors = []
    with uproot.open(path1) as f1, uproot.open(path2) as f2:
        keys1 = {k for k in f1.keys(cycle=False, recursive=True) if k not in ignore}
        keys2 = {k for k in f2.keys(cycle=False, recursive=True) if k not in ignore}
        if keys1 != keys2:
            errors.append(f"tree mismatch: {sorted(keys1)} vs {sorted(keys2)}")
            return errors
        for key in sorted(keys1):
            if isinstance(f1[key], uproot.TTree):
                a = f1[key].arrays(library="ak")
                b = f2[key].arrays(library="ak")
                errors.extend(compare_ak(a, b, key, rtol, atol, ignore))
    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Compare event output files between two directories"
    )
    parser.add_argument("folder1", help="First output directory")
    parser.add_argument("folder2", help="Second output directory")
    parser.add_argument(
        "--rtol", type=float, default=1e-5,
        help="Relative tolerance for float comparison (default: 1e-5)",
    )
    parser.add_argument(
        "--atol", type=float, default=0.0,
        help="Absolute tolerance for float comparison (default: 0, exact)",
    )
    parser.add_argument(
        "--ignore", metavar="FIELD", action="append", default=[],
        help="Field path to ignore, e.g. events/nJet (can be given multiple times)",
    )
    args = parser.parse_args()

    ignore = set(args.ignore)

    subs1 = get_subfolders(args.folder1)
    subs2 = get_subfolders(args.folder2)

    if set(subs1) != set(subs2):
        only1 = sorted(set(subs1) - set(subs2))
        only2 = sorted(set(subs2) - set(subs1))
        msg = "Sample mismatch:"
        if only1:
            msg += f"\n  Only in {args.folder1}: {only1}"
        if only2:
            msg += f"\n  Only in {args.folder2}: {only2}"
        print(msg, file=sys.stderr)
        sys.exit(1)

    all_errors = []

    for sample in sorted(subs1):
        files1 = get_event_files(subs1[sample])
        files2 = get_event_files(subs2[sample])

        if len(files1) != len(files2):
            all_errors.append(
                f"Sample '{sample}': file count mismatch "
                f"({len(files1)} vs {len(files2)})"
            )
            continue

        if len(files1) == 0:
            print(f"Warning: sample '{sample}' has no event files", file=sys.stderr)
            continue

        if len(files1) > 1:
            print(
                f"Warning: sample '{sample}' has {len(files1)} files, "
                "only comparing first",
                file=sys.stderr,
            )

        f1, f2 = files1[0], files2[0]
        ext1, ext2 = os.path.splitext(f1)[1], os.path.splitext(f2)[1]

        if ext1 != ext2:
            all_errors.append(
                f"Sample '{sample}': file type mismatch ({ext1} vs {ext2})"
            )
            continue

        if ext1 in (".h5", ".hdf5"):
            diffs = compare_hdf5(f1, f2, args.rtol, args.atol, ignore)
        elif ext1 == ".root":
            diffs = compare_root(f1, f2, args.rtol, args.atol, ignore)
        else:
            all_errors.append(f"Sample '{sample}': unknown file type '{ext1}'")
            continue

        for diff in diffs:
            all_errors.append(f"Sample '{sample}': {diff}")

    if all_errors:
        print("Differences found:", file=sys.stderr)
        for err in all_errors:
            print(f"  {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
