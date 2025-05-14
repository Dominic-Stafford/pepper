from tqdm import tqdm
import os
import argparse

from pepper import HistCollection

parser = argparse.ArgumentParser(
    description="""Merge output histograms from different pepper
runs. Can be used to combine runs with different datasets, or
to merge the same datasets across multiple eras.
Note that categorizations and cuts need to be identical.
For now only saving hist histograms are supported.""")
parser.add_argument("output_folder", type=str,
                    help="Output folder for the merged histograms.")
parser.add_argument("hist_jsons", nargs="+",
                    help="The hists.json files to merge.")

args = parser.parse_args()

out_folder = args.output_folder

if not os.path.isdir(out_folder):
    os.makedirs(out_folder)

jsons = args.hist_jsons

histcolls = {}
for j in jsons:
    with open(j) as f:
        histcolls[j] = HistCollection.from_json(f)


allkeys = set(k for hc in histcolls.values() for k in hc.keys())

key_fields = histcolls[jsons[0]].key_fields
userdata = histcolls[jsons[0]].userdata
if userdata is not None and "cuts" in userdata:
    userdata["cuts"] = list(set(k[0] for k in allkeys))

histcol_out = HistCollection(out_folder, key_fields, userdata=userdata)

for key in tqdm(allkeys):
    hists = []
    filenames = []
    for j, hc in histcolls.items():
        if key not in hc:
            print(f"Histogram {key} missing in json {j}. Skipping")
            continue
        filenames.append(hc[key])
        hists.append(hc.load(key))
    hsum = sum(hists)
    histcol_out.save(key, hsum, filenames[0], 'hist', None)

with open(os.path.join(out_folder, "hists.json"), "w") as f:
    histcol_out.save_metadata_json(f)
