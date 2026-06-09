import os
import logging

import hist as hi
import uproot

import pepper

logger = logging.getLogger(__name__)


class Processor(pepper.Processor):
    def __init__(self, config, eventdir):
        self.data_pu_hist = config["data_pu_hist"]
        self.data_pu_hist_up = config["data_pu_hist_up"]
        self.data_pu_hist_down = config["data_pu_hist_down"]
        self.year = config["year"]
        config.check_integrity = False

        datahist, datahistup, datahistdown = self.load_input_hists()
        if len(datahist.axes) != 1:
            raise pepper.config.ConfigError(
                "data_pu_hist has invalid number of axes. Only one axis is "
                "allowed.")
        if datahist.axes != datahistup.axes:
            raise pepper.config.ConfigError(
                "data_pu_hist_up does not have the same axes as "
                "data_pu_hist.")
        if datahist.axes != datahistdown.axes:
            raise pepper.config.ConfigError(
                "data_pu_hist_down does not have the same axes as "
                "data_pu_hist.")

        self.axisname = datahist.axes[0].name
        hist_config = {
            "bins": [
                {
                    "name": self.axisname,
                    "label": ("True mean number interactions per bunch "
                              "crossing")
                }
            ],
            "fill": {
                self.axisname: [
                    "Pileup",
                    "nTrueInt",
                    {"function": "int"}
                ]
            }
        }
        if (isinstance(datahist.axes[0], hi.axis.Regular)
                or isinstance(datahist.axes[0], hi.axis.Integer)):
            ax = datahist.axes[0]
            lo_edge = ax.value(0)
            hi_edge = ax.value(len(ax))
            hist_config["bins"][0].update({
                "n_or_arr": len(ax),
                "lo": lo_edge,
                "hi": hi_edge
            })
            if (int(lo_edge) == lo_edge) and (int(hi_edge) == hi_edge):
                # Integer bin edges may cause numerical problems for regular
                # axes - use integer axis if possible
                if (hi_edge - lo_edge) == len(ax):
                    hist_config["bins"][0].update({
                        "lo": int(lo_edge),
                        "hi": int(hi_edge),
                        "type": "int"
                    })
                else:
                    logger.warning(
                        "Using regular axis, which may have issues with "
                        "integer bin edges - please check your output")
        else:
            logger.warning(
                "Using variable axis, which may have issues with integer bin "
                "edges - please check your output")
            hist_config["bins"][0]["n_or_arr"] = datahist.axes[0].edges
        config["hists"] = {"pileup": pepper.HistDefinition(hist_config)}
        if "hists_to_do" in config:
            del config["hists_to_do"]
        if "cuts_to_histogram" in config:
            del config["cuts_to_histogram"]
        config["compute_systematics"] = False
        # Treat all datasets as normal datasets, instead of using them as
        # systematics
        config["dataset_for_systematics"] = {}
        # Set mc_lumifactors to true, to prevent pepper trying to calculate these for the case
        # of mc_lumifactors: false. These are irrelevant in any case when computing pileup weights
        config["mc_lumifactors"] = True

        super().__init__(config, eventdir)

    def preprocess(self, datasets):
        # Only run over MC
        processed = {}
        for key, value in datasets.items():
            if key in self.config["mc_datasets"]:
                processed[key] = value
        return processed

    def setup_selection(self, data, dsname, is_mc, filler):
        # Ignore generator weights, because pileup is independent
        return pepper.Selector(data, on_update=filler.get_callbacks())

    def process_selection(self, selector, dsname, is_mc, filler):
        pass

    def load_input_hists(self):
        with uproot.open(self.data_pu_hist) as f:
            datahist = f["pileup"].to_hist()
        with uproot.open(self.data_pu_hist_up) as f:
            datahistup = f["pileup"].to_hist()
        with uproot.open(self.data_pu_hist_down) as f:
            datahistdown = f["pileup"].to_hist()
        return datahist, datahistup, datahistdown

    @staticmethod
    def _save_hists(hist, datahist, datahistup, datahistdown, filename,
                    sum_datasets):
        datasets = ["all_datasets"] if sum_datasets else hist.axes["dataset"]

        with uproot.recreate(filename) as f:
            for dataset in datasets:
                axpos = sum if sum_datasets else dataset
                denom = hist[{"dataset": axpos}].values().copy()
                # Set bins that are zero in MC to 0 in data to get norm right
                is_nonzero = denom != 0
                # Avoid division by zero warning
                denom[~is_nonzero] = 1
                denom /= denom.sum()
                for datahist_i, suffix in [
                        (datahist, ""), (datahistup, "_up"),
                        (datahistdown, "_down")]:
                    datahist_i = datahist_i * is_nonzero
                    if isinstance(datahist_i.storage_type(),
                                  hi.storage.Weight):
                        norm = datahist_i.sum().value
                    else:
                        norm = datahist_i.sum()
                    ratio = datahist_i / norm / denom

                    f[dataset + suffix] = ratio

    def save_output(self, output, dest):
        datahist, datahistup, datahistdown = self.load_input_hists()

        mchist = None
        for dataset, hists in output["hists"].items():
            if mchist is None:
                mchist = hists[("BeforeCuts", "pileup")].copy()
            else:
                mchist += hists[("BeforeCuts", "pileup")]
        # Set underflow and 0 pileup bin to 0, which might be != 0 only for
        # buggy reasons in MC
        axidx = [ax.name for ax in mchist.axes].index(self.axisname)
        slic = (slice(None),) * axidx + (slice(None, 2),)
        mchist.view(flow=True)[slic].fill(0)

        self._save_hists(
            mchist, datahist, datahistup, datahistdown,
            os.path.join(dest, "pileup.root"), True)
        self._save_hists(
            mchist, datahist, datahistup, datahistdown,
            os.path.join(dest, "pileup_perdataset.root"), False)


if __name__ == "__main__":
    from pepper import runproc
    runproc.run_processor(
        Processor, "Create histograms needed for pileup reweighting",
        mconly=True)
