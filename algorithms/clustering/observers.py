# -*- coding: utf-8 -*-

from collections import OrderedDict

from scitbx.array_family import flex
from dials.util.observer import Observer, singleton
from dials.algorithms.clustering import plots


def uc_params_from_experiments(experiments):
    uc_params = [flex.double() for i in range(6)]
    for expt in experiments:
        uc = expt.crystal.get_unit_cell()
        for i in range(6):
            uc_params[i].append(uc.parameters()[i])
    return uc_params


@singleton
class UnitCellAnalysisObserver(Observer):

    """
    Observer to record unit cell clustering data and make plots.
    """

    def update(self, script):
        """Update the data in the observer."""
        self.data["dendrogram"] = script.unit_cell_dendrogram
        self.data["experiments"] = script._experiments

    def make_plots(self):
        """Generate plots of the unit cell clustering."""

        # experiments = self._data_manager.experiments
        uc_params = uc_params_from_experiments(self.data["experiments"])

        d = OrderedDict()
        d.update(plots.plot_uc_histograms(uc_params))
        d["uc_clustering"] = plots.scipy_dendrogram_to_plotly_json(
            self.data["dendrogram"],
            title="Unit cell clustering",
            xtitle="Dataset",
            ytitle=r"Distance (Å^2)",
        )

        graphs = {"unit_cell_graphs": d}

        return graphs
