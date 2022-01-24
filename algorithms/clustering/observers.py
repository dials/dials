from __future__ import annotations

from scitbx.array_family import flex

from dials.algorithms.clustering import plots
from dials.util.observer import Observer, singleton


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
        try:
            self.data["dendrogram"] = script.unit_cell_dendrogram
        except AttributeError:
            pass
        self.data["experiments"] = script._experiments

    def make_plots(self):
        """Generate plots of the unit cell clustering."""

        uc_params = uc_params_from_experiments(self.data["experiments"])
        d = plots.plot_uc_histograms(uc_params)

        if "dendrogram" in self.data:
            d["uc_clustering"] = plots.scipy_dendrogram_to_plotly_json(
                self.data["dendrogram"],
                title="Unit cell clustering",
                xtitle="Dataset",
                ytitle=r"Distance (Å<sup>2</sup>)",
                help="""\
The results of single-linkage hierarchical clustering on the unit cell parameters using
the Andrews-Bernstein NCDist distance metric (Andrews & Bernstein, 2014). The height at
which two clusters are merged in the dendrogram is a measure of the similarity between
the unit cells in each cluster. A larger separation between two clusters may be
indicative of a higher degree of non-isomorphism between the clusters. Conversely, a
small separation between two clusters suggests that their unit cell parameters are
relatively isomorphous.
""",
            )

        graphs = {"unit_cell_graphs": d}

        return graphs
