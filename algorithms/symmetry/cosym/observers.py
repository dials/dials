"""
Observers for the cosym procedure.
"""
from __future__ import absolute_import, division, print_function
from collections import OrderedDict
import json
from dials.util.observer import Observer, singleton
from dials.algorithms.symmetry.cosym.plots import plot_coords, plot_rij_histogram
from dials.algorithms.symmetry.cosym import SymmetryAnalysis
from dials.algorithms.clustering.observers import UnitCellAnalysisObserver
from jinja2 import Environment, ChoiceLoader, PackageLoader


def register_default_cosym_observers(script):
    """Register the standard observers to the cosym script."""
    script.cosym_analysis.register_observer(
        event="analysed_symmetry", observer=SymmetryAnalysisObserver()
    )
    script.cosym_analysis.register_observer(
        event="analysed_clusters", observer=CosymClusterAnalysisObserver()
    )
    script.register_observer(event="run_cosym", observer=UnitCellAnalysisObserver())
    script.register_observer(
        event="run_cosym", observer=CosymHTMLGenerator(), callback="make_html"
    )
    script.register_observer(
        event="run_cosym", observer=CosymJSONGenerator(), callback="make_json"
    )


@singleton
class CosymHTMLGenerator(Observer):

    """
    Observer to make a html report
    """

    def make_html(self, cosym_script):
        """Collect data from the individual observers and write the html."""
        filename = cosym_script.params.output.html
        if not filename:
            return
        self.data.update(CosymClusterAnalysisObserver().make_plots())
        self.data.update(UnitCellAnalysisObserver().make_plots())
        self.data.update(SymmetryAnalysisObserver().make_tables())
        print("Writing html report to: %s" % filename)
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("cosym_report.html")
        html = template.render(
            page_title="DIALS cosym report",
            cosym_graphs=self.data["cosym_graphs"],
            unit_cell_graphs=self.data["unit_cell_graphs"],
            symmetry_analysis=self.data["symmetry_analysis"],
        )
        with open(filename, "wb") as f:
            f.write(html.encode("ascii", "xmlcharrefreplace"))


@singleton
class CosymJSONGenerator(Observer):

    """
    Observer to make a html report
    """

    def make_json(self, cosym_script):
        """Collect data from the individual observers and write the html."""
        filename = cosym_script.params.output.json
        if not filename:
            return
        self.data.update(CosymClusterAnalysisObserver().make_plots())
        self.data.update(UnitCellAnalysisObserver().make_plots())
        self.data.update(SymmetryAnalysisObserver().get_data())
        print("Writing json to: %s" % filename)
        with open(filename, "wb") as f:
            json.dump(self.data, f)


@singleton
class CosymClusterAnalysisObserver(Observer):

    """
    Observer to record cosym cluster analysis data and make model plots.
    """

    def update(self, cosym):
        """Update the data in the observer."""
        self.data["coordinates"] = cosym.coords
        self.data["labels"] = cosym.cluster_labels
        self.data["rij_matrix"] = cosym.target.rij_matrix

    def make_plots(self):
        """Generate cosym cluster analysis plot data."""
        d = OrderedDict()
        d.update(plot_rij_histogram(self.data["rij_matrix"]))
        d.update(plot_coords(self.data["coordinates"], self.data["labels"]))
        graphs = {"cosym_graphs": d}
        return graphs


@singleton
class SymmetryAnalysisObserver(Observer):

    """
    Observer to record symmetry analysis data and make tables.
    """

    def update(self, cosym):
        if cosym._symmetry_analysis is not None:
            self.data.update(cosym._symmetry_analysis.as_dict())

    def make_tables(self):
        """Generate symmetry analysis tables."""
        d = {"symmetry_analysis": {}}
        if self.data:
            d["symmetry_analysis"].update(
                {
                    "sym_ops_table": SymmetryAnalysis.sym_ops_table(self.data),
                    "subgroups_table": SymmetryAnalysis.subgroups_table(self.data),
                    "summary_table": SymmetryAnalysis.summary_table(self.data),
                }
            )
        return d

    def get_data(self):
        return self.data
