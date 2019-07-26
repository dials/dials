"""
Observers for the cosym procedure.
"""
from __future__ import absolute_import, division, print_function
from collections import OrderedDict
from dials.util.observer import Observer, singleton
from dials.algorithms.symmetry.cosym.plots import plot_coords, plot_rij_histogram
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


@singleton
class CosymHTMLGenerator(Observer):

    """
    Observer to make a html report
    """

    def make_html(self, cosym_script):
        """Collect data from the individual observers and write the html."""
        self.data.update(CosymClusterAnalysisObserver().make_plots())
        self.data.update(UnitCellAnalysisObserver().make_plots())
        self.data.update(SymmetryAnalysisObserver().make_tables())
        filename = cosym_script.params.output.html
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
        if cosym._symmetry_analysis is None:
            return
        d = cosym._symmetry_analysis.as_dict()
        self.data["sym_ops_table"] = cosym._symmetry_analysis.sym_ops_table(d)
        self.data["subgroups_table"] = cosym._symmetry_analysis.subgroups_table(d)
        self.data["summary_table"] = cosym._symmetry_analysis.summary_table(d)

    def make_tables(self):
        """Generate symmetry analysis tables."""
        d = {"symmetry_analysis": self.data}

        return d
