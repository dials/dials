import mock

from scitbx.array_family import flex
from dials.algorithms.symmetry.cosym import observers


def test_SymmetryAnalysisObserver():
    # setup script
    script = mock.Mock()

    # test when no symmetry analysis has been performed
    script._symmetry_analysis = None
    observer = observers.SymmetryAnalysisObserver()
    observer.update(script)
    d = observer.make_tables()
    assert d == {"symmetry_analysis": {}}

    script._symmetry_analysis = mock.Mock()
    script._symmetry_analysis.sym_ops_table = mock.Mock()
    script._symmetry_analysis.subgroups_table = mock.Mock()
    script._symmetry_analysis.as_dict = mock.Mock(
        return_value={
            "subgroup_scores": [
                {
                    "patterson_group": "-P 1",
                    "unit_cell": (10, 10, 10, 90, 90, 90),
                    "cb_op": "x,y,z",
                    "likelihood": 0.9,
                    "confidence": 0.9,
                }
            ]
        }
    )

    # test the observer
    observer = observers.SymmetryAnalysisObserver()
    observer.update(script)
    d = observer.make_tables()
    assert "symmetry_analysis" in d
    assert d["symmetry_analysis"].keys() == [
        "summary_table",
        "subgroups_table",
        "sym_ops_table",
    ]


def test_CosymClusterAnalysisObserver():
    rij_matrix = flex.random_double(16)
    rij_matrix.reshape(flex.grid(4, 4))
    coords = flex.random_double(8)
    coords.reshape(flex.grid(4, 2))

    # setup script
    script = mock.Mock()
    script.target = mock.Mock()
    script.target.rij_matrix = rij_matrix
    script.coords = coords
    script.cluster_labels = flex.int(4, 0)

    # test the observer
    observer = observers.CosymClusterAnalysisObserver()
    observer.update(script)
    d = observer.make_plots()
    assert "cosym_graphs" in d


def test_CosymHTMLGenerator():
    pass
    # script = mock.Mock()
    # script.params.output.html = "test.html"

    ## Test that CosymHTMLGenerator works if all data is empty.
    # observer = observers.CosymHTMLGenerator()
    # observer.make_html(script)
    # assert os.path.exists("test.html")
