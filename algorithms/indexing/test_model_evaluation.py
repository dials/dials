from __future__ import absolute_import, division, print_function

import copy
import os
import pytest

from cctbx import sgtbx, uctbx
from dxtbx.model import Crystal
from dxtbx.serialize import load
from dials.array_family import flex
from dials.algorithms.indexing import model_evaluation
from dials.algorithms.indexing.assign_indices import AssignIndicesGlobal
from dials.algorithms.refinement.refiner import phil_scope as refine_phil


def test_Result():

    rmsds = (0.0269546226916153, 0.03159452239902618, 0.004208711548775884)
    fraction_indexed = 0.9683976336148051
    model_likelihood = 0.95846970346759
    n_indexed = 19152
    crystal = Crystal.from_dict(
        dict(
            [
                ("__id__", "crystal"),
                (
                    "real_space_a",
                    (-45.031874961773504, -14.833784919966813, -48.766092343826806),
                ),
                (
                    "real_space_b",
                    (-32.15263553253188, 3.8725711269478085, 59.82290857456796),
                ),
                (
                    "real_space_c",
                    (29.6683003477066, 59.732569113301565, -13.880892871529552),
                ),
                ("space_group_hall_symbol", " P 1"),
            ]
        )
    )

    result = model_evaluation.Result(
        rmsds=rmsds,
        fraction_indexed=fraction_indexed,
        model_likelihood=model_likelihood,
        n_indexed=n_indexed,
        crystal=crystal,
    )

    assert result.rmsds == rmsds
    assert result.fraction_indexed == fraction_indexed
    assert result.model_likelihood == model_likelihood
    assert result.n_indexed == n_indexed
    assert result.crystal == crystal


def test_ModelRank():

    results = [
        model_evaluation.Result(
            rmsds=(0.0269546226916153, 0.03159452239902618, 0.004208711548775884),
            fraction_indexed=0.9683976336148051,
            model_likelihood=0.95846970346759,
            n_indexed=19152,
            crystal=Crystal.from_dict(
                dict(
                    [
                        ("__id__", "crystal"),
                        (
                            "real_space_a",
                            (
                                -45.031874961773504,
                                -14.833784919966813,
                                -48.766092343826806,
                            ),
                        ),
                        (
                            "real_space_b",
                            (-32.15263553253188, 3.8725711269478085, 59.82290857456796),
                        ),
                        (
                            "real_space_c",
                            (29.6683003477066, 59.732569113301565, -13.880892871529552),
                        ),
                        ("space_group_hall_symbol", " P 1"),
                    ]
                )
            ),
        ),
        model_evaluation.Result(
            rmsds=(0.0341397658828684, 0.027401396596305812, 0.00427723439147068),
            fraction_indexed=0.9849825554937554,
            model_likelihood=0.9562237490188447,
            n_indexed=19480,
            crystal=Crystal.from_dict(
                dict(
                    [
                        ("__id__", "crystal"),
                        (
                            "real_space_a",
                            (29.66830034770662, 59.73256911330157, -13.880892871529573),
                        ),
                        (
                            "real_space_b",
                            (47.516210146598816, -48.77135532028254, 2.824076640788394),
                        ),
                        (
                            "real_space_c",
                            (
                                -45.036407933560845,
                                -14.950807536025826,
                                -49.06808637024198,
                            ),
                        ),
                        ("space_group_hall_symbol", " P 1"),
                    ]
                )
            ),
        ),
        model_evaluation.Result(
            rmsds=(0.30456791867888355, 0.15679214175133024, 0.009635577811258947),
            fraction_indexed=0.33629974212469027,
            model_likelihood=0.6574428619874397,
            n_indexed=6651,
            crystal=Crystal.from_dict(
                dict(
                    [
                        ("__id__", "crystal"),
                        (
                            "real_space_a",
                            (
                                -11.907050303571122,
                                34.85499418820148,
                                30.689745759790572,
                            ),
                        ),
                        (
                            "real_space_b",
                            (-56.943458237132, 19.90418665217566, -18.37834061045143),
                        ),
                        (
                            "real_space_c",
                            (-41.63267941211685, -24.82747139443437, 44.5508337593274),
                        ),
                        ("space_group_hall_symbol", " P 1"),
                    ]
                )
            ),
        ),
    ]

    ranker = model_evaluation.ModelRankWeighted()
    ranker.extend(results)
    assert list(ranker.score_by_fraction_indexed()) == pytest.approx(
        [0.02449862002620243, 0.0, 1.550350501381241]
    )
    assert list(ranker.score_by_rmsd_xy()) == pytest.approx(
        [0.0, 0.07598423386955666, 3.044108569529155]
    )
    assert list(ranker.score_by_volume()) == pytest.approx(
        [0.44207271296753703, 0.45026641813391066, 0.0]
    )
    assert list(ranker.combined_scores()) == pytest.approx(
        [0.19602846593366663, 0.20851345109588518, 11.670183660213905]
    )
    assert (
        str(ranker)
        == """\
----------------------------------------------------------------------------------------------------------------------------------------------
unit_cell                           | volume | volume score | #indexed | % indexed | % indexed score | rmsd_xy | rmsd_xy score | overall score
----------------------------------------------------------------------------------------------------------------------------------------------
68.02 68.03 68.12 109.6 109.5 109.3 | 242890 | 0.44         | 19152    | 97        | 0.02            | 0.04    | 0.00          | 0.20
68.12 68.15 68.26 109.5 109.4 109.4 | 244274 | 0.45         | 19480    | 98        | 0.00            | 0.04    | 0.08          | 0.21
47.94 63.06 65.84 75.2 71.6 74.5    | 178786 | 0.00         | 6651     | 34        | 1.55            | 0.34    | 3.04          | 11.67
----------------------------------------------------------------------------------------------------------------------------------------------"""
    )
    best = ranker.best_model()
    assert best.n_indexed == 19152

    ranker = model_evaluation.ModelRankFilter()
    ranker.extend(results)
    best = ranker.best_model()
    assert best.n_indexed == 19152
    assert (
        str(ranker)
        == """\
----------------------------------------------------------------------------------------
unit_cell                           | volume | n_indexed | fraction_indexed | likelihood
----------------------------------------------------------------------------------------
68.02 68.03 68.12 109.6 109.5 109.3 | 242890 | 19152     | 97               | 0.96
68.12 68.15 68.26 109.5 109.4 109.4 | 244274 | 19480     | 98               | 0.96
47.94 63.06 65.84 75.2 71.6 74.5    | 178786 | 6651      | 34               | 0.66
----------------------------------------------------------------------------------------"""
    )


def test_ModelEvaluation(dials_regression, tmpdir):
    # thaumatin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sweep_path = os.path.join(data_dir, "experiments_import.json")

    input_reflections = flex.reflection_table.from_pickle(pickle_path)
    input_experiments = load.experiment_list(sweep_path, check_format=False)

    input_reflections = input_reflections.select(
        input_reflections["xyzobs.px.value"].parts()[2] < 100
    )
    input_reflections.centroid_px_to_mm(
        input_experiments[0].detector, scan=input_experiments[0].scan
    )
    input_reflections.map_centroids_to_reciprocal_space(
        input_experiments[0].detector,
        input_experiments[0].beam,
        goniometer=input_experiments[0].goniometer,
    )
    input_reflections["imageset_id"] = flex.size_t(input_reflections.size(), 0)
    input_reflections["id"] = flex.int(input_reflections.size(), -1)

    refine_params = refine_phil.fetch().extract()
    evaluator = model_evaluation.ModelEvaluation(refinement_params=refine_params)

    experiments = copy.deepcopy(input_experiments)
    reflections = copy.deepcopy(input_reflections)
    experiments[0].crystal = Crystal.from_dict(
        dict(
            [
                ("__id__", "crystal"),
                (
                    "real_space_a",
                    (20.007058080503633, 49.721143642677994, 16.636052132572846),
                ),
                (
                    "real_space_b",
                    (-15.182202482876685, 24.93846318493148, -50.7116866438356),
                ),
                (
                    "real_space_c",
                    (-135.23051191036296, 41.14539066294313, 55.41374425160883),
                ),
                ("space_group_hall_symbol", " P 1"),
            ]
        )
    )

    assign_indices = AssignIndicesGlobal()
    assign_indices(reflections, experiments)
    result = evaluator.evaluate(experiments, reflections)
    assert result is not None
    assert result.n_indexed == 7313
    assert result.fraction_indexed == pytest.approx(0.341155066244)
    assert result.rmsds == pytest.approx(
        (0.10215787570953785, 0.12954412140128563, 0.0010980583382102509)
    )


def test_filter_doubled_cell():
    import scitbx.matrix

    sgi = sgtbx.space_group_info("P1")
    uc1 = sgi.any_compatible_unit_cell(volume=1000)
    params = uc1.parameters()
    for (m1, m2, m3) in (
        (2, 1, 1),
        (1, 2, 1),
        (1, 1, 2),
        (2, 2, 1),
        (2, 1, 2),
        (1, 2, 2),
        (2, 2, 2),
    ):
        uc2 = uctbx.unit_cell(
            (
                params[0] * m1,
                params[1] * m2,
                params[2] * m3,
                params[3],
                params[4],
                params[5],
            )
        )
        B1 = scitbx.matrix.sqr(uc1.fractionalization_matrix()).transpose()
        B2 = scitbx.matrix.sqr(uc2.fractionalization_matrix()).transpose()
        solutions = [
            model_evaluation.Result(crystal=Crystal(B1, sgi.group()), n_indexed=100),
            model_evaluation.Result(crystal=Crystal(B2, sgi.group()), n_indexed=100),
        ]
        filtered = model_evaluation.filter_doubled_cell(solutions)
        assert len(filtered) == 1
        assert filtered[0].crystal == solutions[0].crystal
