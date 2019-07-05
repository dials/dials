from __future__ import absolute_import, division, print_function

import os
import pytest
from dials.algorithms.refinement.outlier_detection import CentroidOutlierFactory
from dials.algorithms.refinement.outlier_detection.outlier_base import phil_scope
from dials.array_family import flex


@pytest.mark.parametrize(
    "method,colnames,expected_nout",
    [
        ("tukey", ("x_resid", "y_resid", "phi_resid"), 34),
        pytest.param("mcd", ("x_resid", "y_resid", "phi_resid"), 35),
        pytest.param("mcd", ("x_resid", "y_resid", "phi_resid"), 35),
        pytest.param(
            "sauter_poon", ("miller_index", "xyzobs.px.value", "xyzcal.px"), 34
        ),
    ],
)
def test_centroid_outlier(dials_regression, method, colnames, expected_nout):

    flex.set_random_seed(42)
    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "centroid_outlier"
    )
    residuals = flex.reflection_table.from_pickle(
        os.path.join(data_dir, "residuals.refl")
    )
    params = phil_scope.extract()
    params.outlier.algorithm = method
    params.outlier.sauter_poon.px_sz = (0.1, 0.1)  # must be set for SauterPoon
    outlier_detector = CentroidOutlierFactory.from_parameters_and_colnames(
        params, colnames
    )
    outlier_detector(residuals)
    outliers = residuals.get_flags(residuals.flags.centroid_outlier)

    assert outliers.count(True) == expected_nout
