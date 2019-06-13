from __future__ import absolute_import, division, print_function

import pytest
from dials.algorithms.refinement.outlier_detection import CentroidOutlierFactory
from dials.algorithms.refinement.outlier_detection.outlier_base import phil_scope
from dials.array_family import flex
from dials.test.algorithms.refinement.test_stills_prediction_parameters import _Test
from scitbx.python_utils import random_transform


@pytest.fixture(scope="session")
def residuals():
    # Get a reflection table with simulated reflections
    r = _Test().reflections
    n = len(r)

    # Add some columns of random residuals
    flex.set_random_seed(42)
    r = flex.reflection_table.empty_standard(n)
    r["x_resid"] = random_transform.normal_variate(mu=0.0, sigma=1.0, N=n)
    r["y_resid"] = random_transform.normal_variate(mu=0.0, sigma=1.0, N=n)
    r["phi_resid"] = random_transform.normal_variate(mu=0.0, sigma=1.0, N=n)

    # Ensure the calculated positions match the residuals
    x, y, z = r["xyzobs.px.value"].parts()
    x += r["x_resid"]
    y += r["y_resid"]
    z += r["phi_resid"]
    r["xyzcal.px"] = flex.vec3_double(x, y, z)

    # Set all reflections as 'predicted'
    r.set_flags(flex.bool(n, True), r.flags.predicted)

    return r


@pytest.mark.parametrize(
    "method,colnames,expected_nout",
    [
        ("tukey", ("x_resid", "y_resid", "phi_resid"), 21),
        ("mcd", ("x_resid", "y_resid", "phi_resid"), 34),
        ("sauter_poon", ("miller_index", "xyzobs.px.value", "xyzcal.px"), 34),
    ],
)
def test_centroid_outlier(residuals, method, colnames, expected_nout):

    params = phil_scope.extract()
    params.outlier.algorithm = method
    params.outlier.sauter_poon.px_sz = (0.1, 0.1)  # must be set for SauterPoon
    outlier_detector = CentroidOutlierFactory.from_parameters_and_colnames(
        params, colnames
    )
    outlier_detector(residuals)
    outliers = residuals.get_flags(residuals.flags.centroid_outlier)

    assert outliers.count(True) == expected_nout
