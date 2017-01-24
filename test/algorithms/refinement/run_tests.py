from __future__ import absolute_import, division
from libtbx import test_utils
import libtbx.load_env

param_tsts = (
    "$D/test/algorithms/refinement/tst_beam_parameters.py",
    "$D/test/algorithms/refinement/tst_crystal_parameters.py",
    "$D/test/algorithms/refinement/tst_detector_parameters.py",
    "$D/test/algorithms/refinement/tst_prediction_parameters.py",
    "$D/test/algorithms/refinement/tst_scan_varying_model_parameters.py",
    "$D/test/algorithms/refinement/tst_scan_varying_prediction_parameters.py",
    "$D/test/algorithms/refinement/tst_multi_panel_detector_parameterisation.py",
    "$D/test/algorithms/refinement/tst_finite_diffs.py",
    )

refinement_tsts = (
    "$D/test/command_line/tst_refine.py",
    "$D/test/algorithms/refinement/tst_orientation_refinement.py",
    "$D/test/algorithms/refinement/tst_refinement_regression.py",
    "$D/test/algorithms/refinement/tst_multi_experiment_refinement.py",
    "$D/test/algorithms/refinement/tst_stills_refinement.py",
    "$D/test/algorithms/refinement/tst_hierarchical_detector_refinement.py",
    )

misc_tsts = (
    "$D/test/algorithms/refinement/tst_dials_models.py",
    "$D/test/algorithms/refinement/tst_ref_passage_categorisation.py",
    )

def run (*args) :
  tst_list = reduce(lambda x, y: x+y, args)
  build_dir = libtbx.env.under_build("dials")
  dist_dir = libtbx.env.dist_path("dials")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if __name__ == "__main__":
  run(param_tsts, refinement_tsts, misc_tsts)
