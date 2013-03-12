from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
    "$D/test/algorithms/spot_prediction/tst_index_generator.py",
    "$D/test/algorithms/spot_prediction/tst_ray_predictor.py",
    "$D/test/algorithms/spot_prediction/tst_rotation_angles.py",
    "$D/test/algorithms/spot_prediction/tst_spot_prediction.py",
    "$D/test/algorithms/integration/tst_from_beam_vector_to_xds.py",
    "$D/test/algorithms/integration/tst_from_xds_e3_to_phi.py",
    "$D/test/algorithms/integration/tst_from_xds_to_beam_vector.py",
    "$D/test/algorithms/integration/tst_from_xds_coordinate_system.py",
    "$B/test/algorithms/spatial_indexing/tst_quadtree",
    "$B/test/algorithms/spatial_indexing/tst_octree",
    )

def run () :
    build_dir = libtbx.env.under_build("dials")
    dist_dir = libtbx.env.dist_path("dials")
    test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
    run()
