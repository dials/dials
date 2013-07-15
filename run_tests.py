from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
    "$D/test/tst_spots_xds.py",
    "$D/test/model/data/tst_reflection_pickle.py",
    "$D/test/algorithms/spot_prediction/tst_index_generator.py",
    "$D/test/algorithms/spot_prediction/tst_ray_predictor.py",
    "$D/test/algorithms/spot_prediction/tst_rotation_angles.py",
    "$D/test/algorithms/spot_prediction/tst_spot_prediction.py",
    "$D/test/algorithms/integration/tst_from_beam_vector_to_xds.py",
    "$D/test/algorithms/integration/tst_from_xds_e3_to_phi.py",
    "$D/test/algorithms/integration/tst_from_xds_to_beam_vector.py",
    "$D/test/algorithms/integration/tst_xds_coordinate_system.py",
    "$D/test/algorithms/integration/tst_bbox_calculator.py",
    "$D/test/algorithms/integration/profile/tst_grid_sampler.py",
    "$D/test/algorithms/integration/profile/tst_xds_circle_sampler.py",
    "$D/test/algorithms/integration/profile/tst_reference_locator.py",
    "$D/test/algorithms/integration/profile/tst_reference_learner.py",
    "$D/test/algorithms/integration/profile/tst_profile_fitting.py",
    "$D/test/algorithms/centroid/tst_filtered_centroid.py",
    "$D/test/algorithms/centroid/tst_lui_centroid.py",
    "$D/test/algorithms/centroid/tst_toy_centroid.py",
    "$D/test/algorithms/centroid/tst_toy_centroid_helpers.py",
    "$B/test/algorithms/spatial_indexing/tst_quadtree",
    "$B/test/algorithms/spatial_indexing/tst_octree",
    "$B/test/algorithms/spatial_indexing/tst_collision_detection",
    "$D/test/algorithms/refinement/tst_refinement_regression.py",
    "$D/test/algorithms/image/tst_centroid.py",
    "$D/test/algorithms/polygon/clip/tst_clipping.py",
    "$D/test/algorithms/reflexion_basis/tst_coordinate_system.py",
    "$D/test/algorithms/reflexion_basis/tst_map_frames.py",
    "$D/test/algorithms/reflexion_basis/tst_rebin_pixels.py",
    "$D/test/algorithms/reflexion_basis/tst_beam_vector_map.py",
    "$D/test/algorithms/peak_finding/tst_spotfinder.py",
    )

def run () :
    build_dir = libtbx.env.under_build("dials")
    dist_dir = libtbx.env.dist_path("dials")
    test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
    run()
