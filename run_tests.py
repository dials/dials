from __future__ import absolute_import, division, print_function

from libtbx.test_utils.pytest import discover

tst_list = [
    "$B/test/algorithms/spatial_indexing/tst_collision_detection",
    "$B/test/algorithms/spatial_indexing/tst_octree",
    "$B/test/algorithms/spatial_indexing/tst_quadtree",
    "$B/test/algorithms/spot_prediction/tst_reeke_model",
    ] + discover()
