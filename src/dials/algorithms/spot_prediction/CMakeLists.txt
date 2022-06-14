Python_add_library(
    dials_algorithms_spot_prediction_ext
    MODULE
    boost_python/index_generator.cc
    boost_python/reeke_index_generator.cc
    boost_python/scan_varying_helpers.cc
    boost_python/rotation_angles.cc
    boost_python/ray_predictor.cc
    boost_python/scan_varying_ray_predictor.cc
    boost_python/stills_ray_predictor.cc
    boost_python/ray_intersection.cc
    boost_python/reflection_predictor.cc
    boost_python/pixel_labeller.cc
    boost_python/pixel_to_miller_index.cc
    boost_python/spot_prediction_ext.cc
)
target_link_libraries( dials_algorithms_spot_prediction_ext PUBLIC CCTBX::cctbx Boost::python )