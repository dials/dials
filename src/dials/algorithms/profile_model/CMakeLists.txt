Python_add_library(
    dials_algorithms_profile_model_modeller_ext
    MODULE
    modeller/boost_python/sampler.cc
    modeller/boost_python/modeller.cc
    modeller/boost_python/ext.cc
)
target_link_libraries( dials_algorithms_profile_model_modeller_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_profile_model_gaussian_rs_ext
    MODULE
    gaussian_rs/boost_python/gaussian_rs_ext.cc
)
target_link_libraries( dials_algorithms_profile_model_gaussian_rs_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_profile_model_gaussian_rs_transform_ext
    MODULE
    gaussian_rs/transform/boost_python/transform_ext.cc
)
target_link_libraries( dials_algorithms_profile_model_gaussian_rs_transform_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_profile_model_ellipsoid_ext
    ellipsoid/boost_python/ext.cc
)
target_link_libraries(
    dials_algorithms_profile_model_ellipsoid_ext
    PUBLIC
    CCTBX::cctbx Boost::python
)
