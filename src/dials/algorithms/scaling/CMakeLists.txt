Python_add_library(
    dials_scaling_ext
    MODULE
    boost_python/scaling_helper.cc
    boost_python/scaling_ext.cc
)
target_link_libraries( dials_scaling_ext PUBLIC CCTBX::cctbx Boost::python )