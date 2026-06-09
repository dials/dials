Python_add_library(
    dials_algorithms_centroid_simple_ext
    MODULE
    simple/boost_python/simple_ext.cc
)
target_link_libraries( dials_algorithms_centroid_simple_ext PUBLIC CCTBX::cctbx Boost::python )