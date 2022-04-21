Python_add_library(
    dials_algorithms_image_centroid_ext
    MODULE
    boost_python/centroid_ext.cc
)
target_link_libraries( dials_algorithms_image_centroid_ext PUBLIC CCTBX::cctbx Boost::python )