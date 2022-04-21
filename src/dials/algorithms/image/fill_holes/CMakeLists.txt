Python_add_library(
    dials_algorithms_image_fill_holes_ext
    MODULE
    boost_python/fill_holes_ext.cc
)
target_link_libraries( dials_algorithms_image_fill_holes_ext PUBLIC CCTBX::cctbx Boost::python )