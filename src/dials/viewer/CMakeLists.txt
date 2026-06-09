Python_add_library(
    dials_viewer_ext
    MODULE
    boost_python/RGB_2d.cc
    boost_python/RGB_2d_ext.cc
)
target_link_libraries( dials_viewer_ext PUBLIC CCTBX::cctbx Boost::python )