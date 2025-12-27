Python_add_library(
    dials_cosym_ext
    MODULE
    boost_python/cosym_helper.cc
    boost_python/cosym_ext.cc
)
target_link_libraries( dials_cosym_ext PUBLIC CCTBX::cctbx Boost::python )