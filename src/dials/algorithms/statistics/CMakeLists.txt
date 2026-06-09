Python_add_library(
    dials_algorithms_statistics_ext
    MODULE
    boost_python/statistics_ext.cc
)
target_link_libraries( dials_algorithms_statistics_ext PUBLIC CCTBX::cctbx Boost::python )