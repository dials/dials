Python_add_library(
    dials_algorithms_spot_finding_ext
    MODULE
    boost_python/spot_finding_ext.cc
)
target_link_libraries( dials_algorithms_spot_finding_ext PUBLIC CCTBX::cctbx Boost::python )