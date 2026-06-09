Python_add_library(
    dials_algorithms_simulation_ext
    MODULE
    boost_python/simulation_ext.cc
)
target_link_libraries( dials_algorithms_simulation_ext PUBLIC CCTBX::cctbx Boost::python )