Python_add_library(
    dials_algorithms_symmetry_refstat_ext
    MODULE
    boost_python/refstat.cpp
)
target_link_libraries( dials_algorithms_symmetry_refstat_ext PUBLIC CCTBX::cctbx Boost::python )
