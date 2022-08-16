Python_add_library(
    dials_algorithms_shoebox_ext
    MODULE
    boost_python/mask_code.cc
    boost_python/find_overlapping.cc
    boost_python/overload_checker.cc
    boost_python/mask_empirical.cc
    boost_python/mask_overlapping.cc
    boost_python/mask_builder.cc
    boost_python/shoebox_ext.cc
)
target_link_libraries( dials_algorithms_shoebox_ext PUBLIC Boost::python CCTBX::cctbx CCTBX::annlib )