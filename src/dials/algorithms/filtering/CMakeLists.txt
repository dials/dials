Python_add_library( dials_algorithms_filter_ext MODULE boost_python/filter_ext.cc )
target_link_libraries( dials_algorithms_filter_ext PUBLIC CCTBX::cctbx Boost::python )