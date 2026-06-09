Python_add_library( dials_util_ext MODULE boost_python/ext.cc )
target_link_libraries( dials_util_ext PUBLIC Boost::python CCTBX::cctbx CCTBX::ccp4io )

Python_add_library(
    dials_util_streambuf_test_ext
    MODULE
    boost_python/streambuf_test_ext.cpp
)
target_link_libraries( dials_util_streambuf_test_ext PUBLIC Boost::python CCTBX::cctbx )