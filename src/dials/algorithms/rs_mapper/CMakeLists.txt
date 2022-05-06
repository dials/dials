Python_add_library( recviewer_ext MODULE ext.cpp )
target_link_libraries( recviewer_ext PUBLIC CCTBX::cctbx Boost::python )