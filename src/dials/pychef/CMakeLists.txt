Python_add_library( dials_pychef_ext MODULE ext.cc )
target_link_libraries( dials_pychef_ext PUBLIC CCTBX::cctbx Boost::python )