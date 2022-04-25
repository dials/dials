Python_add_library( dials_model_serialize_ext MODULE boost_python/serialize_ext.cc )
target_link_libraries( dials_model_serialize_ext CCTBX::cctbx Boost::python )