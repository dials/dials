Python_add_library( dials_algorithms_polygon_clip_ext MODULE boost_python/clip_ext.cc )
target_link_libraries( dials_algorithms_polygon_clip_ext PUBLIC CCTBX::cctbx Boost::python )