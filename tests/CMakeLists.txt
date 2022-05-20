include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)


add_executable( tst_reeke_model algorithms/spot_prediction/tst_reeke_model.cc )
set_target_properties( tst_reeke_model PROPERTIES RUNTIME_OUTPUT_DIRECTORY "algorithms/spot_prediction")
target_link_libraries( tst_reeke_model PUBLIC Boost::boost CCTBX::scitbx )

add_executable( tst_collision_detection algorithms/spatial_indexing/tst_collision_detection.cc )
set_target_properties( tst_collision_detection PROPERTIES RUNTIME_OUTPUT_DIRECTORY "algorithms/spatial_indexing")
target_link_libraries( tst_collision_detection PUBLIC Boost::boost CCTBX::scitbx )
