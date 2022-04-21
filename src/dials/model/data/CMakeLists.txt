Python_add_library(
    dials_model_data_ext
    MODULE
    boost_python/image_volume.cc
    boost_python/adjacency_list.cc
    boost_python/shoebox.cc
    boost_python/observation.cc
    boost_python/prediction.cc
    boost_python/pixel_list.cc
    boost_python/ray.cc
    boost_python/image.cc
    boost_python/data_ext.cc
)
target_link_libraries( dials_model_data_ext PUBLIC CCTBX::cctbx Boost::python )