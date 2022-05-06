Python_add_library(
    dials_algorithms_image_filter_ext
    MODULE
    boost_python/summed_area.cc
    boost_python/mean_and_variance.cc
    boost_python/index_of_dispersion_filter.cc
    boost_python/convolve.cc
    boost_python/median.cc
    boost_python/distance.cc
    boost_python/anisotropic_diffusion.cc
    boost_python/filter_ext.cc
)
target_link_libraries( dials_algorithms_image_filter_ext PUBLIC CCTBX::cctbx Boost::python )