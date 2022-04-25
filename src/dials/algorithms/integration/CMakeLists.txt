Python_add_library(
    dials_algorithms_integration_ext
    MODULE
    boost_python/corrections.cc
    boost_python/integration_ext.cc
)
target_link_libraries( dials_algorithms_integration_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_integration_integrator_ext
    MODULE
    boost_python/integration_integrator_ext.cc
)
target_link_libraries( dials_algorithms_integration_integrator_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_integration_parallel_integrator_ext
    MODULE
    boost_python/parallel_integrator_ext.cc
)
target_link_libraries( dials_algorithms_integration_parallel_integrator_ext PUBLIC Boost::thread Boost::python CCTBX::cctbx )

Python_add_library(
    dials_algorithms_integration_sum_ext
    MODULE
    sum/boost_python/summation.cc
    sum/boost_python/integration_sum_ext.cc
)
target_link_libraries( dials_algorithms_integration_sum_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library( dials_algorithms_integration_fit_ext MODULE fit/boost_python/ext.cc )
target_link_libraries( dials_algorithms_integration_fit_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library( dials_algorithms_integration_bayes_ext MODULE bayes/boost_python/ext.cc )
target_link_libraries( dials_algorithms_integration_bayes_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library( dials_algorithms_integration_kapton_ext MODULE boost_python/kapton_ext.cc )
target_link_libraries( dials_algorithms_integration_kapton_ext PUBLIC CCTBX::cctbx Boost::python )