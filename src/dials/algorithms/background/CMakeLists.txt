Python_add_library(
    dials_algorithms_background_ext
    MODULE
    boost_python/helpers.cc
    boost_python/background_ext.cc
)
target_link_libraries( dials_algorithms_background_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library( dials_algorithms_background_modeller_ext MODULE boost_python/modeller_ext.cc )
target_link_libraries( dials_algorithms_background_modeller_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_background_simple_ext
    MODULE
    simple/boost_python/outlier_rejector.cc
    simple/boost_python/modeller.cc
    simple/boost_python/creator.cc
    simple/boost_python/background_simple_ext.cc
)
target_link_libraries( dials_algorithms_background_simple_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library(
    dials_algorithms_background_glm_ext
    MODULE
    glm/boost_python/ext.cc
)
target_link_libraries( dials_algorithms_background_glm_ext PUBLIC CCTBX::cctbx Boost::python )

Python_add_library( dials_algorithms_background_gmodel_ext MODULE gmodel/boost_python/ext.cc )
target_link_libraries( dials_algorithms_background_gmodel_ext PUBLIC CCTBX::cctbx Boost::python )