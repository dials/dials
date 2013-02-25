/*
 * detector_helpers.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/format.hpp>
#include <string>
#include <dials/model/experiment/detector.h>
#include <dials/model/experiment/detector_helpers.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  template <typename DetectorType>
  is_coordinate_valid <DetectorType> make_is_coordinate_valid(
      const DetectorType& detector) {
    return is_coordinate_valid<DetectorType>(detector);
  }

  template <typename DetectorType>
  diffracted_beam_intersection <DetectorType> make_diffracted_beam_intersection(
      const DetectorType& detector) {
    return diffracted_beam_intersection<DetectorType>(detector);   
  }

  void export_is_coordinate_valid()
  {
    class_<is_coordinate_valid<FlatPanelDetector> >(
        "FlatPanelDetector_is_coordinate_valid", no_init)
    .def("__call__", 
      &is_coordinate_valid<FlatPanelDetector>::operator() <vec2 <double> >);
    
    class_<is_coordinate_valid<MultiFlatPanelDetector> >(
        "MultiFlatPanelDetector_is_coordinate_valid", no_init)
    .def("__call__", 
      &is_coordinate_valid<MultiFlatPanelDetector>::operator() <vec3 <double> >);    
    
    def("is_coordinate_valid", 
      &make_is_coordinate_valid<FlatPanelDetector>);
    def("is_coordinate_valid", 
      &make_is_coordinate_valid<MultiFlatPanelDetector>);
  }

  void export_diffracted_beam_intersection()
  {
    class_<diffracted_beam_intersection<FlatPanelDetector> >(
        "FlatPanelDetector_diffracted_beam_intersection", no_init)
    .def("__call__", 
      &diffracted_beam_intersection<FlatPanelDetector>::operator());
    
    vec3<double> (diffracted_beam_intersection<MultiFlatPanelDetector>::*func)(vec3<double>) const =
      &diffracted_beam_intersection<MultiFlatPanelDetector>::operator();
    
    vec3<double> (diffracted_beam_intersection<MultiFlatPanelDetector>::*func2)(vec3<double>, std::size_t) const =
      &diffracted_beam_intersection<MultiFlatPanelDetector>::operator();    
    
    class_<diffracted_beam_intersection<MultiFlatPanelDetector> >(
        "MultiFlatPanelDetector_diffracted_beam_intersection", no_init)
    .def("__call__", func)
    .def("__call__", func2);
    
    def("diffracted_beam_intersection",
      &make_diffracted_beam_intersection<FlatPanelDetector>);      
    def("diffracted_beam_intersection",
      &make_diffracted_beam_intersection<MultiFlatPanelDetector>);
  }

  void export_detector_helpers()
  {
    export_is_coordinate_valid();
    export_diffracted_beam_intersection();
  }
  

}}} // namespace = dials::model::boost_python
