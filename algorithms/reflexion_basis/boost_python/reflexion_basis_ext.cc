/*
 * reflexion_basis_ext.cc
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
#include <dials/algorithms/reflexion_basis/coordinate_system.h>

namespace dials { namespace algorithms { namespace reflexion_basis {
  namespace boost_python {

  using namespace boost::python;

  void export_coordinate_system() 
  {
    // Export zeta factor functions       
    def("zeta_factor", 
      (double (*)(vec3<double>, vec3<double>, vec3<double>))&zeta_factor, 
      (arg("m2"), arg("s0"), arg("s1")));
    def("zeta_factor", 
      (double (*)(vec3<double>, vec3<double>))&zeta_factor, 
      (arg("m2"), arg("e1")));
  
    // Export coordinate system
    class_<CoordinateSystem>("CoordinateSystem", no_init)
      .def(init<vec3<double>,
                vec3<double>,
                vec3<double>,
                double>((
        arg("m2"),
        arg("s0"),
        arg("s1"),
        arg("phi"))))
      .def("m2", &CoordinateSystem::m2)
      .def("s0", &CoordinateSystem::s0)
      .def("s1", &CoordinateSystem::s1)
      .def("phi", &CoordinateSystem::phi)
      .def("p_star", &CoordinateSystem::p_star)
      .def("e1_axis", &CoordinateSystem::e1_axis)
      .def("e2_axis", &CoordinateSystem::e2_axis)
      .def("e3_axis", &CoordinateSystem::e3_axis)
      .def("zeta", &CoordinateSystem::zeta)
      .def("path_length_increase", &CoordinateSystem::path_length_increase)
      .def("limits", &CoordinateSystem::limits);


    // Export From BeamVector/Rotation angle
    class_<FromBeamVector>("FromBeamVector", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", &FromBeamVector::operator());
      
    class_<FromRotationAngleFast>("FromRotationAngleFast", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", &FromRotationAngleFast::operator());
      
    class_<FromRotationAngleAccurate>("FromRotationAngleAccurate", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", &FromRotationAngleAccurate::operator());
      
    class_<FromBeamVectorAndRotationAngle>(
        "FromBeamVectorAndRotationAngle", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", 
        (double(FromBeamVectorAndRotationAngle::*)(double)const)
          &FromBeamVectorAndRotationAngle::operator())
      .def("__call__", 
        (vec2<double>(FromBeamVectorAndRotationAngle::*)(vec3<double>)const)
          &FromBeamVectorAndRotationAngle::operator())      
      .def("__call__", 
        (vec3<double>(FromBeamVectorAndRotationAngle::*)(
        vec3<double>, double)const)
          &FromBeamVectorAndRotationAngle::operator());
      
    // Export to beamvector/rotation angle
    class_<ToBeamVector>("ToBeamVector", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", &ToBeamVector::operator());
      
    class_<ToRotationAngleFast>("ToRotationAngleFast", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", &ToRotationAngleFast::operator());
      
    class_<ToRotationAngleAccurate>("ToRotationAngleAccurate", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", &ToRotationAngleAccurate::operator());
      
    class_<ToBeamVectorAndRotationAngle>(
        "ToBeamVectorAndRotationAngle", no_init)
      .def(init<const CoordinateSystem&>())
      .def("__call__", 
        (double(ToBeamVectorAndRotationAngle::*)(double)const)
          &ToBeamVectorAndRotationAngle::operator())
      .def("__call__", 
        (vec3<double>(ToBeamVectorAndRotationAngle::*)(vec2<double>)const)
          &ToBeamVectorAndRotationAngle::operator())  
      .def("__call__", 
        (std::pair<vec3<double>, double>(ToBeamVectorAndRotationAngle::*)(
          vec3<double>)const)&ToBeamVectorAndRotationAngle::operator()); 
        
    // Export transforms from/to detector plane
    class_<FromDetector>("FromDetector", no_init)
      .def(init<const CoordinateSystem&, mat3<double> >())
      .def("__call__", 
        (vec2<double>(FromDetector::*)(vec2<double>)const)
          &FromDetector::operator())
      .def("__call__", 
        (vec3<double>(FromDetector::*)(vec2<double>, double)const)
          &FromDetector::operator());   
          
    class_<ToDetector>("ToDetector", no_init)
      .def(init<const CoordinateSystem&, mat3<double> >())
      .def("__call__",
        (vec2<double>(ToDetector::*)(vec2<double>)const)
          &ToDetector::operator())
      .def("__call__",
        (std::pair<vec2<double>, double>(ToDetector::*)(vec3<double>)const)
          &ToDetector::operator());   
  }

  BOOST_PYTHON_MODULE(dials_algorithms_reflexion_basis_ext)
  {
    export_coordinate_system();
  }

}}}} // namespace = dials::algorithms::reflexion_basis::boost_python
