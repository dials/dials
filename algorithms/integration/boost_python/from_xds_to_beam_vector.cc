/*
 * from_xds_to_beam_vector.cc
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
#include <dials/algorithms/integration/from_xds_to_beam_vector.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_from_xds_to_beam_vector() 
  {
    class_ <FromXdsToBeamVector> ("FromXdsToBeamVector", no_init)
      .def(init <XdsCoordinateSystem, 
                 vec3 <double> > ((
          arg("xcs"), 
          arg("s1"))))
      .def("apply", 
        &FromXdsToBeamVector::operator(), (
          arg("c")));
  }

}}} // namespace = dials::algorithms::boost_python
