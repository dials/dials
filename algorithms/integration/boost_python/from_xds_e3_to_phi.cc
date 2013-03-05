/*
 * from_xds_e3_to_phi.h
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
#include <dials/algorithms/integration/from_xds_e3_to_phi.h>

namespace dials { namespace algorithms { namespace boost_python { 

  using namespace boost::python;
      
  void export_from_xds_e3_to_phi() 
  {
    class_ <FromXdsE3ToPhi> ("FromXdsE3ToPhi", no_init)
      .def(init <double,
                 double> ((
          arg("zeta"), 
          arg("phi"))))
      .def("__call__", 
        &FromXdsE3ToPhi::operator(), (
          arg("e3")));
  }

}}} // namespace = dials::algorithms::boost_python

