/*
 * subtract_background.cc
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
#include <dials/algorithms/integration/subtract_background.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_subtract_background()
  {
    void (SubtractBackground::*subtract_single)(Reflection&) =
      &SubtractBackground::operator();

    flex_bool (SubtractBackground::*subtract_array)(ReflectionList&) = 
      &SubtractBackground::operator();

    class_ <SubtractBackground> ("SubtractBackground", no_init)
        .def(init <int,
                   double> ((
            arg("min_pixels") = 10,
            arg("n_sigma") = -1)))
        .def("__call__", subtract_single, (
            arg("reflection")))
        .def("__call__", subtract_array, (
            arg("reflections")));
  }

}}} // namespace = dials::algorithms::boost_python
