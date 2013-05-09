/*
 * mean_subtractor.cc
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
#include <dials/algorithms/background/mean_subtractor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_mean_subtractor()
  {
    class_<MeanSubtractor, bases<SubtractorStrategy> >(
        "MeanSubtractor", no_init)
      .def(init<std::size_t>((
        arg("min_pixels") = 10)));
  }

}}} // namespace = dials::algorithms::boost_python
