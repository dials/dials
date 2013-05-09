/*
 * gaussian_corrected_mean_subtractor.cc
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
#include <dials/algorithms/background/gaussian_corrected_mean_subtractor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_gaussian_corrected_mean_subtractor()
  {
    class_<GaussianCorrectedMeanSubtractor, bases<MeanSubtractor> >(
        "GaussianCorrectedMeanSubtractor", no_init)
      .def(init<>());
  }

}}} // namespace = dials::algorithms::boost_python
