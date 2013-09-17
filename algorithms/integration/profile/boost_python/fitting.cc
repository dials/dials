/*
 * fitting.cc
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
#include <boost/python/iterator.hpp>
#include <dials/algorithms/integration/profile/fitting.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_fitting()
  {
    class_<ProfileFitting>("ProfileFitting", no_init)
      .def(init<const af::const_ref<double, af::c_grid<3> >&,
                const af::const_ref<double, af::c_grid<3> >&,
                const af::const_ref<double, af::c_grid<3> >&,
                double,
                std::size_t>((
        arg("profile"),
        arg("contents"),
        arg("background"),
        arg("bits") = 1e-3,
        arg("max_iter") = 10)))
      .def("intensity", &ProfileFitting::intensity)
      .def("variance", &ProfileFitting::variance)
      .def("niter", &ProfileFitting::niter)
      .def("error", &ProfileFitting::error);
  }

}}} // namespace = dials::algorithms::boost_python
