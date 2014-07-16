/*
 * reference_learner2d.cc
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
#include <dials/algorithms/integration/profile/reference_learner2d.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;


  void export_reference_learner_2d()
  {
    class_<ReferenceLearner2D>("ReferenceLearner2D", no_init)
      .def(init<
          const GridSampler2D&,
          const af::const_ref<std::size_t>&,
          const af::const_ref<std::size_t>&,
          std::size_t>())
      .def("add", &ReferenceLearner2D::add)
      .def("finish", &ReferenceLearner2D::finish)
      .def("counts", &ReferenceLearner2D::counts)
      .def("profile", &ReferenceLearner2D::profile)
      .def("mask", &ReferenceLearner2D::mask);
  }

}}} // namespace = dials::algorithms::boost_python
