/*
 * reference_learner.cc
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
#include <dials/algorithms/integration/profile/reference_learner.h>
#include <dials/algorithms/integration/profile/grid_sampler.h>
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename ImageSampler>
  void reference_learner_wrapper(const char *name) {

    typedef ImageSampler sampler_type;
    typedef ReferenceLearner<ImageSampler> learner_type;

    class_<learner_type>(name, no_init)
      .def(init<const sampler_type&, int3, double>((
        arg("sampler"),
        arg("grid_size"),
        arg("threshold"))))
      .def("learn", &learner_type::learn)
      .def("counts", &learner_type::counts)
      .def("locate", &learner_type::locate);
  }

  void export_reference_learner()
  {
    reference_learner_wrapper<GridSampler>("GridReferenceLearner");
    reference_learner_wrapper<XdsCircleSampler>("XdsCircleReferenceLearner");
  }

}}} // namespace = dials::algorithms::boost_python
