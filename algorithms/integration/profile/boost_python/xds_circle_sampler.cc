/*
 * xds_circle_sampler.cc
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
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_xds_circle_sampler()
  {
    class_<XdsCircleSampler>("XdsCircleSampler", no_init)
      .def(init<int2>((
        arg("image_size"))))
      .def("image_size", &XdsCircleSampler::image_size)
      .def("image_centre", &XdsCircleSampler::image_centre)
      .def("r0", &XdsCircleSampler::r0)
      .def("r1", &XdsCircleSampler::r1)
      .def("r2", &XdsCircleSampler::r2)
      .def("nearest", &XdsCircleSampler::nearest)
      .def("size", &XdsCircleSampler::size)
      .def("__getitem__", &XdsCircleSampler::operator[])
      .def("__len__", &XdsCircleSampler::size);
  }

}}} // namespace = dials::algorithms::boost_python
