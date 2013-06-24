/*
 * grid_sampler.cc
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
#include <dials/algorithms/integration/profile/grid_sampler.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_grid_sampler()
  {
    class_<GridSampler>("GridSampler", no_init)
      .def(init<int2, int2>((
        arg("image_size"),
        arg("ngrid"))))
      .def("image_size", &GridSampler::image_size)
      .def("ngrid", &GridSampler::ngrid)
      .def("step", &GridSampler::step)
      .def("nearest", &GridSampler::nearest)
      .def("size", &GridSampler::size)
      .def("__getitem__", &GridSampler::operator[])
      .def("__len__", &GridSampler::size);
  }

}}} // namespace = dials::algorithms::boost_python
