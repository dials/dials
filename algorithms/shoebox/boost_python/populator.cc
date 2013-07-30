/*
 * populator.cc
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
#include <dials/algorithms/shoebox/populator.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_populator()
  {
    class_ <Populator> ("Populator", no_init)
      .def(init<ReflectionList&, const flex_double&, const flex_double&>((
        arg("reflection_list"),
        arg("gain_map"),
        arg("dark_map"))))
      .def("add_image", &Populator::add_image, (arg("image"), arg("index")))
      .def("allocate", &Populator::allocate)
      .def("deallocate", &Populator::deallocate)
      .def("initialize", &Populator::initialize);
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
