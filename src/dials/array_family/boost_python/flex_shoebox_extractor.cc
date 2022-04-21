/*
 * flex_shoebox_extractor.cc
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
#include <dials/array_family/shoebox_extractor.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;

  void export_flex_shoebox_extractor() {
    class_<ShoeboxExtractor>("ShoeboxExtractor", no_init)
      .def(init<af::reflection_table, std::size_t, int, int>(
        (boost::python::arg("data"),
         boost::python::arg("npanels"),
         boost::python::arg("frame0"),
         boost::python::arg("frame1"))))
      .def("next", &ShoeboxExtractor::next<int>)
      .def("next", &ShoeboxExtractor::next<float>)
      .def("next", &ShoeboxExtractor::next<double>)
      .def("finished", &ShoeboxExtractor::finished)
      .def("frame0", &ShoeboxExtractor::frame0)
      .def("frame1", &ShoeboxExtractor::frame1)
      .def("frame", &ShoeboxExtractor::frame)
      .def("nframes", &ShoeboxExtractor::nframes)
      .def("npanals", &ShoeboxExtractor::npanels);
  }

}}}  // namespace dials::af::boost_python
