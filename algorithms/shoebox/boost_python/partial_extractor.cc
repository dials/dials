/*
 * partial_extractor.cc
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
#include <dials/algorithms/shoebox/partial_extractor.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  class_<PartialExtractor> extractor_wrapper(const char *name)
  {
    return class_ <PartialExtractor>(name, no_init)
      .def(init<const af::const_ref<std::size_t>&,
                const af::const_ref<int6>&,
                int2, std::size_t>((
          arg("panels"),
          arg("bboxes"),
          arg("zrange"),
          arg("npanels"))))
      .def(init<const af::const_ref<int6>&,
                int2>((
          arg("bboxes"),
          arg("zrange"))))
      .def("add_image", &PartialExtractor::add_image, (
        arg("panel"), 
        arg("frame"), 
        arg("image")))
      .def("indices", &PartialExtractor::indices)
      .def("shoeboxes", &PartialExtractor::shoeboxes)
      .def("shoebox_indices", &PartialExtractor::shoebox_indices);
  }
  
  void export_partial_extractor() {
    extractor_wrapper("PartialExtractor");
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
