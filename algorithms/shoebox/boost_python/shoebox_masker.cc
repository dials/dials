/*
 * shoebox_masker.cc
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
#include <dials/algorithms/shoebox/shoebox_masker.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_shoebox_masker()
  {
    class_ <ShoeboxMasker> ("ShoeboxMasker", no_init)
      .def(init<const flex_int&>((
          arg("detector_mask"))))
      .def("__call__", 
        &ShoeboxMasker::operator(), (
          arg("reflection_list"), 
          arg("adjacency_list")));
  }

}}} // namespace = dials::algorithms::boost_python
