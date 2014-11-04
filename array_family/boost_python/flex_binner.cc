/*
 * flex_binner.cc
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
#include <dials/array_family/binner.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;

  void export_flex_binner() {

    class_<BinIndexer>("BinIndexer", no_init)
      .def("count", &BinIndexer::count)
      .def("sum", &BinIndexer::sum)
      .def("mean", &BinIndexer::mean)
      ;

    class_<Binner>("Binner", no_init)
      .def(init<const af::const_ref<double>&>())
      .def("bins", &Binner::bins)
      .def("indexer", &Binner::indexer)
      ;

  }

}}} // namespace dials::af::boost_python

