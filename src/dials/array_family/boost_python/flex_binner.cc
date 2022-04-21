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

  af::shared<double> sum_double(const BinIndexer &self,
                                const af::const_ref<double> &data) {
    return self.sum(data);
  }

  af::shared<int> sum_int(const BinIndexer &self, const af::const_ref<int> &data) {
    return self.sum(data);
  }

  af::shared<int> sum_bool(const BinIndexer &self, const af::const_ref<bool> &data) {
    return self.sum(data);
  }

  void export_flex_binner() {
    class_<BinIndexer>("BinIndexer", no_init)
      .def("indices", &BinIndexer::indices)
      .def("count", &BinIndexer::count)
      .def("sum", &sum_double)
      .def("sum", &sum_int)
      .def("sum", &sum_bool)
      .def("mean", &BinIndexer::mean);

    class_<Binner>("Binner", no_init)
      .def(init<const af::const_ref<double> &>())
      .def("bins", &Binner::bins)
      .def("indexer", &Binner::indexer)
      .def("__len__", &Binner::size);
  }

}}}  // namespace dials::af::boost_python
