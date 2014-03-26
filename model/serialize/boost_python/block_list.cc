/*
 * block_list.cc
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
#include <dials/model/serialize/block_list.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  af::shared<BasicShoebox> block_list_next_single(
      BlockList &self,
      const af::const_ref<int, af::c_grid<2> > &image) {
    af::shared< af::const_ref<int, af::c_grid<2> > > data(1);
    data[0] = image;
    return self.next(data.const_ref());
  }

  af::shared<BasicShoebox> block_list_next_many(
      BlockList &self,
      boost::python::tuple image) {
    af::shared< af::const_ref<int, af::c_grid<2> > > data(len(image));
    for (std::size_t i = 0; i < data.size(); ++i) {
      data[i] = extract<af::const_ref< int, af::c_grid<2> > >(image[i]);
    }
    return self.next(data.const_ref());
  }

  void export_basic_shoebox()
  {
    class_<BlockList>("BlockList", no_init)
      .def(init<
          const af::const_ref<std::size_t>&,
          const af::const_ref<int6>&,
          int>())
      .def("next", &block_list_next_single)
      .def("next", &block_list_next_many)
      ;
  }

}}} // namespace dials::model::boost_python
