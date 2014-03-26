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

  boost::python::tuple block_list_next_single(
      BlockList &self,
      const af::const_ref<int, af::c_grid<2> > &image) {
    af::shared< af::const_ref<int, af::c_grid<2> > > data(1);
    data[0] = image;
    BlockList::Block block = self.next(data.const_ref());
    return boost::python::make_tuple(
        block.zrange,
        block.index,
        block.shoebox);
  }

  boost::python::tuple block_list_next_many(
      BlockList &self,
      boost::python::tuple image) {
    af::shared< af::const_ref<int, af::c_grid<2> > > data(len(image));
    for (std::size_t i = 0; i < data.size(); ++i) {
      data[i] = extract<af::const_ref< int, af::c_grid<2> > >(image[i]);
    }
    BlockList::Block block = self.next(data.const_ref());
    return boost::python::make_tuple(
        block.zrange,
        block.index,
        block.shoebox);
  }

  void export_block_list()
  {
    class_<BlockList>("BlockList", no_init)
      .def(init<
          const af::const_ref<std::size_t>&,
          const af::const_ref<int6>&,
          int2>())
      .def("next", &block_list_next_single)
      .def("next", &block_list_next_many)
      .def("empty", &BlockList::empty)
      .def("z", &BlockList::z)
      ;
  }

}}} // namespace dials::model::boost_python
