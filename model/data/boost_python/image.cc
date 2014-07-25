/*
 * observation.cc
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
#include <boost/make_shared.hpp>
#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/image.h>
#include <dials/error.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  boost::shared_ptr<Image> make_from_single(af::flex_int data) {
    DIALS_ASSERT(data.accessor().all().size() == 2);
    return boost::make_shared<Image>(
        af::versa<int, af::c_grid<2> >(
          data.handle(), 
          af::c_grid<2>(data.accessor())));
  }

  boost::shared_ptr<Image> make_from_tuple(boost::python::tuple data) {
    af::shared<Image::image_type> d(boost::python::len(data));
    for (std::size_t i = 0; i < d.size(); ++i) {
      af::flex_int dd = boost::python::extract<af::flex_int>(data[i]);
      DIALS_ASSERT(dd.accessor().all().size() == 2);
      d[i] = af::versa<int, af::c_grid<2> >(
              dd.handle(), 
              af::c_grid<2>(dd.accessor()));

    }
    return boost::make_shared<Image>(d.const_ref());
  }

  void export_image() {
    
    class_<Image>("Image", no_init)
      .def("__init__", make_constructor(make_from_single))
      .def("__init__", make_constructor(make_from_tuple))
      .def("__len__", &Image::npanels)
      ;
  }

}}} // namespace dials::model::boost_python

