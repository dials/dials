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

  boost::shared_ptr<Image>
  make_from_single(af::flex_int data, af::flex_bool mask) {
    DIALS_ASSERT(data.accessor().all().size() == 2);
    DIALS_ASSERT(mask.accessor().all().size() == 2);
    return boost::make_shared<Image>(
        af::versa<int, af::c_grid<2> >(
          data.handle(),
          af::c_grid<2>(data.accessor())),
        af::versa<bool, af::c_grid<2> >(
          mask.handle(),
          af::c_grid<2>(data.accessor())));
  }

  boost::shared_ptr<Image>
  make_from_tuple(boost::python::tuple data, boost::python::tuple mask) {
    DIALS_ASSERT(len(data) == len(mask));
    af::shared<Image::int_type> d(boost::python::len(data));
    af::shared<Image::bool_type> m(boost::python::len(mask));
    for (std::size_t i = 0; i < d.size(); ++i) {
      af::flex_int dd = boost::python::extract<af::flex_int>(data[i]);
      DIALS_ASSERT(dd.accessor().all().size() == 2);
      d[i] = af::versa<int, af::c_grid<2> >(
              dd.handle(),
              af::c_grid<2>(dd.accessor()));
      af::flex_bool mm = boost::python::extract<af::flex_bool>(mask[i]);
      DIALS_ASSERT(mm.accessor().all().size() == 2);
      m[i] = af::versa<bool, af::c_grid<2> >(
              mm.handle(),
              af::c_grid<2>(mm.accessor()));

    }
    return boost::make_shared<Image>(d.const_ref(), m.const_ref());
  }

  void export_image() {

    class_<Image>("Image", no_init)
      .def("__init__", make_constructor(make_from_single))
      .def("__init__", make_constructor(make_from_tuple))
      .def("__len__", &Image::npanels)
      ;
  }

}}} // namespace dials::model::boost_python
