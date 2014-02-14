/*
 * transformed_shoebox.cc
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
#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/transformed_shoebox.h>
#include <dials/config.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::flex;
  using scitbx::af::flex_int;

  /** Set the data array as a flex array */
  static
  void set_data(TransformedShoebox &obj,
      typename flex<double>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.data = af::versa<double, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor()));
  }

  /** Set the bgrd array as a flex array */
  static
  void set_background(TransformedShoebox &obj,
      typename flex<double>::type &background) {
    DIALS_ASSERT(background.accessor().all().size() == 3);
    obj.background = af::versa<double, af::c_grid<3> >(
      background.handle(), af::c_grid<3>(background.accessor()));
  }

  class_<TransformedShoebox> transformed_shoebox_wrapper(const char *name)
  {
    class_<TransformedShoebox> shoebox(name);
    shoebox
      .def(init<const TransformedShoebox&>())
      .add_property("data",
        make_getter(&TransformedShoebox::data,
          return_value_policy<return_by_value>()),
        &set_data)
      .add_property("background",
        make_getter(&TransformedShoebox::background,
          return_value_policy<return_by_value>()),
        &set_background)
      .def("__eq__", &TransformedShoebox::operator==)
      .def("__ne__", &TransformedShoebox::operator!=);

    return shoebox;
  }

  void export_transformed_shoebox()
  {
    transformed_shoebox_wrapper("TransformedShoebox");
  }

}}} // namespace dials::model::boost_python
