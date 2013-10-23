/*
 * partial_shoebox.cc
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
#include <dials/model/data/partial_shoebox.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::flex_int;
   
  /** Set the data array as a flex array */
  static
  void set_data(PartialShoebox &obj, flex_int data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.data = af::versa<int, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor()));
  }

  static
  class_<PartialShoebox> partial_shoebox_wrapper(const char *name)
  {
    return class_<PartialShoebox>(name)
      .def(init<const PartialShoebox&>())
      .def(init<std::size_t, const int6&, const int2&>())
      .def(init<const int6&, const int2&>())
      .add_property("data", 
        make_getter(&PartialShoebox::data, 
          return_value_policy<return_by_value>()),
        &set_data)        
      .add_property("bbox", 
        make_getter(&PartialShoebox::bbox, 
          return_value_policy<return_by_value>()),
        make_setter(&PartialShoebox::bbox, 
          return_value_policy<return_by_value>())) 
      .add_property("zrange", 
        make_getter(&PartialShoebox::zrange, 
          return_value_policy<return_by_value>()),
        make_setter(&PartialShoebox::zrange, 
          return_value_policy<return_by_value>())) 
      .def_readwrite("panel", &PartialShoebox::panel)       
      .def("allocate", &PartialShoebox::allocate)
      .def("deallocate", &PartialShoebox::deallocate)
      .def("xoffset", &PartialShoebox::xoffset)
      .def("yoffset", &PartialShoebox::yoffset)
      .def("zoffset", &PartialShoebox::zoffset)
      .def("partial_zoffset", &PartialShoebox::partial_zoffset)
      .def("offset", &PartialShoebox::offset)
      .def("partial_offset", &PartialShoebox::partial_offset)
      .def("xsize", &PartialShoebox::xsize)
      .def("ysize", &PartialShoebox::ysize)
      .def("zsize", &PartialShoebox::zsize)
      .def("size", &PartialShoebox::size)
      .def("partial_zsize", &PartialShoebox::partial_zsize)
      .def("partial_size", &PartialShoebox::partial_size)
      .def("is_consistent", &PartialShoebox::is_consistent)
      .def("is_complete", &PartialShoebox::is_complete)
      .def("__eq__", &PartialShoebox::operator==)
      .def("__ne__", &PartialShoebox::operator!=);
  }
  
  void export_partial_shoebox()
  {
    partial_shoebox_wrapper("PartialShoebox");
  }

}}} // namespace dials::model::boost_python
